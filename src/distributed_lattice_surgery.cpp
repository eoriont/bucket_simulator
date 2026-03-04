#include "distributed_lattice_surgery.hpp"
#include <stim/circuit/gate_target.h>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <set>
#include <unordered_set>

// Helper to create record target for DETECTOR/OBSERVABLE_INCLUDE
static inline uint32_t drec(int32_t offset) {
    return static_cast<uint32_t>(-offset) | stim::TARGET_RECORD_BIT;
}

namespace bucket_sim {

DistributedLatticeSurgeryCircuit::DistributedLatticeSurgeryCircuit(const Config& config)
    : config_(config),
      distance_(config.code_distance),
      merge_rounds_(config.merge_rounds > 0 ? config.merge_rounds : config.code_distance),
      merge_type_(config.merge_type)
{
    if (merge_type_ != MergeType::XX_MERGE_DISTRIBUTED && merge_type_ != MergeType::XX_MERGE && merge_type_ != MergeType::ZZ_MERGE) {
        throw std::invalid_argument("DistributedLatticeSurgeryCircuit requires a merge type (xx, zz, or distributed_xx)");
    }
    if (distance_ < 3 || distance_ % 2 == 0) {
        throw std::invalid_argument("Code distance must be odd and >= 3");
    }
    if (!config_.superstab_ys.empty() && merge_type_ != MergeType::XX_MERGE_DISTRIBUTED) {
        throw std::invalid_argument("superstab_ys is only supported for XX_MERGE_DISTRIBUTED");
    }

    initialize_layout();
    build_stabilizers();
}

uint32_t DistributedLatticeSurgeryCircuit::add_qubit(DQubitType type, DPatchID patch, double x, double y) {
    uint32_t idx = static_cast<uint32_t>(qubits_.size());
    DQubitInfo info{idx, type, patch, x, y};
    qubits_.push_back(info);
    return idx;
}

void DistributedLatticeSurgeryCircuit::initialize_layout() {
    // Distributed lattice surgery layout for XX merge of two d×d rotated patches.
    //
    // Key differences from standard lattice surgery:
    //   1) NO new data qubits in the merge region
    //   2) Only ONE column of merge ancillas bridging the two patches
    //   3) Patch B's checkerboard is FLIPPED so stabilizer types complement
    //      at the boundary, enabling valid weight-4 merge stabilizers
    //
    // Convention (same as original):
    //   INTEGER coords  = ANCILLA qubits (measured each round)
    //   HALF-INT coords = DATA qubits (persist)
    //
    // Layout (d=5 example):
    //   Patch A data: x = 0.5, 1.5, 2.5, 3.5, 4.5  (5 columns)
    //   Patch A anc:  x = 0, 1, 2, 3, 4             (corner at 0 + 4 interior)
    //   Merge anc:    x = 5                          (1 column, single row of ancillas)
    //   Patch B data: x = 5.5, 6.5, 7.5, 8.5, 9.5   (5 columns)
    //   Patch B anc:  x = 6, 7, 8, 9, 10             (4 interior + corner at 10)
    //
    // Merge ancilla at (5, y) touches:
    //   - Patch A data at (4.5, y±0.5) — local CNOTs
    //   - Patch B data at (5.5, y±0.5) — remote CNOTs over interconnect

    uint32_t d = distance_;

    // ========== PATCH A (left side) ==========
    // Standard rotated surface code

    // Left boundary X ancillas at x=0: y=2,4,...,d-1
    for (uint32_t y = 2; y <= d - 1; y += 2) {
        add_qubit(DQubitType::X_ANCILLA, DPatchID::PATCH_A, 0, y);
    }

    for (uint32_t ix = 0; ix < d - 1; ix++) {
        double x_data = ix + 0.5;
        int x_anc = ix + 1;

        // Data qubits
        for (uint32_t j = 0; j < d; j++) {
            add_qubit(DQubitType::DATA, DPatchID::PATCH_A, x_data, j + 0.5);
        }

        // Ancilla qubits
        uint32_t y_min = (x_anc % 2 == 0) ? 1 : 0;
        for (uint32_t j = 0; j < d; j++) {
            uint32_t y = y_min + j;
            bool is_x = ((y - y_min) % 2 == 1);
            add_qubit(is_x ? DQubitType::X_ANCILLA : DQubitType::Z_ANCILLA,
                     DPatchID::PATCH_A, static_cast<double>(x_anc), y);
        }
    }

    // Last column of data for Patch A
    for (uint32_t j = 0; j < d; j++) {
        add_qubit(DQubitType::DATA, DPatchID::PATCH_A, d - 0.5, j + 0.5);
    }

    // ========== MERGE ANCILLAS (single column at x = d) ==========
    // These ancillas bridge the boundary. Each merge ancilla at (d, y) connects
    // to data qubits at (d-0.5, y±0.5) from Patch A and (d+0.5, y±0.5) from Patch B.
    //
    // The merge column simply continues the global checkerboard pattern.
    // Same-type neighboring stabilizers always commute (XX=I, ZZ=I on shared qubits),
    // and different-type neighbors share 0 or 2 data qubits (also commute).
    // So NO flip of Patch B is needed — both patches use the same orientation.

    int merge_x = static_cast<int>(d);
    uint32_t merge_y_min = (merge_x % 2 == 0) ? 1 : 0;
    for (uint32_t j = 0; j < d; j++) {
        uint32_t y = merge_y_min + j;
        bool is_x = ((y - merge_y_min) % 2 == 1);
        add_qubit(is_x ? DQubitType::X_ANCILLA : DQubitType::Z_ANCILLA,
                 DPatchID::MERGE, static_cast<double>(merge_x), y);
    }

    // ========== SEAM ANCILLAS (at x = d, same coords as merge X positions) ==========
    // Seam A: X ancillas owned by Patch A, couple only to Patch A data at (d-0.5, y±0.5)
    // Seam B: X ancillas owned by Patch B, couple only to Patch B data at (d+0.5, y±0.5)
    // y = 1, 3, ..., d-2 (the X-type positions in the merge column)
    for (uint32_t y = 1; y <= d - 2; y += 2) {
        add_qubit(DQubitType::X_ANCILLA, DPatchID::SEAM_A, static_cast<double>(d), y);
    }
    for (uint32_t y = 1; y <= d - 2; y += 2) {
        add_qubit(DQubitType::X_ANCILLA, DPatchID::SEAM_B, static_cast<double>(d), y);
    }

    // ========== PATCH B (right side, SAME checkerboard as Patch A) ==========
    // Patch B's first data column is at x = d + 0.5 (directly adjacent to merge).
    // The checkerboard follows the same global pattern as Patch A — no flip needed.
    // This works because same-type adjacent stabilizers always commute,
    // and different-type neighbors share an even number of data qubits.

    double patch_b_start = d + 0.5;

    for (uint32_t ix = 0; ix < d - 1; ix++) {
        double x_data = patch_b_start + ix;
        int x_anc_int = static_cast<int>(d + 1 + ix);
        double x_anc = static_cast<double>(x_anc_int);

        // Data qubits
        for (uint32_t j = 0; j < d; j++) {
            add_qubit(DQubitType::DATA, DPatchID::PATCH_B, x_data, j + 0.5);
        }

        // Ancilla qubits — SAME pattern as Patch A (no flip)
        uint32_t y_min = (x_anc_int % 2 == 0) ? 1 : 0;
        for (uint32_t j = 0; j < d; j++) {
            uint32_t y = y_min + j;
            bool is_x = ((y - y_min) % 2 == 1);
            add_qubit(is_x ? DQubitType::X_ANCILLA : DQubitType::Z_ANCILLA,
                     DPatchID::PATCH_B, x_anc, y);
        }
    }

    // Last data column for Patch B
    for (uint32_t j = 0; j < d; j++) {
        add_qubit(DQubitType::DATA, DPatchID::PATCH_B, patch_b_start + d - 1, j + 0.5);
    }

    // Right boundary X ancillas for Patch B at x=2d: y=2,4,...,d-1
    int right_x = static_cast<int>(d + d);
    for (uint32_t y = 2; y <= d - 1; y += 2) {
        add_qubit(DQubitType::X_ANCILLA, DPatchID::PATCH_B, static_cast<double>(right_x), y);
    }

    std::cerr << "Distributed lattice surgery layout:" << std::endl;
    std::cerr << "  Distance: " << d << std::endl;
    std::cerr << "  Total qubits: " << qubits_.size() << std::endl;

    size_t n_data = 0, n_anc = 0, n_merge = 0, n_seam = 0;
    for (const auto& q : qubits_) {
        if (q.type == DQubitType::DATA) n_data++;
        else if (q.patch == DPatchID::MERGE) n_merge++;
        else if (q.patch == DPatchID::SEAM_A || q.patch == DPatchID::SEAM_B) n_seam++;
        else n_anc++;
    }
    std::cerr << "  Data qubits: " << n_data << " (2 × " << d*d << " = " << 2*d*d << ")" << std::endl;
    std::cerr << "  Patch ancillas: " << n_anc << std::endl;
    std::cerr << "  Seam ancillas: " << n_seam << " (" << n_seam/2 << " per patch)" << std::endl;
    std::cerr << "  Merge ancillas: " << n_merge << " (single column)" << std::endl;
    std::cerr << "  No merge data qubits (remote CNOTs only)" << std::endl;
}

void DistributedLatticeSurgeryCircuit::build_stabilizers() {
    // Build stabilizers by finding adjacent data qubits for each ancilla.
    // Each ancilla at (x, y) looks for data at (x±0.5, y±0.5).
    //
    // Patch filtering:
    //   SEAM_A only sees PATCH_A data
    //   SEAM_B only sees PATCH_B data
    //   MERGE sees both patches
    //   PATCH_A/PATCH_B see their own data (naturally, by coordinate)

    auto find_data_qubit = [this](double x, double y, DPatchID owner) -> int32_t {
        for (const auto& q : qubits_) {
            if (q.type != DQubitType::DATA) continue;
            if (std::abs(q.x - x) > 0.1 || std::abs(q.y - y) > 0.1) continue;
            // Patch filter for seam ancillas
            if (owner == DPatchID::SEAM_A && q.patch != DPatchID::PATCH_A) continue;
            if (owner == DPatchID::SEAM_B && q.patch != DPatchID::PATCH_B) continue;
            return static_cast<int32_t>(q.index);
        }
        return -1;
    };

    for (const auto& q : qubits_) {
        if (q.type == DQubitType::DATA) continue;

        DStabilizer stab;
        stab.is_x_type = (q.type == DQubitType::X_ANCILLA);
        stab.ancilla = q.index;
        stab.patch = q.patch;
        stab.crosses_boundary = false;

        std::vector<std::pair<double, double>> offsets = {
            {-0.5, -0.5}, {0.5, -0.5}, {-0.5, 0.5}, {0.5, 0.5}
        };

        bool has_a = false, has_b = false;
        for (const auto& [dx, dy] : offsets) {
            int32_t data_idx = find_data_qubit(q.x + dx, q.y + dy, q.patch);
            if (data_idx >= 0) {
                stab.data_qubits.push_back(static_cast<uint32_t>(data_idx));
                if (qubits_[data_idx].patch == DPatchID::PATCH_A) has_a = true;
                if (qubits_[data_idx].patch == DPatchID::PATCH_B) has_b = true;
            }
        }
        stab.crosses_boundary = (has_a && has_b);

        if (!stab.data_qubits.empty()) {
            if (q.patch == DPatchID::PATCH_A) {
                patch_a_stabilizers_.push_back(stab);
            } else if (q.patch == DPatchID::PATCH_B) {
                patch_b_stabilizers_.push_back(stab);
            } else if (q.patch == DPatchID::SEAM_A) {
                seam_a_stabilizers_.push_back(stab);
            } else if (q.patch == DPatchID::SEAM_B) {
                seam_b_stabilizers_.push_back(stab);
            } else {
                merge_stabilizers_.push_back(stab);
            }
        }
    }

    std::cerr << "  Patch A stabilizers: " << patch_a_stabilizers_.size() << std::endl;
    std::cerr << "  Patch B stabilizers: " << patch_b_stabilizers_.size() << std::endl;
    std::cerr << "  Seam A stabilizers: " << seam_a_stabilizers_.size() << std::endl;
    std::cerr << "  Seam B stabilizers: " << seam_b_stabilizers_.size() << std::endl;
    std::cerr << "  Merge stabilizers: " << merge_stabilizers_.size() << std::endl;

    size_t cross_count = 0;
    for (const auto& s : merge_stabilizers_) {
        if (s.crosses_boundary) cross_count++;
        std::cerr << "    Merge " << (s.is_x_type ? "X" : "Z") << " ancilla "
                  << s.ancilla << " at (" << qubits_[s.ancilla].x << ", " << qubits_[s.ancilla].y
                  << ") weight=" << s.data_qubits.size()
                  << (s.crosses_boundary ? " [CROSS-BOUNDARY]" : "") << std::endl;
    }
    std::cerr << "  Cross-boundary stabilizers: " << cross_count << std::endl;
}

void DistributedLatticeSurgeryCircuit::generate_circuit() {
    circuit_ = stim::Circuit();

    // Add qubit coordinates
    // Seam ancillas share (x=d, y) with merge ancillas — offset slightly for Crumble display
    for (const auto& q : qubits_) {
        double display_x = q.x;
        if (q.patch == DPatchID::SEAM_A) display_x -= 0.125;
        else if (q.patch == DPatchID::SEAM_B) display_x += 0.125;
        circuit_.safe_append_u("QUBIT_COORDS", {q.index}, {display_x, q.y});
    }

    generate_general_circuit();
}

void DistributedLatticeSurgeryCircuit::generate_general_circuit() {
    // General distributed lattice surgery circuit for any odd distance >= 3.
    //
    // Structure:
    //   1. Initialize all data qubits
    //   2. Pre-merge stabilizer rounds (interior patches + seam boundary X, no merge)
    //   3. Merge rounds (interior patches + merge, no seam)
    //   4. Post-merge stabilizer round (interior patches + seam boundary X, no merge)
    //   5. Final data qubit measurements
    //
    // Ancilla groups:
    //   z_ancillas / x_ancillas: interior patch Z/X (all rounds)
    //   seam_a_x_ancillas / seam_b_x_ancillas: seam X (pre/post-merge only)
    //   merge_z_ancillas / merge_x_ancillas: merge (merge rounds only)

    // Collect qubits by type and role
    std::vector<uint32_t> all_data;
    std::vector<uint32_t> patch_a_data;
    std::vector<uint32_t> patch_b_data;
    std::vector<uint32_t> z_ancillas;           // Interior patch Z (all rounds)
    std::vector<uint32_t> x_ancillas;           // Interior patch X (all rounds)
    std::vector<uint32_t> seam_a_x_ancillas;    // Seam X owned by Patch A (pre/post-merge)
    std::vector<uint32_t> seam_b_x_ancillas;    // Seam X owned by Patch B (pre/post-merge)
    std::vector<uint32_t> merge_z_ancillas;     // Merge Z (merge rounds only)
    std::vector<uint32_t> merge_x_ancillas;     // Merge X (merge rounds only)

    for (const auto& q : qubits_) {
        if (q.type == DQubitType::DATA) {
            all_data.push_back(q.index);
            if (q.patch == DPatchID::PATCH_A) patch_a_data.push_back(q.index);
            else patch_b_data.push_back(q.index);
        } else if (q.type == DQubitType::Z_ANCILLA) {
            if (q.patch == DPatchID::MERGE) merge_z_ancillas.push_back(q.index);
            else z_ancillas.push_back(q.index);
        } else if (q.type == DQubitType::X_ANCILLA) {
            if (q.patch == DPatchID::MERGE) merge_x_ancillas.push_back(q.index);
            else if (q.patch == DPatchID::SEAM_A) seam_a_x_ancillas.push_back(q.index);
            else if (q.patch == DPatchID::SEAM_B) seam_b_x_ancillas.push_back(q.index);
            else x_ancillas.push_back(q.index);
        }
    }

    // Build CNOT layers with patch-aware data lookup
    // Three groups: cx_layer (interior, all rounds), seam_cx (pre/post), merge_cx (merge)
    //
    // Each ancilla connects to data qubits in 4 directions:
    //   Layer 1: NE (+0.5, -0.5)
    //   Layer 2: SE (+0.5, +0.5)
    //   Layer 3: NW (-0.5, -0.5)
    //   Layer 4: SW (-0.5, +0.5)

    // === Superstabilizer classification ===
    // For each X-type merge-column position listed in config_.superstab_ys, replace the remote
    // merge-X measurement with local seam_a[j] × seam_b[j] (the superstabilizer).
    // seam_a_x_ancillas[j] and merge_x_ancillas[j] share the same y-coordinate (both ordered
    // y=1,3,...,d-2), so a direct index correspondence holds.
    std::unordered_set<uint32_t> ss_y_set(
        config_.superstab_ys.begin(), config_.superstab_ys.end());

    std::vector<size_t> ss_seam_indices;      // indices into seam_a/b arrays (superstab positions)
    std::vector<size_t> normal_seam_indices;  // indices into seam_a/b arrays (non-superstab positions)
    std::vector<uint32_t> ss_seam_a, ss_seam_b;
    std::vector<uint32_t> normal_merge_x;
    std::unordered_set<uint32_t> ss_merge_x_anc_set; // merge-X ancilla indices to suppress in merge CX
    std::unordered_set<uint32_t> ss_seam_anc_set;    // seam ancilla indices to activate in merge CX

    for (size_t j = 0; j < seam_a_x_ancillas.size(); j++) {
        uint32_t y = static_cast<uint32_t>(std::round(qubits_[seam_a_x_ancillas[j]].y));
        if (ss_y_set.count(y)) {
            ss_seam_indices.push_back(j);
            ss_seam_a.push_back(seam_a_x_ancillas[j]);
            ss_seam_b.push_back(seam_b_x_ancillas[j]);
            ss_merge_x_anc_set.insert(merge_x_ancillas[j]);
            ss_seam_anc_set.insert(seam_a_x_ancillas[j]);
            ss_seam_anc_set.insert(seam_b_x_ancillas[j]);
        } else {
            normal_seam_indices.push_back(j);
            normal_merge_x.push_back(merge_x_ancillas[j]);
        }
    }

    const size_t ss_count = ss_seam_a.size();

    std::vector<uint32_t> cx_layer1, cx_layer2, cx_layer3, cx_layer4;
    std::vector<uint32_t> seam_cx1, seam_cx2, seam_cx3, seam_cx4;
    std::vector<uint32_t> merge_cx1, merge_cx2, merge_cx3, merge_cx4;

    auto find_data_at = [this](double ax, double ay, double dx, double dy, DPatchID owner) -> int32_t {
        double tx = ax + dx, ty = ay + dy;
        for (const auto& q : qubits_) {
            if (q.type != DQubitType::DATA) continue;
            if (std::abs(q.x - tx) > 0.1 || std::abs(q.y - ty) > 0.1) continue;
            if (owner == DPatchID::SEAM_A && q.patch != DPatchID::PATCH_A) continue;
            if (owner == DPatchID::SEAM_B && q.patch != DPatchID::PATCH_B) continue;
            return static_cast<int32_t>(q.index);
        }
        return -1;
    };

    for (const auto& q : qubits_) {
        if (q.type == DQubitType::DATA) continue;

        bool is_merge = (q.patch == DPatchID::MERGE);
        bool is_seam = (q.patch == DPatchID::SEAM_A || q.patch == DPatchID::SEAM_B);
        bool is_x = (q.type == DQubitType::X_ANCILLA);

        std::vector<std::pair<double, double>> dirs = {
            {0.5, -0.5},   // Layer 1: NE
            {0.5, 0.5},    // Layer 2: SE
            {-0.5, -0.5},  // Layer 3: NW
            {-0.5, 0.5}    // Layer 4: SW
        };

        for (int layer = 0; layer < 4; layer++) {
            int32_t data_idx = find_data_at(q.x, q.y, dirs[layer].first, dirs[layer].second, q.patch);
            if (data_idx < 0) continue;

            uint32_t ctrl, tgt;
            if (is_x) {
                ctrl = q.index;
                tgt = static_cast<uint32_t>(data_idx);
            } else {
                ctrl = static_cast<uint32_t>(data_idx);
                tgt = q.index;
            }

            auto& target_layer = is_merge ?
                (layer == 0 ? merge_cx1 : layer == 1 ? merge_cx2 : layer == 2 ? merge_cx3 : merge_cx4) :
                is_seam ?
                (layer == 0 ? seam_cx1 : layer == 1 ? seam_cx2 : layer == 2 ? seam_cx3 : seam_cx4) :
                (layer == 0 ? cx_layer1 : layer == 1 ? cx_layer2 : layer == 2 ? cx_layer3 : cx_layer4);

            target_layer.push_back(ctrl);
            target_layer.push_back(tgt);
        }
    }

    // === CX layer filtering for superstabilizer positions ===
    // normal_mcx*: merge-column CX pairs with ss_merge_x ancillas removed (remote ops)
    // ss_scx*:     seam CX pairs for ss positions only (local ops replacing remote merge-X)
    auto filter_cx_excl = [](const std::vector<uint32_t>& cx,
                              const std::unordered_set<uint32_t>& excl) {
        std::vector<uint32_t> out;
        for (size_t i = 0; i < cx.size(); i += 2) {
            if (!excl.count(cx[i]) && !excl.count(cx[i + 1]))
                out.push_back(cx[i]), out.push_back(cx[i + 1]);
        }
        return out;
    };
    auto filter_cx_keep = [](const std::vector<uint32_t>& cx,
                              const std::unordered_set<uint32_t>& keep) {
        std::vector<uint32_t> out;
        for (size_t i = 0; i < cx.size(); i += 2) {
            if (keep.count(cx[i]) || keep.count(cx[i + 1]))
                out.push_back(cx[i]), out.push_back(cx[i + 1]);
        }
        return out;
    };

    auto norm_mcx1 = filter_cx_excl(merge_cx1, ss_merge_x_anc_set);
    auto norm_mcx2 = filter_cx_excl(merge_cx2, ss_merge_x_anc_set);
    auto norm_mcx3 = filter_cx_excl(merge_cx3, ss_merge_x_anc_set);
    auto norm_mcx4 = filter_cx_excl(merge_cx4, ss_merge_x_anc_set);

    auto ss_scx1 = filter_cx_keep(seam_cx1, ss_seam_anc_set);
    auto ss_scx2 = filter_cx_keep(seam_cx2, ss_seam_anc_set);
    auto ss_scx3 = filter_cx_keep(seam_cx3, ss_seam_anc_set);
    auto ss_scx4 = filter_cx_keep(seam_cx4, ss_seam_anc_set);

    // Counts for offset arithmetic
    size_t patch_z_count = z_ancillas.size();
    size_t patch_x_count = x_ancillas.size();
    size_t seam_a_count = seam_a_x_ancillas.size();
    size_t seam_b_count = seam_b_x_ancillas.size();
    size_t seam_count = seam_a_count + seam_b_count;
    size_t merge_z_count = merge_z_ancillas.size();
    size_t merge_x_count = merge_x_ancillas.size();
    size_t normal_merge_x_count = normal_merge_x.size(); // = merge_x_count - ss_count

    // Pre-merge round measurement order: M(z_ancillas), MX(x_ancillas), MX(seam_a), MX(seam_b)
    size_t pre_total = patch_z_count + patch_x_count + seam_count;
    // Merge round measurement order:
    //   M(z), M(merge_z), MX(x), MX(normal_merge_x), MX(ss_seam_a), MX(ss_seam_b)
    // ss_merge_total = merge_total + ss_count (when ss_count==0 reduces to merge_total exactly)
    size_t merge_total = patch_z_count + merge_z_count + patch_x_count + merge_x_count;
    size_t ss_merge_total = merge_total + ss_count;
    // Post-merge = same as pre-merge
    size_t post_total = pre_total;

    // ====== Initial reset of all data qubits ======
    circuit_.safe_append_u("TICK", {}, {});
    circuit_.safe_append_u("R", all_data, {});

    // ====== Round 0: Pre-merge round (interior patches + seam, no merge) ======
    {
        // Combine seam ancillas for reset/measure
        std::vector<uint32_t> all_seam_x;
        all_seam_x.insert(all_seam_x.end(), seam_a_x_ancillas.begin(), seam_a_x_ancillas.end());
        all_seam_x.insert(all_seam_x.end(), seam_b_x_ancillas.begin(), seam_b_x_ancillas.end());

        std::vector<uint32_t> all_x_pre = x_ancillas;
        all_x_pre.insert(all_x_pre.end(), all_seam_x.begin(), all_seam_x.end());

        circuit_.safe_append_u("TICK", {}, {});
        circuit_.safe_append_u("R", z_ancillas, {});
        circuit_.safe_append_u("RX", all_x_pre, {});

        // CNOT layers: interior + seam
        auto emit_pre_layer = [&](const std::vector<uint32_t>& interior_cx,
                                   const std::vector<uint32_t>& seam_cx_layer) {
            std::vector<uint32_t> combined = interior_cx;
            combined.insert(combined.end(), seam_cx_layer.begin(), seam_cx_layer.end());
            circuit_.safe_append_u("TICK", {}, {});
            if (!combined.empty()) circuit_.safe_append_u("CX", combined, {});
        };

        emit_pre_layer(cx_layer1, seam_cx1);
        emit_pre_layer(cx_layer2, seam_cx2);
        emit_pre_layer(cx_layer3, seam_cx3);
        emit_pre_layer(cx_layer4, seam_cx4);

        // Measurements: M(z), MX(x), MX(seam_a), MX(seam_b)
        circuit_.safe_append_u("TICK", {}, {});
        circuit_.safe_append_u("M", z_ancillas, {});
        circuit_.safe_append_u("MX", all_x_pre, {});

        // Round 0 detectors: Z-ancillas only (first measurement after R|0⟩, deterministic)
        for (size_t i = 0; i < patch_z_count; i++) {
            const auto& q = qubits_[z_ancillas[i]];
            circuit_.safe_append_u("DETECTOR", {drec(-(int32_t)(pre_total - i))}, {q.x, q.y, 0});
        }

        circuit_.safe_append_u("TICK", {}, {});
    }

    // ====== Merge rounds ======
    // Block layout per round (ss_merge_total measurements):
    //   M(z[0..pz-1])             offsets: -(ss_merge_total - i)
    //   M(merge_z[0..mz-1])       offsets: -(ss_merge_total - pz - i)
    //   MX(x[0..px-1])            offsets: -(ss_merge_total - pz - mz - i)
    //   MX(normal_merge_x[0..nm]) offsets: -(2*ss + nm - m)
    //   MX(ss_seam_a[0..ss-1])    offsets: -(2*ss - k)
    //   MX(ss_seam_b[0..ss-1])    offsets: -(ss - k)          [last in block]
    // When ss_count==0, ss_merge_total==merge_total and all formulas reduce to the original.
    for (uint32_t r = 0; r < merge_rounds_; r++) {
        if (r < merge_rounds_ - 1) {
            circuit_.safe_append_u("TICK", {}, {});
        }

        // Reset: interior Z + merge_z; interior X + normal_merge_x + ss_seam_a + ss_seam_b
        std::vector<uint32_t> all_z_reset = z_ancillas;
        all_z_reset.insert(all_z_reset.end(), merge_z_ancillas.begin(), merge_z_ancillas.end());
        circuit_.safe_append_u("R", all_z_reset, {});

        std::vector<uint32_t> all_x_reset = x_ancillas;
        all_x_reset.insert(all_x_reset.end(), normal_merge_x.begin(), normal_merge_x.end());
        all_x_reset.insert(all_x_reset.end(), ss_seam_a.begin(), ss_seam_a.end());
        all_x_reset.insert(all_x_reset.end(), ss_seam_b.begin(), ss_seam_b.end());
        circuit_.safe_append_u("RX", all_x_reset, {});

        circuit_.safe_append_u("TICK", {}, {});

        // CNOT layers: interior + normal remote merge + local ss seam
        // Noise is injected later by inject_interconnect_noise() in the simulator.
        auto emit_ss_merge_layer = [&](const std::vector<uint32_t>& int_cx,
                                       const std::vector<uint32_t>& nm_cx,
                                       const std::vector<uint32_t>& ss_cx) {
            std::vector<uint32_t> combined = int_cx;
            combined.insert(combined.end(), nm_cx.begin(), nm_cx.end());
            combined.insert(combined.end(), ss_cx.begin(), ss_cx.end());
            if (!combined.empty()) circuit_.safe_append_u("CX", combined, {});
            circuit_.safe_append_u("TICK", {}, {});
        };

        emit_ss_merge_layer(cx_layer1, norm_mcx1, ss_scx1);
        emit_ss_merge_layer(cx_layer2, norm_mcx2, ss_scx2);
        emit_ss_merge_layer(cx_layer3, norm_mcx3, ss_scx3);
        emit_ss_merge_layer(cx_layer4, norm_mcx4, ss_scx4);

        // Measurements: M(z, merge_z), MX(x, normal_merge_x, ss_seam_a, ss_seam_b)
        std::vector<uint32_t> all_z_meas = z_ancillas;
        all_z_meas.insert(all_z_meas.end(), merge_z_ancillas.begin(), merge_z_ancillas.end());
        circuit_.safe_append_u("M", all_z_meas, {});

        std::vector<uint32_t> all_x_meas = x_ancillas;
        all_x_meas.insert(all_x_meas.end(), normal_merge_x.begin(), normal_merge_x.end());
        all_x_meas.insert(all_x_meas.end(), ss_seam_a.begin(), ss_seam_a.end());
        all_x_meas.insert(all_x_meas.end(), ss_seam_b.begin(), ss_seam_b.end());
        circuit_.safe_append_u("MX", all_x_meas, {});

        // ---- Detectors ----
        const size_t pz = patch_z_count, mz = merge_z_count, px = patch_x_count;
        const size_t nm = normal_merge_x_count, ss = ss_count;

        // Patch Z: compare with previous round
        for (size_t i = 0; i < pz; i++) {
            const auto& q = qubits_[z_ancillas[i]];
            int32_t curr = -(int32_t)(ss_merge_total - i);
            int32_t prev = (r == 0)
                ? -(int32_t)(ss_merge_total + pre_total - i)
                : -(int32_t)(2 * ss_merge_total - i);
            circuit_.safe_append_u("DETECTOR", {drec(curr), drec(prev)},
                {q.x, q.y, static_cast<double>(r + 1)});
        }

        // Merge Z: no detector on round 0 (anti-commutes with seam X → random first outcome)
        if (r >= 1) {
            for (size_t i = 0; i < mz; i++) {
                const auto& q = qubits_[merge_z_ancillas[i]];
                int32_t curr = -(int32_t)(ss_merge_total - pz - i);
                int32_t prev = -(int32_t)(2 * ss_merge_total - pz - i);
                circuit_.safe_append_u("DETECTOR", {drec(curr), drec(prev)},
                    {q.x, q.y, static_cast<double>(r + 1)});
            }
        }

        // Patch X: compare with previous round
        for (size_t i = 0; i < px; i++) {
            const auto& q = qubits_[x_ancillas[i]];
            int32_t curr = -(int32_t)(ss_merge_total - pz - mz - i);
            int32_t prev = (r == 0)
                ? -(int32_t)(ss_merge_total + pre_total - pz - i)
                : -(int32_t)(2 * ss_merge_total - pz - mz - i);
            circuit_.safe_append_u("DETECTOR", {drec(curr), drec(prev)},
                {q.x, q.y, static_cast<double>(r + 1)});
        }

        // Normal merge-X positions: 3-term transition (r==0) or 2-term steady-state (r>0)
        for (size_t m = 0; m < nm; m++) {
            const auto& q = qubits_[normal_merge_x[m]];
            int32_t curr = -(int32_t)(2 * ss + nm - m);
            size_t j = normal_seam_indices[m]; // index in seam_a/seam_b arrays
            if (r == 0) {
                // merge_x ⊕ seam_a_prev ⊕ seam_b_prev  (transition from pre-merge)
                int32_t prev_a = -(int32_t)(ss_merge_total + seam_count - j);
                int32_t prev_b = -(int32_t)(ss_merge_total + seam_b_count - j);
                circuit_.safe_append_u("DETECTOR", {drec(curr), drec(prev_a), drec(prev_b)},
                    {q.x, q.y, static_cast<double>(r + 1)});
            } else {
                int32_t prev = -(int32_t)(ss_merge_total + 2 * ss + nm - m);
                circuit_.safe_append_u("DETECTOR", {drec(curr), drec(prev)},
                    {q.x, q.y, static_cast<double>(r + 1)});
            }
        }

        // Superstabilizer positions: always 4-term seam_a ⊕ seam_b ⊕ prev_seam_a ⊕ prev_seam_b
        for (size_t k = 0; k < ss; k++) {
            const auto& q = qubits_[ss_seam_a[k]];
            int32_t curr_a = -(int32_t)(2 * ss - k);
            int32_t curr_b = -(int32_t)(ss - k);
            int32_t prev_a, prev_b;
            if (r == 0) {
                // Compare with seam_a/b from pre-merge round
                size_t j = ss_seam_indices[k];
                prev_a = -(int32_t)(ss_merge_total + seam_count - j);
                prev_b = -(int32_t)(ss_merge_total + seam_b_count - j);
            } else {
                prev_a = -(int32_t)(ss_merge_total + 2 * ss - k);
                prev_b = -(int32_t)(ss_merge_total + ss - k);
            }
            circuit_.safe_append_u("DETECTOR",
                {drec(curr_a), drec(curr_b), drec(prev_a), drec(prev_b)},
                {q.x, q.y, static_cast<double>(r + 1)});
        }

        circuit_.safe_append_u("TICK", {}, {});
    }

    // ====== Final round: Post-merge (interior patches + seam, no merge) ======
    {
        std::vector<uint32_t> all_seam_x;
        all_seam_x.insert(all_seam_x.end(), seam_a_x_ancillas.begin(), seam_a_x_ancillas.end());
        all_seam_x.insert(all_seam_x.end(), seam_b_x_ancillas.begin(), seam_b_x_ancillas.end());

        std::vector<uint32_t> all_x_post = x_ancillas;
        all_x_post.insert(all_x_post.end(), all_seam_x.begin(), all_seam_x.end());

        circuit_.safe_append_u("TICK", {}, {});
        circuit_.safe_append_u("R", z_ancillas, {});
        circuit_.safe_append_u("RX", all_x_post, {});

        auto emit_post_layer = [&](const std::vector<uint32_t>& interior_cx,
                                    const std::vector<uint32_t>& seam_cx_layer) {
            std::vector<uint32_t> combined = interior_cx;
            combined.insert(combined.end(), seam_cx_layer.begin(), seam_cx_layer.end());
            circuit_.safe_append_u("TICK", {}, {});
            if (!combined.empty()) circuit_.safe_append_u("CX", combined, {});
        };

        emit_post_layer(cx_layer1, seam_cx1);
        emit_post_layer(cx_layer2, seam_cx2);
        emit_post_layer(cx_layer3, seam_cx3);
        emit_post_layer(cx_layer4, seam_cx4);

        // Measurements: M(z), MX(x), MX(seam_a), MX(seam_b)
        circuit_.safe_append_u("TICK", {}, {});
        circuit_.safe_append_u("M", z_ancillas, {});
        circuit_.safe_append_u("MX", all_x_post, {});

        // Post-merge detectors: previous round was the last ss_merge round
        uint32_t final_round = merge_rounds_ + 1;
        const size_t pz = patch_z_count, mz = merge_z_count;
        const size_t nm = normal_merge_x_count, ss = ss_count;

        // Patch Z: compare with last ss_merge round
        for (size_t i = 0; i < patch_z_count; i++) {
            const auto& q = qubits_[z_ancillas[i]];
            int32_t curr = -(int32_t)(post_total - i);
            int32_t prev = -(int32_t)(post_total + ss_merge_total - i);
            circuit_.safe_append_u("DETECTOR", {drec(curr), drec(prev)},
                {q.x, q.y, static_cast<double>(final_round)});
        }

        // Patch X: compare with last ss_merge round
        for (size_t i = 0; i < patch_x_count; i++) {
            const auto& q = qubits_[x_ancillas[i]];
            int32_t curr = -(int32_t)(post_total - pz - i);
            // x[i] in ss_merge sits at -(ss_merge_total - pz - mz - i) from that block's end
            int32_t prev = -(int32_t)(post_total + ss_merge_total - pz - mz - i);
            circuit_.safe_append_u("DETECTOR", {drec(curr), drec(prev)},
                {q.x, q.y, static_cast<double>(final_round)});
        }

        // Seam X exit: inverse of the merge-entry transition.
        // Normal positions (j ∈ normal_seam_indices): 3-term seam_a ⊕ seam_b ⊕ last_nmx
        for (size_t m = 0; m < nm; m++) {
            size_t j = normal_seam_indices[m];
            const auto& q = qubits_[seam_a_x_ancillas[j]];
            int32_t curr_a = -(int32_t)(seam_count - j);
            int32_t curr_b = -(int32_t)(seam_b_count - j);
            // normal_merge_x[m] in ss_merge: -(2*ss + nm - m) from that block's end
            int32_t prev_nmx = -(int32_t)(post_total + 2 * ss + nm - m);
            circuit_.safe_append_u("DETECTOR",
                {drec(curr_a), drec(curr_b), drec(prev_nmx)},
                {q.x, q.y, static_cast<double>(final_round)});
        }
        // Superstabilizer positions (k ∈ ss_seam_indices): 4-term seam_a ⊕ seam_b ⊕ prev_a ⊕ prev_b
        for (size_t k = 0; k < ss; k++) {
            size_t j = ss_seam_indices[k];
            const auto& q = qubits_[ss_seam_a[k]];
            int32_t curr_a = -(int32_t)(seam_count - j);
            int32_t curr_b = -(int32_t)(seam_b_count - j);
            // ss_seam_a[k]/ss_seam_b[k] in ss_merge: -(2*ss-k) and -(ss-k) from that block's end
            int32_t prev_a = -(int32_t)(post_total + 2 * ss - k);
            int32_t prev_b = -(int32_t)(post_total + ss - k);
            circuit_.safe_append_u("DETECTOR",
                {drec(curr_a), drec(curr_b), drec(prev_a), drec(prev_b)},
                {q.x, q.y, static_cast<double>(final_round)});
        }

        circuit_.safe_append_u("TICK", {}, {});
    }

    // ====== Final measurement of all data qubits ======
    circuit_.safe_append_u("M", all_data, {});

    // Final data detectors: for each Z-stabilizer ancilla in the patches,
    // combine data measurements with the last ancilla measurement
    uint32_t data_round = merge_rounds_ + 2;

    for (const auto& stab : patch_a_stabilizers_) {
        if (stab.is_x_type) continue;
        const auto& aq = qubits_[stab.ancilla];

        std::vector<uint32_t> det_targets;
        for (uint32_t d : stab.data_qubits) {
            for (size_t k = 0; k < all_data.size(); k++) {
                if (all_data[k] == d) {
                    det_targets.push_back(drec(-(int32_t)(all_data.size() - k)));
                    break;
                }
            }
        }
        for (size_t k = 0; k < z_ancillas.size(); k++) {
            if (z_ancillas[k] == stab.ancilla) {
                det_targets.push_back(drec(-(int32_t)(all_data.size() + post_total - k)));
                break;
            }
        }
        circuit_.safe_append_u("DETECTOR", det_targets, {aq.x, aq.y, static_cast<double>(data_round)});
    }

    for (const auto& stab : patch_b_stabilizers_) {
        if (stab.is_x_type) continue;
        const auto& aq = qubits_[stab.ancilla];

        std::vector<uint32_t> det_targets;
        for (uint32_t d : stab.data_qubits) {
            for (size_t k = 0; k < all_data.size(); k++) {
                if (all_data[k] == d) {
                    det_targets.push_back(drec(-(int32_t)(all_data.size() - k)));
                    break;
                }
            }
        }
        for (size_t k = 0; k < z_ancillas.size(); k++) {
            if (z_ancillas[k] == stab.ancilla) {
                det_targets.push_back(drec(-(int32_t)(all_data.size() + post_total - k)));
                break;
            }
        }
        circuit_.safe_append_u("DETECTOR", det_targets, {aq.x, aq.y, static_cast<double>(data_round)});
    }

    // Logical observables
    std::vector<uint32_t> obs_a_targets;
    for (size_t k = 0; k < all_data.size(); k++) {
        const auto& q = qubits_[all_data[k]];
        if (q.patch == DPatchID::PATCH_A && std::abs(q.x - 0.5) < 0.1) {
            obs_a_targets.push_back(drec(-(int32_t)(all_data.size() - k)));
        }
    }
    if (!obs_a_targets.empty()) {
        circuit_.safe_append_u("OBSERVABLE_INCLUDE", obs_a_targets, {0});
    }

    double max_data_x = 0;
    for (const auto& q : qubits_) {
        if (q.type == DQubitType::DATA && q.patch == DPatchID::PATCH_B) {
            max_data_x = std::max(max_data_x, q.x);
        }
    }
    std::vector<uint32_t> obs_b_targets;
    for (size_t k = 0; k < all_data.size(); k++) {
        const auto& q = qubits_[all_data[k]];
        if (q.patch == DPatchID::PATCH_B && std::abs(q.x - max_data_x) < 0.1) {
            obs_b_targets.push_back(drec(-(int32_t)(all_data.size() - k)));
        }
    }
    if (!obs_b_targets.empty()) {
        circuit_.safe_append_u("OBSERVABLE_INCLUDE", obs_b_targets, {1});
    }
}

stim::Circuit DistributedLatticeSurgeryCircuit::generate() {
    generate_circuit();
    return circuit_;
}

size_t DistributedLatticeSurgeryCircuit::num_data_qubits() const {
    size_t count = 0;
    for (const auto& q : qubits_) {
        if (q.type == DQubitType::DATA) count++;
    }
    return count;
}

size_t DistributedLatticeSurgeryCircuit::num_detectors() const {
    return circuit_.count_detectors();
}

std::string DistributedLatticeSurgeryCircuit::annotated_stim_str() const {
    // Post-process circuit_.str() to insert #!pragma POLYGON annotations
    // for stabilizer coloring in the Stim interactive viewer.
    //
    // Colors:
    //   Blue    (0,0,1,0.25)   = Patch X-stabilizer
    //   Red     (1,0,0,0.25)   = Patch Z-stabilizer
    //   Green   (0,0.7,0,0.3)  = Merge X-stabilizer
    //   Orange  (1,0.5,0,0.3)  = Merge Z-stabilizer
    //   Magenta (1,0,1,0.5)    = Remote CNOT edge (one per ancilla↔Patch-B data qubit)

    // Helper: generate a POLYGON pragma for one stabilizer
    auto make_polygon = [this](const DStabilizer& stab, bool is_merge) -> std::string {
        const char* color;
        if (is_merge) {
            color = stab.is_x_type ? "0,0.7,0,0.3" : "1,0.5,0,0.3";
        } else {
            color = stab.is_x_type ? "0,0,1,0.25" : "1,0,0,0.25";
        }

        double ax = qubits_[stab.ancilla].x;
        double ay = qubits_[stab.ancilla].y;

        // Map data qubits by their relative position to the ancilla
        std::map<std::pair<int,int>, uint32_t> pos_map;
        for (uint32_t q : stab.data_qubits) {
            int dx = static_cast<int>(std::round((qubits_[q].x - ax) * 2));
            int dy = static_cast<int>(std::round((qubits_[q].y - ay) * 2));
            pos_map[{dx, dy}] = q;
        }

        // Order clockwise: SW, SE, NE, NW
        std::vector<std::pair<int,int>> cw = {{-1,-1}, {1,-1}, {1,1}, {-1,1}};
        std::ostringstream oss;
        oss << "#!pragma POLYGON(" << color << ")";
        for (auto& [ox, oy] : cw) {
            auto it = pos_map.find({ox, oy});
            if (it != pos_map.end()) {
                oss << " " << it->second;
            }
        }
        return oss.str();
    };

    // Build pragma blocks
    std::string patch_pragmas;
    for (const auto& stab : patch_a_stabilizers_) {
        patch_pragmas += make_polygon(stab, false) + "\n";
    }
    for (const auto& stab : patch_b_stabilizers_) {
        patch_pragmas += make_polygon(stab, false) + "\n";
    }

    // Seam stabilizers use patch color (blue for X)
    std::string seam_pragmas;
    for (const auto& stab : seam_a_stabilizers_) {
        seam_pragmas += make_polygon(stab, false) + "\n";
    }
    for (const auto& stab : seam_b_stabilizers_) {
        seam_pragmas += make_polygon(stab, false) + "\n";
    }

    std::string merge_pragmas;
    for (const auto& stab : merge_stabilizers_) {
        merge_pragmas += make_polygon(stab, true) + "\n";
    }

    // Remote CNOT edges: one 2-qubit polygon per (merge ancilla, Patch B data) pair
    std::string remote_cx_pragmas;
    for (const auto& stab : merge_stabilizers_) {
        for (uint32_t dq : stab.data_qubits) {
            if (qubits_[dq].patch == DPatchID::PATCH_B) {
                std::ostringstream oss;
                oss << "#!pragma POLYGON(1,0,1,0.5) " << stab.ancilla << " " << dq;
                remote_cx_pragmas += oss.str() + "\n";
            }
        }
    }

    // Split circuit text into lines
    std::string base = circuit_.str();
    std::vector<std::string> lines;
    std::istringstream iss(base);
    std::string line;
    while (std::getline(iss, line)) {
        lines.push_back(line);
    }

    // Walk through lines and insert pragmas at round boundaries.
    // Detection: a line starting with "R " or "RX " that follows a "TICK" line
    // marks the start of a new reset group:
    //   Group 0 = data reset (no pragmas)
    //   Group 1 = pre-merge round (patch pragmas already shown in initial block)
    //   Groups 2..1+merge_rounds = merge rounds (patch + merge pragmas)
    //   Group 2+merge_rounds = post-merge round (patch pragmas)
    std::ostringstream out;
    int reset_group = -1;
    bool prev_was_tick = false;
    bool in_qubit_coords = true;

    for (size_t i = 0; i < lines.size(); i++) {
        const auto& l = lines[i];

        // After QUBIT_COORDS block ends, insert initial patch + seam pragmas
        if (in_qubit_coords &&
            (l.size() < 12 || l.substr(0, 12) != "QUBIT_COORDS")) {
            in_qubit_coords = false;
            out << patch_pragmas << seam_pragmas;
        }

        // Detect start of a reset group
        bool is_reset = (l.size() >= 2 && l[0] == 'R' && l[1] == ' ') ||
                         (l.size() >= 3 && l[0] == 'R' && l[1] == 'X' && l[2] == ' ');
        if (is_reset && prev_was_tick) {
            reset_group++;
            if (reset_group >= 2 &&
                reset_group <= 1 + static_cast<int>(merge_rounds_)) {
                // Before a merge round: show patch + merge stabilizers + remote CX edges
                out << patch_pragmas << merge_pragmas << remote_cx_pragmas;
            } else if (reset_group == 2 + static_cast<int>(merge_rounds_)) {
                // Before post-merge round: show patch + seam stabilizers (no merge)
                out << patch_pragmas << seam_pragmas;
            }
        }

        prev_was_tick = (l == "TICK");
        out << l << "\n";
    }

    return out.str();
}

} // namespace bucket_sim
