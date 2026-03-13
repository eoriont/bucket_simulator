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
    if (!config_.superstabilizers.empty() && merge_type_ != MergeType::XX_MERGE_DISTRIBUTED) {
        throw std::invalid_argument("superstabilizers is only supported for XX_MERGE_DISTRIBUTED");
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

    // ========== DATA QUBIT REMOVAL for superstabilizer positions ==========
    // Disable data qubits at the exact (x, y) positions listed in superstabilizers config.
    for (const auto& [tx, ty] : config_.superstabilizers) {
        bool found = false;
        for (const auto& q : qubits_) {
            if (q.type == DQubitType::DATA &&
                std::abs(q.x - tx) < 0.1 && std::abs(q.y - ty) < 0.1) {
                removed_data_qubits_.insert(q.index);
                std::cerr << "  Disabling data qubit " << q.index
                          << " at (" << q.x << ", " << q.y << ") for superstabilizer" << std::endl;
                found = true;
                break;
            }
        }
        if (!found) {
            std::cerr << "  Warning: no data qubit found at superstabilizer position ("
                      << tx << ", " << ty << ")" << std::endl;
        }
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
            if (removed_data_qubits_.count(q.index)) continue;
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

        // hi
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

    // Detect weight-1 merge-Z ancillas (become weight-1 after data qubit removal → suppress)
    for (const auto& stab : merge_stabilizers_) {
        if (!stab.is_x_type && stab.data_qubits.size() <= 1) {
            merge_suppressed_z_ancillas_.insert(stab.ancilla);
            std::cerr << "  Suppressing weight-1 merge-Z ancilla " << stab.ancilla
                      << " at (" << qubits_[stab.ancilla].x << ", " << qubits_[stab.ancilla].y << ")" << std::endl;
        }
    }

    // Collect partner data qubits of suppressed merge-Z ancillas.
    // When a merge-Z becomes weight-1, its surviving data qubit (the Patch A side)
    // loses all Z-type syndrome coverage during merge rounds. Exclude it from merge
    // CX layers so no stabilizer incorrectly couples to it during merge.
    for (const auto& stab : merge_stabilizers_) {
        if (!stab.is_x_type && merge_suppressed_z_ancillas_.count(stab.ancilla)) {
            for (uint32_t dq : stab.data_qubits) {
                merge_excluded_data_qubits_.insert(dq);
                std::cerr << "  Excluding partner data qubit " << dq
                          << " at (" << qubits_[dq].x << ", " << qubits_[dq].y
                          << ") from merge rounds (lost Z coverage)" << std::endl;
            }
        }
    }

    // Populate the general suppressed-ancilla set.
    // Seed with all suppressed Z ancillas, then add any merge ancilla (X or Z)
    // that touches a merge_excluded data qubit — these are the weight-4 ancillas
    // that touched both deleted qubits in the boundary cascade (step 3 of
    // surface_general_defect's handle_boundary_data logic).
    merge_suppressed_ancillas_ = merge_suppressed_z_ancillas_;
    for (const auto& stab : merge_stabilizers_) {
        if (merge_suppressed_ancillas_.count(stab.ancilla)) continue;
        for (uint32_t dq : stab.data_qubits) {
            if (merge_excluded_data_qubits_.count(dq)) {
                merge_suppressed_ancillas_.insert(stab.ancilla);
                std::cerr << "  Suppressing merge " << (stab.is_x_type ? "X" : "Z")
                          << " ancilla " << stab.ancilla
                          << " at (" << qubits_[stab.ancilla].x << ", " << qubits_[stab.ancilla].y
                          << ") (touched excluded data qubit " << dq << ")" << std::endl;
                break;
            }
        }
    }

    // Extended cascade for patch ancillas.
    // A patch ancilla adjacent to a removed boundary data qubit can become weight-1 in
    // merge rounds (build_stabilizers uses find_data_qubit which already excludes
    // removed_data_qubits_). Suppress it in merge rounds and cascade its surviving
    // partner data qubit into merge_excluded_data_qubits_, then suppress any patch
    // ancilla touching that excluded qubit (same 3-step logic as for merge ancillas).
    std::vector<std::vector<DStabilizer>*> patch_lists = {&patch_a_stabilizers_, &patch_b_stabilizers_};

    // Step 1 extended: weight-1 patch ancillas in merge rounds
    for (auto* lst : patch_lists) {
        for (const auto& stab : *lst) {
            if (stab.data_qubits.size() != 1) continue;
            if (merge_suppressed_ancillas_.count(stab.ancilla)) continue;
            merge_suppressed_ancillas_.insert(stab.ancilla);
            merge_excluded_data_qubits_.insert(stab.data_qubits[0]);
            std::cerr << "  Suppressing weight-1 patch " << (stab.is_x_type ? "X" : "Z")
                      << " ancilla " << stab.ancilla
                      << " at (" << qubits_[stab.ancilla].x << ", " << qubits_[stab.ancilla].y
                      << "), excluding partner data qubit " << stab.data_qubits[0]
                      << " at (" << qubits_[stab.data_qubits[0]].x << ", " << qubits_[stab.data_qubits[0]].y << ")" << std::endl;
        }
    }
    // Step 3 extended: suppress patch ancillas touching BOTH a removed qubit AND a
    // merge-excluded qubit. This mirrors the reference cascade (handle_boundary_data
    // step 3): only the weight-4 syndrome that touches both deleted data qubits is
    // deleted. An ancilla touching only the excluded partner (but not the originally
    // removed qubit) is left active.
    // Because find_data_qubit already strips removed qubits from stab.data_qubits,
    // we check the 4 diagonal neighbors directly for removed qubit membership.
    const std::pair<double,double> diag_dirs[4] = {{0.5,-0.5},{0.5,0.5},{-0.5,-0.5},{-0.5,0.5}};
    auto touches_removed = [&](const DStabilizer& stab) -> bool {
        const auto& aq = qubits_[stab.ancilla];
        for (auto [dx, dy] : diag_dirs) {
            double tx = aq.x + dx, ty = aq.y + dy;
            for (uint32_t ridx : removed_data_qubits_) {
                const auto& rq = qubits_[ridx];
                if (std::abs(rq.x - tx) < 0.1 && std::abs(rq.y - ty) < 0.1) return true;
            }
        }
        return false;
    };
    for (auto* lst : patch_lists) {
        for (const auto& stab : *lst) {
            if (merge_suppressed_ancillas_.count(stab.ancilla)) continue;
            bool hits_excluded = false;
            for (uint32_t dq : stab.data_qubits)
                if (merge_excluded_data_qubits_.count(dq)) { hits_excluded = true; break; }
            if (!hits_excluded) continue;
            if (!touches_removed(stab)) continue;  // must also touch the originally removed qubit
            merge_suppressed_ancillas_.insert(stab.ancilla);
            std::cerr << "  Suppressing patch " << (stab.is_x_type ? "X" : "Z")
                      << " ancilla " << stab.ancilla
                      << " at (" << qubits_[stab.ancilla].x << ", " << qubits_[stab.ancilla].y
                      << ") (touches both removed and excluded data qubits)" << std::endl;
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
                  << (s.crosses_boundary ? " [CROSS-BOUNDARY]" : "")
                  << (merge_suppressed_ancillas_.count(s.ancilla) ? " [SUPPRESSED]" : "") << std::endl;
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

    std::vector<uint32_t> ss_data; // Superstabilizer data qubits: active in pre/post, idle in merge
    for (const auto& q : qubits_) {
        if (q.type == DQubitType::DATA) {
            if (removed_data_qubits_.count(q.index)) {
                ss_data.push_back(q.index);
                all_data.push_back(q.index);
                continue;
            }
            all_data.push_back(q.index);
            if (q.patch == DPatchID::PATCH_A) patch_a_data.push_back(q.index);
            else patch_b_data.push_back(q.index);
        } else if (q.type == DQubitType::Z_ANCILLA) {
            if (q.patch == DPatchID::MERGE) {
                if (!merge_suppressed_ancillas_.count(q.index))
                    merge_z_ancillas.push_back(q.index);
            } else {
                z_ancillas.push_back(q.index);
            }
        } else if (q.type == DQubitType::X_ANCILLA) {
            if (q.patch == DPatchID::MERGE) {
                if (!merge_suppressed_ancillas_.count(q.index)) merge_x_ancillas.push_back(q.index);
                // suppressed merge X ancillas are dropped entirely
            } else if (q.patch == DPatchID::SEAM_A) seam_a_x_ancillas.push_back(q.index);
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

    // === Superstabilizer: data qubit removal approach ===
    // Patch B data qubits at SS positions are removed in initialize_layout().
    // The merge-X ancilla at each SS y still fires normally (now weight-3 with 1 remote CNOT).
    // No gauge splitting — all merge-X ancillas and all seam positions are treated uniformly.
    // (Gauge-split approach commented out; see git history if needed.)
    const size_t ss_count = 0;
    std::vector<uint32_t> normal_merge_x = merge_x_ancillas;

    // Match each surviving merge-X ancilla to its seam counterpart by y-coordinate.
    // When a merge-X is suppressed its y position is absent from normal_merge_x,
    // so seam positions at that y become "orphans" with no merge-X partner.
    std::vector<size_t> normal_seam_indices(normal_merge_x.size(), SIZE_MAX);
    std::vector<size_t> orphan_seam_indices;
    for (size_t j = 0; j < seam_a_x_ancillas.size(); j++) {
        double seam_y = qubits_[seam_a_x_ancillas[j]].y;
        bool matched = false;
        for (size_t m = 0; m < normal_merge_x.size(); m++) {
            if (std::abs(qubits_[normal_merge_x[m]].y - seam_y) < 0.1) {
                normal_seam_indices[m] = j;
                matched = true;
                break;
            }
        }
        if (!matched) orphan_seam_indices.push_back(j);
    }

    std::vector<uint32_t> cx_layer1, cx_layer2, cx_layer3, cx_layer4;           // merge rounds (excludes ss_data)
    std::vector<uint32_t> pre_cx_layer1, pre_cx_layer2, pre_cx_layer3, pre_cx_layer4; // pre/post rounds (includes ss_data)
    std::vector<uint32_t> seam_cx1, seam_cx2, seam_cx3, seam_cx4;
    std::vector<uint32_t> pre_seam_cx1, pre_seam_cx2, pre_seam_cx3, pre_seam_cx4;   // seam CXs including ss_data neighbors
    std::vector<uint32_t> merge_cx1, merge_cx2, merge_cx3, merge_cx4;
    // Gauge CX layers: for interior superstabilizer ancillas (weight-3 due to a
    // disabled-but-not-removed data qubit that is NOT on the code boundary).
    // X-gauges fire in the same half-round as normal X ancillas.
    // Z-gauges fire in the same half-round as normal Z ancillas.
    std::vector<uint32_t> x_gauge_cx1, x_gauge_cx2, x_gauge_cx3, x_gauge_cx4;
    std::vector<uint32_t> z_gauge_cx1, z_gauge_cx2, z_gauge_cx3, z_gauge_cx4;

    auto find_data_at = [this](double ax, double ay, double dx, double dy, DPatchID owner) -> int32_t {
        double tx = ax + dx, ty = ay + dy;
        for (const auto& q : qubits_) {
            if (q.type != DQubitType::DATA) continue;
            if (removed_data_qubits_.count(q.index)) continue;
            if (merge_excluded_data_qubits_.count(q.index)) continue;
            if (std::abs(q.x - tx) > 0.1 || std::abs(q.y - ty) > 0.1) continue;
            if (owner == DPatchID::SEAM_A && q.patch != DPatchID::PATCH_A) continue;
            if (owner == DPatchID::SEAM_B && q.patch != DPatchID::PATCH_B) continue;
            return static_cast<int32_t>(q.index);
        }
        return -1;
    };

    // Same as find_data_at but includes removed (superstabilizer) data qubits — used for pre/post rounds.
    auto find_data_at_pre = [this](double ax, double ay, double dx, double dy, DPatchID owner) -> int32_t {
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

    // Returns true if any diagonal neighbor of ancilla (ax, ay) is a removed or
    // merge-excluded data qubit that lies on the code boundary. These ancillas change
    // weight between pre/post rounds (full weight) and merge rounds (reduced weight),
    // so their transition detectors are non-deterministic and must be skipped.
    auto has_boundary_ss_neighbor = [this](double ax, double ay) -> bool {
        const std::pair<double,double> dirs[4] = {{0.5,-0.5},{0.5,0.5},{-0.5,-0.5},{-0.5,0.5}};
        for (auto [dx, dy] : dirs) {
            double tx = ax + dx, ty = ay + dy;
            for (const auto& q : qubits_) {
                if (q.type != DQubitType::DATA) continue;
                if (std::abs(q.x - tx) > 0.1 || std::abs(q.y - ty) > 0.1) continue;
                if (removed_data_qubits_.count(q.index) || merge_excluded_data_qubits_.count(q.index)) {
                    bool on_boundary = (std::abs(q.y - 0.5) < 0.1 ||
                                        std::abs(q.y - (distance_ - 0.5)) < 0.1);
                    if (on_boundary) return true;
                }
            }
        }
        return false;
    };

    // Returns true if any diagonal neighbor of ancilla (ax, ay) is an interior
    // disabled data qubit — i.e., a superstabilizer position not on the code boundary.
    // "Interior" means the qubit's y coordinate is strictly between the top and bottom
    // data rows of its patch (not 0.5 and not d-0.5).
    auto has_interior_ss_neighbor = [this](double ax, double ay) -> bool {
        const std::pair<double,double> dirs[4] = {{0.5,-0.5},{0.5,0.5},{-0.5,-0.5},{-0.5,0.5}};
        for (auto [dx, dy] : dirs) {
            double tx = ax + dx, ty = ay + dy;
            for (uint32_t idx : removed_data_qubits_) {
                const auto& dq = qubits_[idx];
                if (std::abs(dq.x - tx) > 0.1 || std::abs(dq.y - ty) > 0.1) continue;
                bool on_boundary = (std::abs(dq.y - 0.5) < 0.1 ||
                                    std::abs(dq.y - (distance_ - 0.5)) < 0.1);
                if (!on_boundary) return true;
            }
        }
        return false;
    };

    for (const auto& q : qubits_) {
        if (q.type == DQubitType::DATA) continue;

        bool is_merge = (q.patch == DPatchID::MERGE);
        bool is_seam = (q.patch == DPatchID::SEAM_A || q.patch == DPatchID::SEAM_B);
        bool is_x = (q.type == DQubitType::X_ANCILLA);

        // Skip suppressed ancillas in merge rounds: applies to merge ancillas and any
        // patch ancilla suppressed by the extended boundary cascade.
        if (merge_suppressed_ancillas_.count(q.index)) continue;

        // Ancillas (patch or merge) adjacent to a disabled (non-boundary) SS data qubit
        // are gauges: their CX pairs go into separate gauge layers.
        bool is_gauge = !is_seam && has_interior_ss_neighbor(q.x, q.y);

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

            auto& target_layer = is_seam ?
                (layer == 0 ? seam_cx1 : layer == 1 ? seam_cx2 : layer == 2 ? seam_cx3 : seam_cx4) :
                is_gauge && is_x ?
                (layer == 0 ? x_gauge_cx1 : layer == 1 ? x_gauge_cx2 : layer == 2 ? x_gauge_cx3 : x_gauge_cx4) :
                is_gauge ?
                (layer == 0 ? z_gauge_cx1 : layer == 1 ? z_gauge_cx2 : layer == 2 ? z_gauge_cx3 : z_gauge_cx4) :
                is_merge ?
                (layer == 0 ? merge_cx1 : layer == 1 ? merge_cx2 : layer == 2 ? merge_cx3 : merge_cx4) :
                (layer == 0 ? cx_layer1 : layer == 1 ? cx_layer2 : layer == 2 ? cx_layer3 : cx_layer4);

            target_layer.push_back(ctrl);
            target_layer.push_back(tgt);
        }
    }

    // Build pre/post CX layers: same as above but using find_data_at_pre (includes ss_data qubits).
    // Only patch and seam ancillas participate in pre/post rounds; merge ancillas are excluded.
    for (const auto& q : qubits_) {
        if (q.type == DQubitType::DATA) continue;
        bool is_merge = (q.patch == DPatchID::MERGE);
        bool is_seam  = (q.patch == DPatchID::SEAM_A || q.patch == DPatchID::SEAM_B);
        bool is_x     = (q.type == DQubitType::X_ANCILLA);
        if (is_merge) continue; // merge ancillas not active in pre/post rounds

        std::vector<std::pair<double, double>> dirs = {
            {0.5, -0.5}, {0.5, 0.5}, {-0.5, -0.5}, {-0.5, 0.5}
        };
        for (int layer = 0; layer < 4; layer++) {
            int32_t data_idx = find_data_at_pre(q.x, q.y, dirs[layer].first, dirs[layer].second, q.patch);
            if (data_idx < 0) continue;
            uint32_t ctrl, tgt;
            if (is_x) { ctrl = q.index; tgt = static_cast<uint32_t>(data_idx); }
            else       { ctrl = static_cast<uint32_t>(data_idx); tgt = q.index; }
            auto& tl = is_seam ?
                (layer == 0 ? pre_seam_cx1 : layer == 1 ? pre_seam_cx2 : layer == 2 ? pre_seam_cx3 : pre_seam_cx4) :
                (layer == 0 ? pre_cx_layer1 : layer == 1 ? pre_cx_layer2 : layer == 2 ? pre_cx_layer3 : pre_cx_layer4);
            tl.push_back(ctrl);
            tl.push_back(tgt);
        }
    }

    // No CX filtering needed: all merge-X ancillas are active (ss_count=0).
    auto norm_mcx1 = merge_cx1;
    auto norm_mcx2 = merge_cx2;
    auto norm_mcx3 = merge_cx3;
    auto norm_mcx4 = merge_cx4;

    // Extract gauge ancillas from patch and merge ancilla lists.
    // A gauge ancilla (patch or merge) is adjacent to a disabled (non-boundary) SS qubit.
    // It is removed from the normal lists and placed in x_gauge_ancillas / z_gauge_ancillas.
    std::vector<uint32_t> x_gauge_ancillas, z_gauge_ancillas;
    {
        std::vector<uint32_t> new_x, new_z;
        for (uint32_t idx : x_ancillas) {
            if (has_interior_ss_neighbor(qubits_[idx].x, qubits_[idx].y))
                x_gauge_ancillas.push_back(idx);
            else
                new_x.push_back(idx);
        }
        for (uint32_t idx : z_ancillas) {
            if (has_interior_ss_neighbor(qubits_[idx].x, qubits_[idx].y))
                z_gauge_ancillas.push_back(idx);
            else
                new_z.push_back(idx);
        }
        x_ancillas = std::move(new_x);
        z_ancillas = std::move(new_z);

        // Also extract merge ancillas that neighbor an SS position.
        std::vector<uint32_t> new_mz, new_mx;
        for (uint32_t idx : merge_z_ancillas) {
            if (has_interior_ss_neighbor(qubits_[idx].x, qubits_[idx].y))
                z_gauge_ancillas.push_back(idx);
            else
                new_mz.push_back(idx);
        }
        for (uint32_t idx : normal_merge_x) {
            if (has_interior_ss_neighbor(qubits_[idx].x, qubits_[idx].y))
                x_gauge_ancillas.push_back(idx);
            else
                new_mx.push_back(idx);
        }
        merge_z_ancillas = std::move(new_mz);
        normal_merge_x   = std::move(new_mx);
    }
    bool has_gauges = !x_gauge_ancillas.empty() || !z_gauge_ancillas.empty();
    has_gauges_ = has_gauges;

    // Boundary-weight-changing ancillas: patch ancillas (not gauge, not seam, not merge)
    // that neighbor a boundary SS qubit. Their stabilizer weight differs between pre/post
    // rounds and merge rounds, so the transition detectors (r==0 entry, post-merge exit)
    // are non-deterministic and must be skipped.
    std::unordered_set<uint32_t> boundary_changing_x_set, boundary_changing_z_set;
    for (uint32_t idx : x_ancillas) {
        if (has_boundary_ss_neighbor(qubits_[idx].x, qubits_[idx].y))
            boundary_changing_x_set.insert(idx);
    }
    for (uint32_t idx : z_ancillas) {
        if (has_boundary_ss_neighbor(qubits_[idx].x, qubits_[idx].y))
            boundary_changing_z_set.insert(idx);
    }
    if (!boundary_changing_x_set.empty() || !boundary_changing_z_set.empty()) {
        std::cerr << "  Boundary-weight-changing X ancillas (" << boundary_changing_x_set.size() << "):";
        for (uint32_t idx : boundary_changing_x_set)
            std::cerr << " " << idx << "(" << qubits_[idx].x << "," << qubits_[idx].y << ")";
        std::cerr << "\n";
        std::cerr << "  Boundary-weight-changing Z ancillas (" << boundary_changing_z_set.size() << "):";
        for (uint32_t idx : boundary_changing_z_set)
            std::cerr << " " << idx << "(" << qubits_[idx].x << "," << qubits_[idx].y << ")";
        std::cerr << "\n";
    }
    std::cerr << "  X gauge ancillas (" << x_gauge_ancillas.size() << "):" << std::endl;
    for (uint32_t idx : x_gauge_ancillas)
        std::cerr << "    idx=" << idx << " coords=(" << qubits_[idx].x << "," << qubits_[idx].y << ")" << std::endl;
    std::cerr << "  Z gauge ancillas (" << z_gauge_ancillas.size() << "):" << std::endl;
    for (uint32_t idx : z_gauge_ancillas)
        std::cerr << "    idx=" << idx << " coords=(" << qubits_[idx].x << "," << qubits_[idx].y << ")" << std::endl;

    // Separate patch gauge ancillas (used in pre/post rounds as normal ancillas)
    // from merge gauge ancillas (only active during merge rounds).
    std::vector<uint32_t> patch_x_gauge, patch_z_gauge;
    for (uint32_t idx : x_gauge_ancillas)
        if (qubits_[idx].patch != DPatchID::MERGE) patch_x_gauge.push_back(idx);
    for (uint32_t idx : z_gauge_ancillas)
        if (qubits_[idx].patch != DPatchID::MERGE) patch_z_gauge.push_back(idx);

    // Counts for offset arithmetic
    size_t patch_z_count = z_ancillas.size();
    size_t patch_x_count = x_ancillas.size();
    size_t seam_a_count = seam_a_x_ancillas.size();
    size_t seam_b_count = seam_b_x_ancillas.size();
    size_t seam_count = seam_a_count + seam_b_count;
    size_t merge_z_count = merge_z_ancillas.size();
    size_t merge_x_count = merge_x_ancillas.size();
    size_t normal_merge_x_count = normal_merge_x.size(); // = merge_x_count - ss_count
    size_t xg_count = x_gauge_ancillas.size();
    size_t zg_count = z_gauge_ancillas.size();
    size_t patch_xg_count = patch_x_gauge.size();
    size_t patch_zg_count = patch_z_gauge.size();

    // Pre-merge round measurement order:
    //   M(z_ancillas + patch_z_gauges), MX(x_ancillas + patch_x_gauges + seam)
    // Patch gauge ancillas participate as full-weight stabilizers (SS data qubits are
    // active in pre/post rounds) — no half-round splitting needed.
    size_t pre_total = patch_z_count + patch_zg_count + patch_x_count + patch_xg_count + seam_count;
    // When gauges are present each merge round is split into two half-rounds:
    //   half1: M(z + merge_z), MX(x + merge_x + x_gauges)  — x-gauges fire and are measured
    //   half2: M(z + merge_z + z_gauges), MX(x + merge_x)  — z-gauges fire and are measured
    size_t pre_h1  = patch_z_count + patch_x_count + seam_count + xg_count;
    size_t pre_h2  = patch_z_count + zg_count       + patch_x_count + seam_count;
    size_t merge_h1 = patch_z_count + merge_z_count + patch_x_count + normal_merge_x_count + xg_count;
    size_t merge_h2 = patch_z_count + merge_z_count + zg_count + patch_x_count + normal_merge_x_count;
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
        std::vector<uint32_t> all_seam_x;
        all_seam_x.insert(all_seam_x.end(), seam_a_x_ancillas.begin(), seam_a_x_ancillas.end());
        all_seam_x.insert(all_seam_x.end(), seam_b_x_ancillas.begin(), seam_b_x_ancillas.end());

        auto emit_pre_layer = [&](const std::vector<uint32_t>& int_cx,
                                   const std::vector<uint32_t>& seam_cx_layer,
                                   const std::vector<uint32_t>& gauge_cx) {
            std::vector<uint32_t> combined = int_cx;
            combined.insert(combined.end(), seam_cx_layer.begin(), seam_cx_layer.end());
            combined.insert(combined.end(), gauge_cx.begin(), gauge_cx.end());
            circuit_.safe_append_u("TICK", {}, {});
            if (!combined.empty()) circuit_.safe_append_u("CX", combined, {});
        };

        {
            // In pre/post rounds, SS data qubits are active, so patch gauge ancillas
            // behave as normal full-weight stabilizers — include them without splitting.
            std::vector<uint32_t> all_z_pre = z_ancillas;
            all_z_pre.insert(all_z_pre.end(), patch_z_gauge.begin(), patch_z_gauge.end());

            std::vector<uint32_t> all_x_pre = x_ancillas;
            all_x_pre.insert(all_x_pre.end(), patch_x_gauge.begin(), patch_x_gauge.end());
            all_x_pre.insert(all_x_pre.end(), all_seam_x.begin(), all_seam_x.end());

            circuit_.safe_append_u("TICK", {}, {});
            circuit_.safe_append_u("R", all_z_pre, {});
            circuit_.safe_append_u("RX", all_x_pre, {});

            emit_pre_layer(pre_cx_layer1, pre_seam_cx1, {});
            emit_pre_layer(pre_cx_layer2, pre_seam_cx2, {});
            emit_pre_layer(pre_cx_layer3, pre_seam_cx3, {});
            emit_pre_layer(pre_cx_layer4, pre_seam_cx4, {});

            circuit_.safe_append_u("TICK", {}, {});
            circuit_.safe_append_u("M", all_z_pre, {});
            circuit_.safe_append_u("MX", all_x_pre, {});

            for (size_t i = 0; i < patch_z_count; i++) {
                const auto& q = qubits_[z_ancillas[i]];
                circuit_.safe_append_u("DETECTOR", {drec(-(int32_t)(pre_total - i))}, {q.x, q.y, 0});
            }
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

        auto emit_merge_layer = [&](const std::vector<uint32_t>& int_cx,
                                    const std::vector<uint32_t>& nm_cx,
                                    const std::vector<uint32_t>& gauge_cx) {
            std::vector<uint32_t> combined = int_cx;
            combined.insert(combined.end(), nm_cx.begin(), nm_cx.end());
            combined.insert(combined.end(), gauge_cx.begin(), gauge_cx.end());
            if (!combined.empty()) circuit_.safe_append_u("CX", combined, {});
            circuit_.safe_append_u("TICK", {}, {});
        };

        if (!has_gauges) {
            // ---- No gauges: original single-half structure ----
            std::vector<uint32_t> all_z_reset = z_ancillas;
            all_z_reset.insert(all_z_reset.end(), merge_z_ancillas.begin(), merge_z_ancillas.end());
            circuit_.safe_append_u("R", all_z_reset, {});

            std::vector<uint32_t> all_x_reset = x_ancillas;
            all_x_reset.insert(all_x_reset.end(), normal_merge_x.begin(), normal_merge_x.end());
            circuit_.safe_append_u("RX", all_x_reset, {});
            circuit_.safe_append_u("TICK", {}, {});

            emit_merge_layer(cx_layer1, norm_mcx1, {});
            emit_merge_layer(cx_layer2, norm_mcx2, {});
            emit_merge_layer(cx_layer3, norm_mcx3, {});
            emit_merge_layer(cx_layer4, norm_mcx4, {});

            std::vector<uint32_t> all_z_meas = z_ancillas;
            all_z_meas.insert(all_z_meas.end(), merge_z_ancillas.begin(), merge_z_ancillas.end());
            circuit_.safe_append_u("M", all_z_meas, {});

            std::vector<uint32_t> all_x_meas = x_ancillas;
            all_x_meas.insert(all_x_meas.end(), normal_merge_x.begin(), normal_merge_x.end());
            circuit_.safe_append_u("MX", all_x_meas, {});
        } else {
            // ---- Gauges present: two half-rounds ----
            // Half 1 (X-gauge half):
            {
                std::vector<uint32_t> z_reset = z_ancillas;
                z_reset.insert(z_reset.end(), merge_z_ancillas.begin(), merge_z_ancillas.end());
                z_reset.insert(z_reset.end(), z_gauge_ancillas.begin(), z_gauge_ancillas.end());

                std::vector<uint32_t> x_reset = x_ancillas;
                x_reset.insert(x_reset.end(), normal_merge_x.begin(), normal_merge_x.end());
                // On the very first merge round, x_gauge ancillas have no prior half 2
                // to carry state from, so reset them here too.
                if (r == 0)
                    x_reset.insert(x_reset.end(), x_gauge_ancillas.begin(), x_gauge_ancillas.end());

                circuit_.safe_append_u("R", z_reset, {});
                circuit_.safe_append_u("RX", x_reset, {});
                circuit_.safe_append_u("TICK", {}, {});

                emit_merge_layer(cx_layer1, norm_mcx1, x_gauge_cx1);
                emit_merge_layer(cx_layer2, norm_mcx2, x_gauge_cx2);
                emit_merge_layer(cx_layer3, norm_mcx3, x_gauge_cx3);
                emit_merge_layer(cx_layer4, norm_mcx4, x_gauge_cx4);

                // Half 1: measure z_anc + merge_z (NOT z_gauges — they carry state into half 2)
                std::vector<uint32_t> z_meas = z_ancillas;
                z_meas.insert(z_meas.end(), merge_z_ancillas.begin(), merge_z_ancillas.end());
                circuit_.safe_append_u("M", z_meas, {});

                std::vector<uint32_t> x_meas = x_ancillas;
                x_meas.insert(x_meas.end(), normal_merge_x.begin(), normal_merge_x.end());
                x_meas.insert(x_meas.end(), x_gauge_ancillas.begin(), x_gauge_ancillas.end());
                circuit_.safe_append_u("MX", x_meas, {});
            }

            // Half 2 (Z-gauge half):
            {
                // Half 2: reset z_anc + merge_z (NOT z_gauges — they already carry state from half 1 reset)
                std::vector<uint32_t> z_reset = z_ancillas;
                z_reset.insert(z_reset.end(), merge_z_ancillas.begin(), merge_z_ancillas.end());

                std::vector<uint32_t> x_reset = x_ancillas;
                x_reset.insert(x_reset.end(), normal_merge_x.begin(), normal_merge_x.end());
                x_reset.insert(x_reset.end(), x_gauge_ancillas.begin(), x_gauge_ancillas.end());

                circuit_.safe_append_u("TICK", {}, {});
                circuit_.safe_append_u("R", z_reset, {});
                circuit_.safe_append_u("RX", x_reset, {});
                circuit_.safe_append_u("TICK", {}, {});

                emit_merge_layer(cx_layer1, norm_mcx1, z_gauge_cx1);
                emit_merge_layer(cx_layer2, norm_mcx2, z_gauge_cx2);
                emit_merge_layer(cx_layer3, norm_mcx3, z_gauge_cx3);
                emit_merge_layer(cx_layer4, norm_mcx4, z_gauge_cx4);

                std::vector<uint32_t> z_meas = z_ancillas;
                z_meas.insert(z_meas.end(), merge_z_ancillas.begin(), merge_z_ancillas.end());
                z_meas.insert(z_meas.end(), z_gauge_ancillas.begin(), z_gauge_ancillas.end());
                circuit_.safe_append_u("M", z_meas, {});

                std::vector<uint32_t> x_meas = x_ancillas;
                x_meas.insert(x_meas.end(), normal_merge_x.begin(), normal_merge_x.end());
                circuit_.safe_append_u("MX", x_meas, {});

            }
        }

        // ---- Detectors ----
        if (!has_gauges) {
            const size_t pz = patch_z_count, mz = merge_z_count, px = patch_x_count;
            const size_t nm = normal_merge_x_count, ss = ss_count;

            // Patch Z: compare with previous round
            // Skip r==0 for boundary-weight-changing ancillas (non-deterministic transition)
            for (size_t i = 0; i < pz; i++) {
                const auto& q = qubits_[z_ancillas[i]];
                if (r == 0 && boundary_changing_z_set.count(z_ancillas[i])) continue;
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
            // Skip r==0 for boundary-weight-changing ancillas (weight-4 pre vs weight-3 merge → non-det)
            for (size_t i = 0; i < px; i++) {
                const auto& q = qubits_[x_ancillas[i]];
                if (r == 0 && boundary_changing_x_set.count(x_ancillas[i])) continue;
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
                size_t j = normal_seam_indices[m];
                if (r == 0) {
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
        } else {
            // Gauge detectors: placed after half2 (end of the full double-half round).
            // Normal Z ancillas (patch + merge): compare current half2 with previous half2.
            // The half2 measurement order is: M(z, merge_z, z_gauges), MX(x, merge_x).
            // merge_full = merge_h1 + merge_h2; pre_full = pre_h1 + pre_h2.
            const size_t merge_full = merge_h1 + merge_h2;
            const size_t pre_full   = pre_h1  + pre_h2;

            // Patch Z: at offset -(merge_h2 - i) from end of half2
            for (size_t i = 0; i < patch_z_count; i++) {
                const auto& q = qubits_[z_ancillas[i]];
                int32_t curr = -(int32_t)(merge_h2 - i);
                int32_t prev = (r == 0)
                    ? -(int32_t)(merge_h2 + merge_h1 + pre_total - i)
                    : -(int32_t)(merge_h2 + merge_h1 + merge_h2 - i);
                circuit_.safe_append_u("DETECTOR", {drec(curr), drec(prev)},
                    {q.x, q.y, static_cast<double>(r + 1)});
            }

            // Merge Z (no detector on r==0, same anti-commutation reason as non-gauge path)
            if (r >= 1) {
                for (size_t i = 0; i < merge_z_count; i++) {
                    const auto& q = qubits_[merge_z_ancillas[i]];
                    int32_t curr = -(int32_t)(merge_h2 - patch_z_count - i);
                    int32_t prev = -(int32_t)(merge_full + merge_h2 - patch_z_count - i);
                    circuit_.safe_append_u("DETECTOR", {drec(curr), drec(prev)},
                        {q.x, q.y, static_cast<double>(r + 1)});
                }
            }

            // Patch X: at offset -(merge_h2 - pz - mz - zg - i) from end of half2
            for (size_t i = 0; i < patch_x_count; i++) {
                const auto& q = qubits_[x_ancillas[i]];
                int32_t curr = -(int32_t)(patch_x_count + normal_merge_x_count - i);
                int32_t prev = (r == 0)
                    ? -(int32_t)(merge_h2 + merge_h1 + pre_total - patch_z_count - patch_zg_count - i)
                    : -(int32_t)(merge_full + patch_x_count + normal_merge_x_count - i);
                circuit_.safe_append_u("DETECTOR", {drec(curr), drec(prev)},
                    {q.x, q.y, static_cast<double>(r + 1)});
            }

            // X-gauge superstabilizer: XOR of all x_gauge measurements in current half1
            // with those in the previous half1.
            // x_gauge[i] in current half1 is at -(merge_h2 + xg_count - i) from end of half2.
            // x_gauge[i] in previous half1: go back another merge_h2 + merge_h1 (prev round).
            // r==0 gauge detectors omitted: first merge half compares weight-3 against
            // the weight-4 pre-round measurement — non-deterministic across the boundary.
            if (!x_gauge_ancillas.empty() && r > 0) {
                std::vector<uint32_t> tgts;
                for (size_t i = 0; i < xg_count; i++) {
                    int32_t curr = -(int32_t)(merge_h2 + xg_count - i);
                    int32_t prev = -(int32_t)(merge_h2 + xg_count - i + merge_full);
                    tgts.push_back(drec(curr));
                    tgts.push_back(drec(prev));
                }
                const auto& q = qubits_[x_gauge_ancillas[0]];
                circuit_.safe_append_u("DETECTOR", tgts, {q.x, q.y, static_cast<double>(r + 1)});
            }

            if (!z_gauge_ancillas.empty() && r > 0) {
                std::vector<uint32_t> tgts;
                for (size_t i = 0; i < zg_count; i++) {
                    int32_t curr = -(int32_t)(zg_count - i);
                    int32_t prev = -(int32_t)(zg_count - i + merge_full);
                    tgts.push_back(drec(curr));
                    tgts.push_back(drec(prev));
                }
                const auto& q = qubits_[z_gauge_ancillas[0]];
                circuit_.safe_append_u("DETECTOR", tgts, {q.x, q.y, static_cast<double>(r + 1)});
            }
        }

        circuit_.safe_append_u("TICK", {}, {});
    }

    // ====== Final round: Post-merge (interior patches + seam, no merge) ======
    {
        std::vector<uint32_t> all_seam_x;
        all_seam_x.insert(all_seam_x.end(), seam_a_x_ancillas.begin(), seam_a_x_ancillas.end());
        all_seam_x.insert(all_seam_x.end(), seam_b_x_ancillas.begin(), seam_b_x_ancillas.end());

        auto emit_post_layer = [&](const std::vector<uint32_t>& int_cx,
                                    const std::vector<uint32_t>& seam_cx_layer,
                                    const std::vector<uint32_t>& gauge_cx) {
            std::vector<uint32_t> combined = int_cx;
            combined.insert(combined.end(), seam_cx_layer.begin(), seam_cx_layer.end());
            combined.insert(combined.end(), gauge_cx.begin(), gauge_cx.end());
            circuit_.safe_append_u("TICK", {}, {});
            if (!combined.empty()) circuit_.safe_append_u("CX", combined, {});
        };

        uint32_t final_round = merge_rounds_ + 1;

        {
            const size_t pz = patch_z_count, mz = merge_z_count;
            const size_t nm = normal_merge_x_count, ss = ss_count;

            // Same as pre-round: patch gauge ancillas are full-weight here.
            std::vector<uint32_t> all_z_post = z_ancillas;
            all_z_post.insert(all_z_post.end(), patch_z_gauge.begin(), patch_z_gauge.end());

            std::vector<uint32_t> all_x_post = x_ancillas;
            all_x_post.insert(all_x_post.end(), patch_x_gauge.begin(), patch_x_gauge.end());
            all_x_post.insert(all_x_post.end(), all_seam_x.begin(), all_seam_x.end());

            circuit_.safe_append_u("TICK", {}, {});
            circuit_.safe_append_u("R", all_z_post, {});
            circuit_.safe_append_u("RX", all_x_post, {});

            emit_post_layer(pre_cx_layer1, pre_seam_cx1, {});
            emit_post_layer(pre_cx_layer2, pre_seam_cx2, {});
            emit_post_layer(pre_cx_layer3, pre_seam_cx3, {});
            emit_post_layer(pre_cx_layer4, pre_seam_cx4, {});

            circuit_.safe_append_u("TICK", {}, {});
            circuit_.safe_append_u("M", all_z_post, {});
            circuit_.safe_append_u("MX", all_x_post, {});

            // Use correct last-merge-round measurement totals for has_gauges vs no-gauges.
            // For has_gauges: z_anc in half2 is at -(post_total + merge_h2 - i) from end of post-round.
            // For no-gauges: z_anc is at -(post_total + ss_merge_total - i).
            // Skip boundary-weight-changing ancillas: their post-round measurement is weight-4 while
            // the last merge measurement was weight-3, making the transition detector non-deterministic.
            for (size_t i = 0; i < patch_z_count; i++) {
                const auto& q = qubits_[z_ancillas[i]];
                if (boundary_changing_z_set.count(z_ancillas[i])) continue;
                int32_t curr = -(int32_t)(post_total - i);
                int32_t prev = has_gauges
                    ? -(int32_t)(post_total + merge_h2 - i)
                    : -(int32_t)(post_total + ss_merge_total - i);
                circuit_.safe_append_u("DETECTOR", {drec(curr), drec(prev)},
                    {q.x, q.y, static_cast<double>(final_round)});
            }
            for (size_t i = 0; i < patch_x_count; i++) {
                const auto& q = qubits_[x_ancillas[i]];
                if (boundary_changing_x_set.count(x_ancillas[i])) continue;
                int32_t curr = has_gauges
                    ? -(int32_t)(post_total - patch_z_count - patch_zg_count - i)
                    : -(int32_t)(post_total - pz - i);
                int32_t prev = has_gauges
                    ? -(int32_t)(post_total + patch_x_count + normal_merge_x_count - i)
                    : -(int32_t)(post_total + ss_merge_total - pz - mz - i);
                circuit_.safe_append_u("DETECTOR", {drec(curr), drec(prev)},
                    {q.x, q.y, static_cast<double>(final_round)});
            }
            for (size_t m = 0; m < nm; m++) {
                size_t j = normal_seam_indices[m];
                const auto& q = qubits_[seam_a_x_ancillas[j]];
                int32_t curr_a = -(int32_t)(seam_count - j);
                int32_t curr_b = -(int32_t)(seam_b_count - j);
                int32_t prev_nmx = -(int32_t)(post_total + 2 * ss + nm - m);
                circuit_.safe_append_u("DETECTOR",
                    {drec(curr_a), drec(curr_b), drec(prev_nmx)},
                    {q.x, q.y, static_cast<double>(final_round)});
            }

            // Orphan seam positions (merge-X suppressed at this y): no valid post-merge
            // detector. The seam X stabilizer anti-commutes with the adjacent merge-round
            // Z stabilizer because the excluded data qubit removes one of the two qubit
            // overlaps that normally cancel — leaving a net anti-commutation. Skipped.
        }

        circuit_.safe_append_u("TICK", {}, {});
    }

    // ====== Final measurement of all data qubits ======
    circuit_.safe_append_u("M", all_data, {});

    // Final data detectors: for each Z-stabilizer ancilla in the patches,
    // combine data measurements with the last ancilla measurement
    uint32_t data_round = merge_rounds_ + 2;

    // Helper: find Z ancilla measurement offset from end of post-round.
    // Searches z_ancillas first, then patch_z_gauge (appended after z_ancillas in M()).
    auto find_z_anc_offset = [&](uint32_t ancilla_idx) -> int32_t {
        for (size_t k = 0; k < z_ancillas.size(); k++) {
            if (z_ancillas[k] == ancilla_idx)
                return (int32_t)(all_data.size() + post_total - k);
        }
        // Gauge ancillas are intentionally excluded: their post-round measurement is
        // weight-4 while the last merge half was weight-3, making cross-boundary
        // detectors non-deterministic.
        return -1; // not found
    };

    // Helper: add data qubit final measurement reference by qubit index.
    auto add_data_meas = [&](std::vector<uint32_t>& tgts, uint32_t qubit_idx) {
        for (size_t k = 0; k < all_data.size(); k++) {
            if (all_data[k] == qubit_idx) {
                tgts.push_back(drec(-(int32_t)(all_data.size() - k)));
                return;
            }
        }
    };

    // Set of gauge ancilla indices — their final data detectors are skipped because
    // the last post-round measurement is weight-4 while the last merge measurement
    // was weight-3, making a cross-boundary detector non-deterministic.
    std::unordered_set<uint32_t> gauge_ancilla_set(
        z_gauge_ancillas.begin(), z_gauge_ancillas.end());
    // Also exclude boundary-weight-changing Z ancillas from final data detectors
    gauge_ancilla_set.insert(boundary_changing_z_set.begin(), boundary_changing_z_set.end());

    for (const auto& stab : patch_a_stabilizers_) {
        if (stab.is_x_type) continue;
        if (gauge_ancilla_set.count(stab.ancilla)) continue;
        const auto& aq = qubits_[stab.ancilla];

        std::vector<uint32_t> det_targets;
        for (uint32_t d : stab.data_qubits) add_data_meas(det_targets, d);
        int32_t off = find_z_anc_offset(stab.ancilla);
        if (off >= 0) det_targets.push_back(drec(-off));
        circuit_.safe_append_u("DETECTOR", det_targets, {aq.x, aq.y, static_cast<double>(data_round)});
    }

    for (const auto& stab : patch_b_stabilizers_) {
        if (stab.is_x_type) continue;
        if (gauge_ancilla_set.count(stab.ancilla)) continue;
        const auto& aq = qubits_[stab.ancilla];

        std::vector<uint32_t> det_targets;
        for (uint32_t d : stab.data_qubits) add_data_meas(det_targets, d);
        int32_t off = find_z_anc_offset(stab.ancilla);
        if (off >= 0) det_targets.push_back(drec(-off));
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

    // Helper: generate a POLYGON pragma for one stabilizer.
    // include_removed=true also adds SS data qubits (for pre/post rounds).
    // exclude_merge_excluded=true drops merge_excluded_data_qubits_ (for merge rounds).
    // Returns empty string if the resulting polygon has no vertices.
    auto make_polygon = [this](const DStabilizer& stab, bool is_merge, bool include_removed = false, bool exclude_merge_excluded = false) -> std::string {
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
            if (exclude_merge_excluded && merge_excluded_data_qubits_.count(q)) continue;
            int dx = static_cast<int>(std::round((qubits_[q].x - ax) * 2));
            int dy = static_cast<int>(std::round((qubits_[q].y - ay) * 2));
            pos_map[{dx, dy}] = q;
        }
        if (include_removed) {
            for (uint32_t ridx : removed_data_qubits_) {
                const auto& dq = qubits_[ridx];
                int dx = static_cast<int>(std::round((dq.x - ax) * 2));
                int dy = static_cast<int>(std::round((dq.y - ay) * 2));
                if (std::abs(dx) != 1 || std::abs(dy) != 1) continue;
                // Patch-filter: seam ancillas only couple to data on their own patch side.
                // SEAM_A → PATCH_A only; SEAM_B → PATCH_B only.
                if (stab.patch == DPatchID::SEAM_A && dq.patch != DPatchID::PATCH_A) continue;
                if (stab.patch == DPatchID::SEAM_B && dq.patch != DPatchID::PATCH_B) continue;
                // Patch ancillas only include removed qubits on their own patch.
                if (stab.patch == DPatchID::PATCH_A && dq.patch != DPatchID::PATCH_A) continue;
                if (stab.patch == DPatchID::PATCH_B && dq.patch != DPatchID::PATCH_B) continue;
                pos_map[{dx, dy}] = ridx;
            }
        }

        if (pos_map.empty()) return "";

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

    // Build pragma blocks:
    //   patch_pragmas_full = full-weight polygons for pre/post rounds (SS qubits active)
    //   patch_pragmas      = reduced-weight polygons for merge rounds (SS qubits removed)
    std::string patch_pragmas_full;
    std::string patch_pragmas;
    for (const auto& stab : patch_a_stabilizers_) {
        patch_pragmas_full += make_polygon(stab, false, true, false) + "\n";
        if (!merge_suppressed_ancillas_.count(stab.ancilla))
            patch_pragmas += make_polygon(stab, false, false, true) + "\n";
    }
    for (const auto& stab : patch_b_stabilizers_) {
        patch_pragmas_full += make_polygon(stab, false, true, false) + "\n";
        if (!merge_suppressed_ancillas_.count(stab.ancilla))
            patch_pragmas += make_polygon(stab, false, false, true) + "\n";
    }

    // Seam stabilizers: full-weight for pre/post rounds, reduced for merge rounds
    std::string seam_pragmas_full;
    std::string seam_pragmas;
    for (const auto& stab : seam_a_stabilizers_) {
        seam_pragmas_full += make_polygon(stab, false, true) + "\n";
        seam_pragmas      += make_polygon(stab, false, false) + "\n";
    }
    for (const auto& stab : seam_b_stabilizers_) {
        seam_pragmas_full += make_polygon(stab, false, true) + "\n";
        seam_pragmas      += make_polygon(stab, false, false) + "\n";
    }

    std::string merge_pragmas;
    for (const auto& stab : merge_stabilizers_) {
        if (merge_suppressed_ancillas_.count(stab.ancilla)) continue;
        merge_pragmas += make_polygon(stab, true, false, true) + "\n";
    }

    // Remote CNOT edges: one 2-qubit polygon per (merge ancilla, Patch B data) pair
    std::string remote_cx_pragmas;
    for (const auto& stab : merge_stabilizers_) {
        if (merge_suppressed_ancillas_.count(stab.ancilla)) continue;
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
    // marks the start of a new reset group.
    //
    // When superstabilizers are present (has_gauges), each merge round is split into
    // two half-rounds, so halves_per_merge=2; otherwise halves_per_merge=1.
    //   Group 0                               = data reset (no pragmas)
    //   Group 1                               = pre-merge round (initial pragmas already shown)
    //   Groups 2 .. 1 + halves*merge_rounds   = merge rounds (patch + merge pragmas)
    //   Group  2 + halves*merge_rounds         = post-merge round (patch + seam pragmas)
    // halves_per_merge matches the circuit: 2 only when gauge ancillas are present
    // (interior SS qubits split each merge round into two half-rounds).
    // Boundary SS qubits suppress ancillas but do not create gauges, so halves=1.
    const int halves_per_merge = has_gauges_ ? 2 : 1;
    const int merge_half_groups = halves_per_merge * static_cast<int>(merge_rounds_);

    std::ostringstream out;
    int reset_group = -1;
    bool prev_was_tick = false;
    bool in_qubit_coords = true;

    for (size_t i = 0; i < lines.size(); i++) {
        const auto& l = lines[i];

        // After QUBIT_COORDS block ends, insert initial full-weight patch + seam pragmas
        if (in_qubit_coords &&
            (l.size() < 12 || l.substr(0, 12) != "QUBIT_COORDS")) {
            in_qubit_coords = false;
            out << patch_pragmas_full << seam_pragmas_full;
        }

        // Detect start of a reset group
        bool is_reset = (l.size() >= 2 && l[0] == 'R' && l[1] == ' ') ||
                         (l.size() >= 3 && l[0] == 'R' && l[1] == 'X' && l[2] == ' ');
        if (is_reset && prev_was_tick) {
            reset_group++;
            if (reset_group >= 2 && reset_group <= 1 + merge_half_groups) {
                // Before any merge half-round: reduced-weight patch + merge + remote CX edges (no seam)
                out << patch_pragmas << merge_pragmas << remote_cx_pragmas;
            } else if (reset_group == 2 + merge_half_groups) {
                // Before post-merge round: full-weight patch + seam (SS qubits active again)
                out << patch_pragmas_full << seam_pragmas_full;
            }
        }

        prev_was_tick = (l == "TICK");
        out << l << "\n";
    }

    return out.str();
}

} // namespace bucket_sim
