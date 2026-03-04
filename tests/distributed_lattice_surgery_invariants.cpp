#include "config.hpp"
#include "distributed_lattice_surgery.hpp"

#include <stim.h>

#include <cmath>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using bucket_sim::Config;
using bucket_sim::DPatchID;
using bucket_sim::DQubitType;
using bucket_sim::DistributedLatticeSurgeryCircuit;
using bucket_sim::MergeType;

namespace {

bool approx_equal(double a, double b, double eps = 1e-9) {
    return std::abs(a - b) < eps;
}

bool crosses_cut(double x1, double x2, double split_x) {
    return (x1 < split_x && x2 >= split_x) || (x1 >= split_x && x2 < split_x);
}

void require(bool cond, const std::string &msg) {
    if (!cond) {
        std::cerr << "[FAIL] " << msg << std::endl;
        std::exit(1);
    }
}

// Throws (via stim) if the circuit has any non-deterministic (gauge) detectors.
void check_det_determinism(const stim::Circuit &circuit, const std::string &label) {
    try {
        stim::circuit_to_dem(circuit, stim::DemOptions{
            .decompose_errors = false,
            .flatten_loops = true,
            .allow_gauge_detectors = false,
        });
    } catch (const std::exception &e) {
        std::cerr << "[FAIL] non-deterministic detector in " << label << ": " << e.what() << std::endl;
        std::exit(1);
    }
}

void run_case(uint32_t d) {
    Config cfg;
    cfg.code_distance = d;
    cfg.rounds = d;
    cfg.merge_type = MergeType::XX_MERGE_DISTRIBUTED;
    cfg.merge_rounds = d;
    cfg.distributed = true;
    cfg.interconnect_error = 0.01;
    cfg.total_shots = 1000;

    DistributedLatticeSurgeryCircuit ls(cfg);
    const auto &qubits = ls.qubits();
    const auto &patch_a_stabs = ls.patch_a_stabs();
    const auto &patch_b_stabs = ls.patch_b_stabs();
    const auto &seam_a_stabs = ls.seam_a_stabs();
    const auto &seam_b_stabs = ls.seam_b_stabs();
    const auto &merge_stabs = ls.merge_stabs();

    size_t data_count = 0;
    size_t merge_data_count = 0;
    size_t merge_ancilla_count = 0;
    for (const auto &q : qubits) {
        if (q.type == DQubitType::DATA) {
            data_count++;
            if (q.patch == DPatchID::MERGE) {
                merge_data_count++;
            }
        }
        if (q.patch == DPatchID::MERGE && q.type != DQubitType::DATA) {
            merge_ancilla_count++;
            require(approx_equal(q.x, static_cast<double>(d)), "merge ancilla not on seam column x=d");
        }
    }

    require(data_count == static_cast<size_t>(2 * d * d), "unexpected total data qubit count");
    require(merge_data_count == 0, "merge data strip should not exist in distributed mode");
    require(merge_ancilla_count == static_cast<size_t>(d), "expected one seam ancilla column of size d");
    require(merge_stabs.size() == static_cast<size_t>(d), "expected exactly d merge stabilizers");

    // Each patch should have d²-1 stabilizers (interior + boundary)
    size_t expected_patch_stabs = static_cast<size_t>(d * d - 1);
    require(patch_a_stabs.size() + seam_a_stabs.size() == expected_patch_stabs,
            "Patch A total stabilizers (interior + seam) should be d²-1, got " +
            std::to_string(patch_a_stabs.size() + seam_a_stabs.size()));
    require(patch_b_stabs.size() + seam_b_stabs.size() == expected_patch_stabs,
            "Patch B total stabilizers (interior + seam) should be d²-1, got " +
            std::to_string(patch_b_stabs.size() + seam_b_stabs.size()));

    // Seam stabilizer count: (d-1)/2 per patch
    size_t expected_seam = static_cast<size_t>((d - 1) / 2);
    require(seam_a_stabs.size() == expected_seam,
            "Seam A should have (d-1)/2 stabilizers, got " + std::to_string(seam_a_stabs.size()));
    require(seam_b_stabs.size() == expected_seam,
            "Seam B should have (d-1)/2 stabilizers, got " + std::to_string(seam_b_stabs.size()));

    // All seam stabilizers should be X-type and weight-2
    for (const auto &stab : seam_a_stabs) {
        require(stab.is_x_type, "seam A stabilizer should be X-type");
        require(stab.data_qubits.size() == 2, "seam A stabilizer should be weight-2");
        // Should only touch Patch A data
        for (uint32_t qidx : stab.data_qubits) {
            require(qubits[qidx].patch == DPatchID::PATCH_A, "seam A stab should only touch Patch A data");
        }
    }
    for (const auto &stab : seam_b_stabs) {
        require(stab.is_x_type, "seam B stabilizer should be X-type");
        require(stab.data_qubits.size() == 2, "seam B stabilizer should be weight-2");
        for (uint32_t qidx : stab.data_qubits) {
            require(qubits[qidx].patch == DPatchID::PATCH_B, "seam B stab should only touch Patch B data");
        }
    }

    for (const auto &stab : merge_stabs) {
        bool has_a = false;
        bool has_b = false;
        for (uint32_t qidx : stab.data_qubits) {
            if (qubits[qidx].patch == DPatchID::PATCH_A) has_a = true;
            if (qubits[qidx].patch == DPatchID::PATCH_B) has_b = true;
        }
        require(has_a && has_b, "every merge stabilizer must bridge both patches");
    }

    // Verify circuit is valid (generates, has merge CX pairs crossing boundary)
    auto circuit = ls.generate();
    check_det_determinism(circuit, "run_case d=" + std::to_string(d));
    auto coords = circuit.get_final_qubit_coords();
    const double split_x = static_cast<double>(d);

    std::istringstream iss(circuit.str());
    std::string line;
    uint64_t crossing_cx_pairs = 0;
    while (std::getline(iss, line)) {
        if (line.empty()) continue;
        std::istringstream ls_line(line);
        std::string gate;
        ls_line >> gate;
        if (gate == "CX" || gate == "CNOT") {
            uint32_t a, b;
            while (ls_line >> a >> b) {
                auto ita = coords.find(a);
                auto itb = coords.find(b);
                require(ita != coords.end() && itb != coords.end(), "CX pair missing coordinates");
                // Remote: one qubit on QPU A (x <= split_x), other on QPU B (x > split_x)
                bool crosses = (ita->second[0] <= split_x && itb->second[0] > split_x) ||
                               (ita->second[0] > split_x && itb->second[0] <= split_x);
                if (crosses) crossing_cx_pairs++;
            }
        }
    }

    require(crossing_cx_pairs > 0, "expected at least one remote CX crossing the cut");
    // Note: DEPOLARIZE2 noise is injected by the simulator (inject_interconnect_noise),
    // not by the circuit generator, so we don't check for it here.
}

// Count remote CX pairs (crossing the split at x=d) in a circuit.
uint64_t count_remote_cx(const stim::Circuit& circuit, double split_x) {
    auto coords = circuit.get_final_qubit_coords();
    std::istringstream iss(circuit.str());
    std::string line;
    uint64_t count = 0;
    while (std::getline(iss, line)) {
        if (line.empty()) continue;
        std::istringstream ls_line(line);
        std::string gate;
        ls_line >> gate;
        if (gate == "CX" || gate == "CNOT") {
            uint32_t a, b;
            while (ls_line >> a >> b) {
                auto ita = coords.find(a);
                auto itb = coords.find(b);
                if (ita == coords.end() || itb == coords.end()) continue;
                bool crosses = (ita->second[0] <= split_x && itb->second[0] > split_x) ||
                               (ita->second[0] > split_x && itb->second[0] <= split_x);
                if (crosses) count++;
            }
        }
    }
    return count;
}

void run_superstab_case(uint32_t d) {
    const double split_x = static_cast<double>(d);

    // Build the list of all X-type merge positions (y = 1,3,...,d-2) for full superstabilization.
    std::vector<uint32_t> all_x_ys;
    for (uint32_t y = 1; y <= d - 2; y += 2) all_x_ys.push_back(y);

    // Normal merge (no superstabilizer) for baseline remote-CX count.
    Config cfg_normal;
    cfg_normal.code_distance = d;
    cfg_normal.rounds = d;
    cfg_normal.merge_type = MergeType::XX_MERGE_DISTRIBUTED;
    cfg_normal.merge_rounds = d;
    cfg_normal.distributed = true;
    DistributedLatticeSurgeryCircuit ls_normal(cfg_normal);
    auto circ_normal = ls_normal.generate();
    check_det_determinism(circ_normal, "normal d=" + std::to_string(d));
    uint64_t remote_cx_normal = count_remote_cx(circ_normal, split_x);

    // Partial superstabilizer: only the first X position.
    if (!all_x_ys.empty()) {
        Config cfg_partial = cfg_normal;
        cfg_partial.superstab_ys = {all_x_ys[0]};
        DistributedLatticeSurgeryCircuit ls_partial(cfg_partial);
        auto circ_partial = ls_partial.generate();
        check_det_determinism(circ_partial, "partial superstab d=" + std::to_string(d));
        uint64_t remote_cx_partial = count_remote_cx(circ_partial, split_x);
        require(remote_cx_partial < remote_cx_normal,
                "partial superstab should reduce remote CX count (d=" + std::to_string(d) + ")");
    }

    // Full superstabilizer: all X positions.
    Config cfg_full = cfg_normal;
    cfg_full.superstab_ys = all_x_ys;
    DistributedLatticeSurgeryCircuit ls_full(cfg_full);
    auto circ_full = ls_full.generate();
    check_det_determinism(circ_full, "full superstab d=" + std::to_string(d));
    uint64_t remote_cx_full = count_remote_cx(circ_full, split_x);
    require(remote_cx_full < remote_cx_normal,
            "full superstab should reduce remote CX count (d=" + std::to_string(d) + ")");
    // Each superstabilized X position eliminates exactly 2 remote CX per round:
    //   - seam_a (at x=d-0.125, A side) only touches Patch A data → local CX
    //   - seam_b (at x=d+0.125, B side) only touches Patch B data → local CX
    //   vs. merge-X (at x=d, A side) which had 2 CX to Patch B data (remote)
    uint64_t expected_reduction = static_cast<uint64_t>(all_x_ys.size()) * 2 * d;
    require(remote_cx_normal - remote_cx_full == expected_reduction,
            "full superstab should eliminate 2 remote CX per superstab position per round "
            "(expected reduction " + std::to_string(expected_reduction) +
            ", got " + std::to_string(remote_cx_normal - remote_cx_full) +
            ", d=" + std::to_string(d) + ")");

    // No-superstab with explicit empty list should behave identically to normal.
    Config cfg_empty = cfg_normal;
    cfg_empty.superstab_ys = {};
    DistributedLatticeSurgeryCircuit ls_empty(cfg_empty);
    auto circ_empty = ls_empty.generate();
    check_det_determinism(circ_empty, "empty superstab_ys d=" + std::to_string(d));
    uint64_t remote_cx_empty = count_remote_cx(circ_empty, split_x);
    require(remote_cx_empty == remote_cx_normal,
            "empty superstab_ys should produce same circuit as no superstab (d=" + std::to_string(d) + ")");
}

}  // namespace

int main() {
    run_case(3);
    run_case(5);
    run_superstab_case(3);
    run_superstab_case(5);
    std::cout << "[PASS] distributed lattice surgery invariants" << std::endl;
    return 0;
}
