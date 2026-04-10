#include "simulator.hpp"

#include <cmath>
#include <cstdint>
#include <iostream>
#include <regex>
#include <string>

using bucket_sim::Config;
using bucket_sim::DistillationProtocol;
using bucket_sim::MergeType;
using bucket_sim::SurfaceCodeSimulator;

namespace {

bool approx_equal(double a, double b, double eps = 1e-9) {
    return std::abs(a - b) <= eps;
}

void require(bool cond, const std::string &msg) {
    if (!cond) {
        std::cerr << "[FAIL] " << msg << std::endl;
        std::exit(1);
    }
}

Config make_base_config() {
    Config cfg;
    cfg.code_distance = 3;
    cfg.rounds = 3;
    cfg.total_shots = 100;
    cfg.distributed = true;
    cfg.merge_type = MergeType::XX_MERGE_DISTRIBUTED;
    cfg.merge_rounds = 3;
    cfg.physical_error = 1e-3;
    cfg.measurement_error = 1e-3;
    cfg.reset_error = 1e-3;
    cfg.interconnect_error = 0.0;
    cfg.channel_depolarization_error = 0.02;
    cfg.entanglement_rate = 1e9;
    cfg.measurement_delay = 1e-3;
    cfg.raw_epr_fidelity = 0.99;
    cfg.distillation_backup_batches = 1;
    return cfg;
}

void test_no_distillation() {
    Config cfg = make_base_config();
    cfg.distillation_protocol = DistillationProtocol::NONE;
    cfg.distillation_rounds = 0;

    SurfaceCodeSimulator legacy(cfg, 0, 1, true);
    const auto &legacy_noise = legacy.get_noise_summary();
    require(!legacy_noise.distillation_feasible, "distillation should be infeasible when disabled");
    require(legacy_noise.effective_distillation_protocol == DistillationProtocol::NONE,
            "effective protocol should remain NONE when distillation is disabled");
    require(approx_equal(legacy_noise.effective_channel_error, cfg.channel_depolarization_error),
            "effective channel error should equal raw channel error without distillation");
    require(approx_equal(legacy_noise.remote_cnot_error, cfg.channel_depolarization_error),
            "remote CNOT error should equal effective channel error when accurate_rcx is off");

    cfg.accurate_rcx = true;
    SurfaceCodeSimulator accurate(cfg, 0, 1, true);
    const auto &accurate_noise = accurate.get_noise_summary();
    double expected_remote = 1.0 - (1.0 - cfg.channel_depolarization_error) *
                                     (1.0 - cfg.physical_error) * (1.0 - cfg.physical_error);
    require(approx_equal(accurate_noise.remote_cnot_error, expected_remote),
            "accurate_rcx should fold in the two local physical-error factors");
}

void test_distillation_feasible(uint32_t backups) {
    Config cfg = make_base_config();
    cfg.distillation_protocol = DistillationProtocol::PUMPING_2TO1;
    cfg.distillation_rounds = 1;
    cfg.distillation_backup_batches = backups;

    SurfaceCodeSimulator simulator(cfg, 0, 1, true);
    const auto &noise = simulator.get_noise_summary();
    require(noise.distillation_feasible, "distillation should be feasible in the high-rate test case");
    require(noise.effective_distillation_protocol == DistillationProtocol::PUMPING_2TO1,
            "effective protocol should match the pinned protocol in the feasible test case");

    double distilled_error = 1.0 - noise.distilled_fidelity;
    double p_all_fail = std::pow(1.0 - noise.distillation_success_probability, static_cast<double>(backups));
    double expected_effective = (1.0 - p_all_fail) * distilled_error +
                                p_all_fail * cfg.channel_depolarization_error;

    require(approx_equal(noise.probability_all_distillation_fail, p_all_fail),
            "summary should report the all-fail probability implied by backup batches");
    require(approx_equal(noise.effective_channel_error, expected_effective),
            "effective channel error should match the distilled/raw mixture");
}

void test_pinned_distillation_can_extend_round() {
    Config cfg = make_base_config();
    cfg.distillation_protocol = DistillationProtocol::PUMPING_2TO1;
    cfg.distillation_rounds = 1;
    cfg.entanglement_rate = 1e3;
    cfg.measurement_delay = 1e-9;

    SurfaceCodeSimulator simulator(cfg, 0, 1, true);
    const auto &noise = simulator.get_noise_summary();
    require(!noise.distillation_feasible, "distillation should be timing-infeasible in the low-rate test case");
    require(noise.effective_distillation_protocol == DistillationProtocol::PUMPING_2TO1,
            "pinned protocol should still be used when it extends the round");
    require(noise.effective_channel_error < cfg.channel_depolarization_error,
            "pinned distillation should still improve the channel even when timing is infeasible");
    require(noise.idling_time_us > cfg.measurement_delay * 1e6,
            "timing-infeasible distillation should extend the syndrome round");
}

void test_radar_can_choose_raw_pairs() {
    Config cfg = make_base_config();
    cfg.distillation_protocol = DistillationProtocol::RADAR;
    cfg.distillation_rounds = 1;
    cfg.entanglement_rate = 1e6;
    cfg.measurement_delay = 7e-6;

    SurfaceCodeSimulator simulator(cfg, 0, 1, true);
    const auto &noise = simulator.get_noise_summary();
    require(noise.effective_distillation_protocol == DistillationProtocol::NONE,
            "radar should choose raw EPR pairs when distillation would be too costly");
    require(approx_equal(noise.effective_channel_error, cfg.channel_depolarization_error),
            "radar should fall back to the raw channel when it selects no distillation");
}

void test_pinned_parameters_report_timing_feasibility() {
    Config cfg = make_base_config();
    cfg.distillation_protocol = DistillationProtocol::PUMPING_3TO1;
    cfg.distillation_rounds = 2;
    cfg.distillation_backup_batches = 4;
    cfg.entanglement_rate = 1e12;
    cfg.measurement_delay = 1e-6;

    SurfaceCodeSimulator feasible(cfg, 0, 1, true);
    const auto &feasible_noise = feasible.get_noise_summary();
    require(feasible_noise.effective_distillation_protocol == DistillationProtocol::PUMPING_3TO1,
            "pinned protocol should be used in the feasibility test");
    require(feasible_noise.distillation_backup_batches == 4,
            "summary should preserve the configured backup batch count");
    require(feasible_noise.raw_pairs_per_distilled == 6,
            "two rounds of 3to1 pumping should consume 6 raw pairs per distilled pair");
    require(feasible_noise.distillation_feasible,
            "high ENR should make the pinned (protocol, k, m) timing-feasible");

    cfg.entanglement_rate = 1e6;
    SurfaceCodeSimulator infeasible(cfg, 0, 1, true);
    const auto &infeasible_noise = infeasible.get_noise_summary();
    require(!infeasible_noise.distillation_feasible,
            "low ENR should make the same pinned (protocol, k, m) timing-infeasible");
}

void test_distillation_overhead_guard_against_monolithic() {
    Config cfg = make_base_config();
    cfg.distillation_protocol = DistillationProtocol::RECURRENCE_3TO1;
    cfg.distillation_rounds = 2;
    cfg.distillation_backup_batches = 5;
    cfg.entanglement_rate = 1e12;
    cfg.measurement_delay = 1e-3;

    SurfaceCodeSimulator simulator(cfg, 0, 1, true);
    const auto &noise = simulator.get_noise_summary();
    require(noise.effective_distillation_protocol == DistillationProtocol::RECURRENCE_3TO1,
            "overhead test should use the pinned recurrence protocol");
    require(noise.raw_pairs_per_distilled == 9,
            "two rounds of 3to1 recurrence should consume 9 raw pairs per distilled pair");
    require(noise.remote_cnots_per_cycle == 5,
            "d=3 distributed layout should have 5 remote CX operations per merge round");
    require(noise.distillation_epr_slots_required == 225,
            "recurrence overhead should scale as remote_cnots * 3^k * m");
    require(noise.distillation_qubits_required == 450,
            "recurrence overhead should use two qubits per concurrent EPR slot");
    require(noise.monolithic_equivalent_qubits == 41,
            "monolithic equivalent for d=3 XX lattice surgery should require 41 qubits");
    require(!noise.distillation_qubit_overhead_feasible,
            "distillation overhead larger than the monolithic footprint should be rejected");
}

void test_pumping_overhead_is_constant_width() {
    Config cfg = make_base_config();
    cfg.distillation_protocol = DistillationProtocol::PUMPING_2TO1;
    cfg.distillation_rounds = 4;
    cfg.distillation_backup_batches = 3;
    cfg.entanglement_rate = 1e12;
    cfg.measurement_delay = 1e-3;

    SurfaceCodeSimulator simulator(cfg, 0, 1, true);
    const auto &noise = simulator.get_noise_summary();
    require(noise.effective_distillation_protocol == DistillationProtocol::PUMPING_2TO1,
            "constant-width overhead test should use pumping");
    require(noise.raw_pairs_per_distilled == 8,
            "four rounds of 2to1 pumping should still consume 8 raw pairs in expectation accounting");
    require(noise.remote_cnots_per_cycle == 5,
            "d=3 distributed layout should have 5 remote CX operations per merge round");
    require(noise.distillation_epr_slots_required == 30,
            "2to1 pumping overhead should use a constant 2-pair width per chain");
    require(noise.distillation_qubits_required == 60,
            "2to1 pumping should need four qubits per chain per backup batch");
}

void test_remote_cx_injection_is_uniform() {
    Config cfg = make_base_config();
    cfg.accurate_rcx = true;
    cfg.distillation_protocol = DistillationProtocol::PUMPING_2TO1;
    cfg.distillation_rounds = 1;
    cfg.distillation_backup_batches = 3;

    SurfaceCodeSimulator simulator(cfg, 0, 1, true);
    const auto &noise = simulator.get_noise_summary();
    const std::string &annotated = simulator.get_annotated_circuit_str();
    require(!annotated.empty(), "annotated circuit should be available after remote-noise injection");

    std::regex depol_line(R"(DEPOLARIZE2\(([0-9eE+\-.]+)\).+# interconnect)");
    std::smatch match;
    std::string::const_iterator search_start = annotated.cbegin();
    size_t seen = 0;
    while (std::regex_search(search_start, annotated.cend(), match, depol_line)) {
        double injected = std::stod(match[1].str());
        require(approx_equal(injected, noise.remote_cnot_error, 1e-6),
                "all remote CX injections should use the same config-level probability");
        search_start = match.suffix().first;
        seen++;
    }
    require(seen > 0, "expected at least one injected remote-CNOT noise line");
}

}  // namespace

int main() {
    test_no_distillation();
    test_distillation_feasible(1);
    test_distillation_feasible(3);
    test_pinned_distillation_can_extend_round();
    test_radar_can_choose_raw_pairs();
    test_pinned_parameters_report_timing_feasibility();
    test_distillation_overhead_guard_against_monolithic();
    test_pumping_overhead_is_constant_width();
    test_remote_cx_injection_is_uniform();
    std::cout << "[PASS] accurate_rcx_formulas" << std::endl;
    return 0;
}
