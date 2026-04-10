#pragma once

#include "config.hpp"
#include "lattice_surgery.hpp"
#include <stim.h>
#include <pymatching/sparse_blossom/matcher/mwpm.h>
#include <cstdint>
#include <chrono>

namespace bucket_sim {

struct SimulationResults {
    uint64_t total_shots;
    uint64_t num_errors;
    double logical_error_rate;
    double runtime_seconds;
    double shots_per_second;
};

struct RankStats {
    int rank;
    uint64_t shots;
    double runtime_seconds;
    double shots_per_second;
};

// Result of a distillation calculation (paper Equations 3-6)
struct DistillationResult {
    DistillationProtocol protocol = DistillationProtocol::NONE;
    double output_fidelity;      // F_distill after k rounds
    double success_probability;  // P_s for the distillation attempt
    uint32_t raw_pairs_consumed; // Number of raw EPR pairs needed per output
    double distillation_time;    // Time required for distillation (seconds)
};

// Summary of computed noise parameters for a distributed circuit run
struct NoiseSummary {
    // Remote CNOT noise
    DistillationProtocol effective_distillation_protocol = DistillationProtocol::NONE;
    bool accurate_rcx_enabled  = false;
    bool distillation_feasible = false;
    bool distillation_qubit_overhead_feasible = true;
    uint32_t distillation_backup_batches = 1;
    double raw_channel_error   = 0.0;
    double effective_channel_error = 0.0;
    double distilled_fidelity   = 0.0;  // EPR fidelity after distillation
    double remote_cnot_error    = 0.0;  // p_cnot from Eq. 1
    double distilled_error      = 0.0;
    double distillation_success_probability = 1.0;
    double probability_all_distillation_fail = 0.0;
    uint32_t raw_pairs_per_distilled = 0; // raw EPR pairs consumed per distilled pair
    uint32_t remote_cnots_per_cycle  = 0; // remote CX pairs in one merge round
    uint32_t epr_pairs_per_round     = 0; // raw_pairs_per_distilled * remote_cnots_per_cycle
    uint32_t distillation_epr_slots_required = 0; // raw EPR slots across all remote CX and backup batches
    uint32_t distillation_qubits_required = 0; // 2 qubits per raw EPR slot
    uint32_t monolithic_equivalent_qubits = 0; // qubits needed to place the same LS experiment on one chip

    // Idling noise
    double distillation_time_ns = 0.0;  // t_dist in nanoseconds
    double idling_time_us       = 0.0;  // t_idle = max(t_meas, t_dist) in microseconds
    double p_X = 0.0;
    double p_Y = 0.0;
    double p_Z = 0.0;
    bool timing_constraint_satisfied = true;
};

struct RemoteCnotNoiseModel {
    DistillationProtocol effective_distillation_protocol = DistillationProtocol::NONE;
    double raw_channel_error = 0.0;
    double effective_channel_error = 0.0;
    double distilled_error = 0.0;
    double distilled_fidelity = 0.0;
    double distillation_success_probability = 1.0;
    double probability_all_distillation_fail = 0.0;
    double remote_cnot_error = 0.0;
    double distillation_time = 0.0;
    bool distillation_feasible = false;
    bool distillation_qubit_overhead_feasible = true;
    uint32_t distillation_backup_batches = 1;
    uint32_t raw_pairs_per_distilled = 1;
    uint32_t distillation_epr_slots_required = 0;
    uint32_t distillation_qubits_required = 0;
    uint32_t monolithic_equivalent_qubits = 0;
};

class SurfaceCodeSimulator {
private:
    Config config_;
    int mpi_rank_;
    int mpi_size_;
    uint64_t shots_this_rank_;

    stim::Circuit circuit_;
    pm::Mwpm mwpm_decoder_;

    uint64_t local_errors_;
    double local_runtime_;
    std::string annotated_circuit_str_;  // Pragma-annotated stim text (if available)

    NoiseSummary noise_summary_;
    uint32_t distributed_remote_cnots_per_cycle_override_ = 0;

    // Bucket mode statistics
    uint32_t num_sampled_buckets_;
    double sampled_probability_mass_;
    double bias_bound_;
    double statistical_error_;
    bool has_bucket_stats_;

    std::chrono::high_resolution_clock::time_point start_time_;
    std::chrono::high_resolution_clock::time_point end_time_;

    void initialize_circuit();
    void initialize_lattice_surgery_circuit();
    void initialize_distributed_lattice_surgery_circuit();
    void initialize_decoder();
    void inject_interconnect_noise();
    void inject_entanglement_idling_noise();
    uint32_t count_remote_cnots_in_cycle();

    // Distillation calculations (paper Equations 3-6)
    DistillationResult compute_distillation(DistillationProtocol protocol, uint32_t rounds) const;
    RemoteCnotNoiseModel evaluate_remote_cnot_noise_model(DistillationProtocol protocol,
                                                          uint32_t rounds,
                                                          uint32_t num_remote_cnots) const;
    DistillationProtocol select_radar_distillation_protocol(uint32_t rounds,
                                                            uint32_t num_remote_cnots) const;
    static double compute_2to1_fidelity(double F1, double F2);
    static double compute_2to1_success_prob(double F1, double F2);
    static double compute_3to1_fidelity(double F1, double F2);
    static double compute_3to1_success_prob(double F1, double F2);
    double compute_raw_channel_error() const;
    RemoteCnotNoiseModel compute_remote_cnot_noise_model(uint32_t num_remote_cnots) const;
    uint32_t compute_monolithic_equivalent_qubits() const;

    // Remote CNOT error (paper Equation 1)
    double compute_cnot_error_from_fidelity(double fidelity) const;
    static double compute_remote_cnot_error_from_effective_channel(double effective_channel_error,
                                                                   double physical_error,
                                                                   bool accurate_rcx);

    // Timing calculations
    double compute_distillation_time(uint32_t num_remote_cnots, const DistillationResult& distill) const;
    bool check_timing_constraint(double distillation_time) const;

    void run_monte_carlo();
    void run_bucket();

public:
    SurfaceCodeSimulator(const Config& config, int rank, int size, bool skip_decoder = false);

    void run();

    uint64_t get_local_errors() const { return local_errors_; }
    double get_local_runtime() const { return local_runtime_; }
    uint64_t get_local_shots() const { return shots_this_rank_; }

    bool has_bucket_stats() const { return has_bucket_stats_; }
    uint32_t get_num_sampled_buckets() const { return num_sampled_buckets_; }
    double get_sampled_probability_mass() const { return sampled_probability_mass_; }
    double get_bias_bound() const { return bias_bound_; }
    double get_statistical_error() const { return statistical_error_; }

    RankStats get_rank_stats() const;

    const NoiseSummary& get_noise_summary() const { return noise_summary_; }

    // Get the generated circuit (for visualization/debugging)
    const stim::Circuit& get_circuit() const { return circuit_; }
    const std::string& get_annotated_circuit_str() const { return annotated_circuit_str_; }
};

} // namespace bucket_sim
