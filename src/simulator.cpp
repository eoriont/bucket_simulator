#include "simulator.hpp"
#include "decoder.hpp"
#include "dem_analyzer.hpp"
#include "bucket_sampler.hpp"
#include "lattice_surgery.hpp"
#include "distributed_lattice_surgery.hpp"
#include <stim/circuit/circuit.h>
#include <stim/simulators/error_analyzer.h>
#include <pymatching/sparse_blossom/driver/user_graph.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>

namespace bucket_sim {

constexpr uint64_t BATCH_SIZE = 1000000; // Process 1M shots per batch

SurfaceCodeSimulator::SurfaceCodeSimulator(const Config& config, int rank, int size, bool skip_decoder)
    : config_(config),
      mpi_rank_(rank),
      mpi_size_(size),
      local_errors_(0),
      local_runtime_(0.0),
      num_sampled_buckets_(0),
      sampled_probability_mass_(0.0),
      bias_bound_(0.0),
      statistical_error_(0.0),
      has_bucket_stats_(false)
{
    // Calculate shots for this rank
    shots_this_rank_ = config_.total_shots / mpi_size_;
    if (rank == mpi_size_ - 1) {
        // Last rank takes any remainder
        shots_this_rank_ += config_.total_shots % mpi_size_;
    }

    if (mpi_rank_ == 0) {
        std::cout << "Initializing simulator with configuration:" << std::endl;
        std::cout << "  Code distance: " << config_.code_distance << std::endl;
        std::cout << "  Rounds: " << config_.rounds << std::endl;
        std::cout << "  Physical error: " << config_.physical_error << std::endl;
        std::cout << "  Measurement error: " << config_.measurement_error << std::endl;
        std::cout << "  Reset error: " << config_.reset_error << std::endl;
        std::cout << "  Total shots: " << config_.total_shots << std::endl;
        std::cout << "  Code type: " << config_.code_type << std::endl;
        std::cout << "  MPI ranks: " << mpi_size_ << std::endl;
    }

    initialize_circuit();
    if (!skip_decoder) {
        initialize_decoder();
    }
}

void SurfaceCodeSimulator::initialize_circuit() {
    // Check if we're doing lattice surgery
    if (config_.merge_type == MergeType::XX_MERGE_DISTRIBUTED) {
        initialize_distributed_lattice_surgery_circuit();
        return;
    }
    if (config_.merge_type != MergeType::NONE) {
        initialize_lattice_surgery_circuit();
        return;
    }

    // Set up circuit generation parameters
    stim::CircuitGenParameters params(config_.rounds, config_.code_distance, config_.code_type);

    // Set error rates
    params.before_round_data_depolarization = config_.physical_error;
    params.after_clifford_depolarization = config_.physical_error;
    params.after_reset_flip_probability = config_.reset_error;
    params.before_measure_flip_probability = config_.measurement_error;

    // Generate surface code circuit
    circuit_ = stim::generate_surface_code_circuit(params).circuit;

    // Inject interconnect noise if distributed
    inject_interconnect_noise();

    // Inject entanglement-limited idling noise if distributed
    inject_entanglement_idling_noise();

    if (mpi_rank_ == 0) {
        std::cout << "Circuit generated:" << std::endl;
        std::cout << "  Qubits: " << circuit_.count_qubits() << std::endl;
        std::cout << "  Detectors: " << circuit_.count_detectors() << std::endl;
        std::cout << "  Observables: " << circuit_.count_observables() << std::endl;
        if (config_.distributed) {
            std::cout << "  Mode: Distributed QEC" << std::endl;
            std::cout << "  Interconnect Error: " << config_.interconnect_error << std::endl;
        }
    }
}

void SurfaceCodeSimulator::initialize_lattice_surgery_circuit() {
    if (mpi_rank_ == 0) {
        std::cout << "Initializing lattice surgery circuit:" << std::endl;
        std::cout << "  Merge type: " << (config_.merge_type == MergeType::XX_MERGE ? "XX" : "ZZ") << std::endl;
        std::cout << "  Code distance: " << config_.code_distance << std::endl;
        std::cout << "  Merge rounds: " << (config_.merge_rounds > 0 ? config_.merge_rounds : config_.code_distance) << std::endl;
        std::cout << "  Distributed: " << (config_.distributed ? "yes" : "no") << std::endl;
    }

    // Create lattice surgery circuit generator
    LatticeSurgeryCircuit ls_circuit(config_);

    // Generate the circuit
    circuit_ = ls_circuit.generate();

    if (mpi_rank_ == 0) {
        std::cout << "Lattice surgery circuit generated:" << std::endl;
        std::cout << "  Total qubits: " << circuit_.count_qubits() << std::endl;
        std::cout << "  Detectors: " << circuit_.count_detectors() << std::endl;
        std::cout << "  Observables: " << circuit_.count_observables() << std::endl;
        std::cout << "  Merge stabilizers cross QPU boundary: "
                  << (config_.distributed ? "yes (remote CNOTs required)" : "n/a") << std::endl;
    }
}

void SurfaceCodeSimulator::initialize_distributed_lattice_surgery_circuit() {
    if (mpi_rank_ == 0) {
        std::cout << "Initializing DISTRIBUTED lattice surgery circuit:" << std::endl;
        std::cout << "  Mode: Remote CNOTs only, no merge data qubits" << std::endl;
        std::cout << "  Code distance: " << config_.code_distance << std::endl;
        std::cout << "  Merge rounds: " << (config_.merge_rounds > 0 ? config_.merge_rounds : config_.code_distance) << std::endl;
    }

    DistributedLatticeSurgeryCircuit dls_circuit(config_);
    circuit_ = dls_circuit.generate();

    // Seed annotated string with pragma polygons so inject_interconnect_noise
    // preserves them when rewriting the circuit text.
    annotated_circuit_str_ = dls_circuit.annotated_stim_str();

    // Inject interconnect noise (remote merge CNOTs get fidelity-derived error,
    // local merge CNOTs get physical_error)
    inject_interconnect_noise();

    // Inject entanglement-limited idling noise if configured
    inject_entanglement_idling_noise();

    if (mpi_rank_ == 0) {
        std::cout << "Distributed lattice surgery circuit generated:" << std::endl;
        std::cout << "  Total qubits: " << circuit_.count_qubits() << std::endl;
        std::cout << "  Data qubits: " << dls_circuit.num_data_qubits() << std::endl;
        std::cout << "  Detectors: " << circuit_.count_detectors() << std::endl;
        std::cout << "  Observables: " << circuit_.count_observables() << std::endl;
    }
}

void SurfaceCodeSimulator::initialize_decoder() {
    // Convert circuit to detector error model
    auto dem = stim::ErrorAnalyzer::circuit_to_detector_error_model(
        circuit_,
        false,  // decompose_errors
        true,   // fold_loops
        false,  // allow_gauge_detectors
        0,      // approximate_disjoint_errors_threshold
        false,  // ignore_decomposition_failures
        false   // block_decomposition_from_introducing_remnant_edges
    );

    // Create MWPM decoder from detector error model
    size_t num_buckets = pm::NUM_DISTINCT_WEIGHTS;
    mwpm_decoder_ = pm::detector_error_model_to_mwpm(dem, num_buckets);

    if (mpi_rank_ == 0) {
        std::cout << "MWPM decoder initialized" << std::endl;
    }
}

void SurfaceCodeSimulator::inject_interconnect_noise() {
    if (!config_.distributed) {
        return;
    }

    // Compute distillation result using accurate formulas (paper Equations 3-6)
    DistillationResult distill = compute_distillation(config_.distillation_rounds);

    // Compute remote CNOT error using paper Equation 1
    double p_cnot = compute_cnot_error_from_fidelity(distill.output_fidelity);

    // If user specified interconnect_error > 0, add it to the distillation error
    // (models additional noise beyond EPR infidelity)
    double total_error = p_cnot;
    if (config_.interconnect_error > 0) {
        // Combine errors (assuming independent):
        // p_total = 1 - (1-p_cnot)(1-p_interconnect)
        total_error = 1.0 - (1.0 - p_cnot) * (1.0 - config_.interconnect_error);
    }

    if (mpi_rank_ == 0) {
        std::cout << "Entanglement distillation:" << std::endl;
        std::cout << "  Protocol: ";
        switch (config_.distillation_protocol) {
            case DistillationProtocol::PUMPING_2TO1: std::cout << "2→1 Pumping"; break;
            case DistillationProtocol::PUMPING_3TO1: std::cout << "3→1 Pumping"; break;
            case DistillationProtocol::RECURRENCE_2TO1: std::cout << "2→1 Recurrence"; break;
            case DistillationProtocol::RECURRENCE_3TO1: std::cout << "3→1 Recurrence"; break;
            default: std::cout << "None"; break;
        }
        std::cout << std::endl;
        std::cout << "  Distillation rounds: " << config_.distillation_rounds << std::endl;
        std::cout << "  Raw EPR fidelity: " << config_.raw_epr_fidelity << std::endl;
        std::cout << "  Distilled EPR fidelity: " << distill.output_fidelity << std::endl;
        std::cout << "  Distillation success prob: " << distill.success_probability << std::endl;
        std::cout << "  Raw pairs per distilled: " << distill.raw_pairs_consumed << std::endl;
        std::cout << "  Remote CNOT error (Eq.1): " << total_error << std::endl;
    }

    if (total_error <= 0) {
        return;  // No error to inject
    }

    // Get qubit coordinate mapping
    auto coords = circuit_.get_final_qubit_coords();

    // Determine split position (boundary between QPU A and QPU B)
    // QPU A owns x <= split_x, QPU B owns x > split_x.
    // Merge ancillas at x = split_x are on QPU A.
    double split_x = config_.code_distance;

    // Local gate error for merge CNOTs that don't cross the boundary
    double local_error = config_.physical_error;

    // Use pragma-annotated string if available, otherwise raw circuit
    std::string circuit_str = annotated_circuit_str_.empty()
        ? circuit_.str() : annotated_circuit_str_;
    std::istringstream iss(circuit_str);
    std::ostringstream oss;

    std::string line;
    while (std::getline(iss, line)) {
        // Check if this line is a CX/CNOT instruction
        if (line.find("CX ") == 0 || line.find("CNOT ") == 0) {
            // Parse qubit pairs from the line
            std::istringstream line_stream(line);
            std::string gate;
            line_stream >> gate;  // Read CX or CNOT

            std::vector<std::pair<uint32_t, uint32_t>> all_pairs;
            std::vector<std::pair<uint32_t, uint32_t>> remote_pairs;
            uint32_t ctrl, tgt;
            while (line_stream >> ctrl >> tgt) {
                all_pairs.push_back({ctrl, tgt});
                if (coords.count(ctrl) && coords.count(tgt)) {
                    double ctrl_x = coords[ctrl][0];
                    double tgt_x = coords[tgt][0];

                    // Remote: one qubit on QPU A (x <= split_x), other on QPU B (x > split_x)
                    bool crosses = (ctrl_x <= split_x && tgt_x > split_x) ||
                                  (ctrl_x > split_x && tgt_x <= split_x);
                    if (crosses) {
                        remote_pairs.push_back({ctrl, tgt});
                    }
                }
            }

            oss << line << "\n";

            // Baseline DEPOLARIZE2(physical_error) on ALL CX pairs
            if (!all_pairs.empty() && local_error > 0) {
                oss << "DEPOLARIZE2(" << local_error << ")";
                for (const auto& pair : all_pairs) {
                    oss << " " << pair.first << " " << pair.second;
                }
                oss << "\n";
            }

            // Additional DEPOLARIZE2(remote_error) on remote pairs only
            if (!remote_pairs.empty()) {
                oss << "DEPOLARIZE2(" << total_error << ")";
                for (const auto& pair : remote_pairs) {
                    oss << " " << pair.first << " " << pair.second;
                }
                oss << "  # interconnect (" << remote_pairs.size() << " remote CX)\n";
            }
        } else {
            oss << line << "\n";
        }
    }

    // Save annotated text (with comments) before parsing, since Stim strips comments
    annotated_circuit_str_ = oss.str();

    // Parse modified circuit
    circuit_ = stim::Circuit(oss.str().c_str());
}

uint32_t SurfaceCodeSimulator::count_remote_cnots_in_cycle() {
    if (!config_.distributed) {
        return 0;
    }

    auto coords = circuit_.get_final_qubit_coords();
    double split_x = config_.code_distance;

    // Find the REPEAT block (main stabilizer cycles)
    for (const auto& op : circuit_.operations) {
        if (op.gate_type == stim::GateType::REPEAT) {
            // Get the body of one cycle
            const stim::Circuit& cycle_body = op.repeat_block_body(circuit_);

            uint32_t remote_count = 0;

            // Count crossing CNOTs in this cycle
            for (const auto& cycle_op : cycle_body.operations) {
                if (cycle_op.gate_type == stim::GateType::CX) {
                    // Check each CNOT pair
                    for (size_t i = 0; i + 1 < cycle_op.targets.size(); i += 2) {
                        uint32_t ctrl = cycle_op.targets[i].qubit_value();
                        uint32_t tgt = cycle_op.targets[i + 1].qubit_value();

                        if (coords.count(ctrl) && coords.count(tgt)) {
                            double ctrl_x = coords[ctrl][0];
                            double tgt_x = coords[tgt][0];

                            bool crosses = (ctrl_x < split_x && tgt_x >= split_x) ||
                                          (ctrl_x >= split_x && tgt_x < split_x);

                            if (crosses) {
                                remote_count++;
                            }
                        }
                    }
                }
            }

            return remote_count;
        }
    }

    return 0;
}

// Paper Equation 3: 2→1 Distillation output fidelity (DEJMPS bilateral protocol)
// Uses the standard bilateral CNOT distillation protocol
double SurfaceCodeSimulator::compute_2to1_fidelity(double F1, double F2) {
    double p_s = compute_2to1_success_prob(F1, F2);
    if (p_s <= 0) return F2;  // Return target fidelity if distillation fails

    // Standard DEJMPS formula: F_out = (F1×F2 + (1-F1)(1-F2)/9) / p_s
    double numerator = (F1 * F2) + ((1.0 - F1) * (1.0 - F2) / 9.0);
    return std::min(1.0, std::max(0.0, numerator / p_s));
}

// Paper Equation 4: 2→1 Distillation success probability
double SurfaceCodeSimulator::compute_2to1_success_prob(double F1, double F2) {
    // Standard DEJMPS success probability
    // p_s = F1×F2 + (1-F1)(1-F2)/9 + 2(F1(1-F2) + (1-F1)F2)/3 × (1/3)
    // Simplified: p_s ≈ (1 + 2F1 - 1)(1 + 2F2 - 1)/9 + rest...
    // Using standard form: p_s = (1 - 2(1-F1)/3)(1 - 2(1-F2)/3) + 4(1-F1)(1-F2)/9
    double term1 = (1.0 - 2.0 * (1.0 - F1) / 3.0) * (1.0 - 2.0 * (1.0 - F2) / 3.0);
    double term2 = 4.0 * (1.0 - F1) * (1.0 - F2) / 9.0;
    return std::max(0.01, term1 + term2);  // Minimum 1% success rate
}

// Paper Equation 5: 3→1 Distillation output fidelity
// Using generalized pumping formula for 3→1 protocol
double SurfaceCodeSimulator::compute_3to1_fidelity(double F1, double F2) {
    double p_s = compute_3to1_success_prob(F1, F2);
    if (p_s <= 0) return F2;

    // 3→1 distillation: uses two auxiliary pairs (F1) to improve target pair (F2)
    // The output fidelity follows a higher-order improvement
    // F_out ≈ (F1² × F2 + correction terms) / p_s
    // Using a standard approximation for high-fidelity regime:
    double F1_sq = F1 * F1;
    double one_minus_F1 = 1.0 - F1;

    // Numerator: weighted contribution from successful branches
    double numerator = F1_sq * F2 + F1_sq * (1.0 - F2) / 27.0 +
                       2.0 * F1 * one_minus_F1 * F2 / 9.0 +
                       one_minus_F1 * one_minus_F1 * F2 / 81.0;

    return std::min(1.0, std::max(0.0, numerator / p_s));
}

// Paper Equation 6: 3→1 Distillation success probability
double SurfaceCodeSimulator::compute_3to1_success_prob(double F1, double F2) {
    // Success probability for 3→1 protocol
    // Approximation based on independent success of two 2→1 stages
    double F1_sq = F1 * F1;
    double one_minus_F1 = 1.0 - F1;

    // p_s = sum of all matching measurement outcome probabilities
    double p_s = F1_sq * F2 + F1_sq * (1.0 - F2) / 9.0 +
                 2.0 * F1 * one_minus_F1 * (F2 / 3.0 + (1.0 - F2) / 9.0) +
                 one_minus_F1 * one_minus_F1 * (F2 / 9.0 + (1.0 - F2) / 81.0);

    return std::max(0.01, p_s);  // Minimum 1% success rate
}

// Compute full distillation result for given number of rounds
DistillationResult SurfaceCodeSimulator::compute_distillation(uint32_t rounds) const {
    DistillationResult result;
    result.output_fidelity = config_.raw_epr_fidelity;
    result.success_probability = 1.0;
    result.raw_pairs_consumed = 1;
    result.distillation_time = 0.0;

    if (!config_.distributed || config_.distillation_protocol == DistillationProtocol::NONE || rounds == 0) {
        return result;
    }

    double F_raw = config_.raw_epr_fidelity;
    double F_current = F_raw;
    double cumulative_success = 1.0;
    uint32_t total_pairs = 0;

    for (uint32_t k = 0; k < rounds; k++) {
        double F_out, p_s;
        uint32_t pairs_this_round;

        switch (config_.distillation_protocol) {
            case DistillationProtocol::PUMPING_2TO1:
                // Pumping: auxiliary pairs are always raw, target improves
                F_out = compute_2to1_fidelity(F_raw, F_current);
                p_s = compute_2to1_success_prob(F_raw, F_current);
                pairs_this_round = 2;  // Linear: 2 pairs per round
                break;

            case DistillationProtocol::PUMPING_3TO1:
                // Pumping: two auxiliary pairs (raw), one target
                F_out = compute_3to1_fidelity(F_raw, F_current);
                p_s = compute_3to1_success_prob(F_raw, F_current);
                pairs_this_round = 3;  // Linear: 3 pairs per round
                break;

            case DistillationProtocol::RECURRENCE_2TO1:
                // Recurrence: both inputs are distilled from previous round
                F_out = compute_2to1_fidelity(F_current, F_current);
                p_s = compute_2to1_success_prob(F_current, F_current);
                pairs_this_round = (k == 0) ? 2 : (1u << (k + 1));  // Exponential: 2^(k+1)
                break;

            case DistillationProtocol::RECURRENCE_3TO1:
                // Recurrence: all three inputs distilled
                F_out = compute_3to1_fidelity(F_current, F_current);
                p_s = compute_3to1_success_prob(F_current, F_current);
                pairs_this_round = (k == 0) ? 3 : static_cast<uint32_t>(std::pow(3, k + 1));  // 3^(k+1)
                break;

            default:
                return result;
        }

        F_current = std::max(0.0, std::min(1.0, F_out));
        cumulative_success *= p_s;
        total_pairs = pairs_this_round;  // For pumping it's cumulative, for recurrence it's the tree size
    }

    // For pumping, total pairs is linear: b*k (where b=2 or 3)
    // For recurrence, total pairs is exponential: b^k
    if (config_.distillation_protocol == DistillationProtocol::PUMPING_2TO1) {
        total_pairs = 2 * rounds;  // O(b*k) for pumping
    } else if (config_.distillation_protocol == DistillationProtocol::PUMPING_3TO1) {
        total_pairs = 3 * rounds;
    }
    // Recurrence already computed as b^k above

    result.output_fidelity = F_current;
    result.success_probability = cumulative_success;
    result.raw_pairs_consumed = std::max(1u, total_pairs);

    // Distillation time: time to generate required raw pairs
    // Time = raw_pairs / (entanglement_rate × success_probability)
    if (config_.entanglement_rate > 0 && cumulative_success > 0) {
        result.distillation_time = static_cast<double>(result.raw_pairs_consumed) /
                                   (config_.entanglement_rate * cumulative_success);
    }

    return result;
}

// Paper Equation 1: Remote CNOT error from EPR fidelity
double SurfaceCodeSimulator::compute_cnot_error_from_fidelity(double fidelity) const {
    // P_RemoteCX = 1 - F × (1 - p)²
    // where p is the local gate error rate
    double p = config_.physical_error;
    double p_remote = 1.0 - fidelity * (1.0 - p) * (1.0 - p);
    return std::max(0.0, p_remote);
}

// Compute time required for distillation
double SurfaceCodeSimulator::compute_distillation_time(uint32_t num_remote_cnots,
                                                        const DistillationResult& distill) const {
    if (config_.entanglement_rate <= 0) return 0.0;

    // Total EPR pairs needed = num_remote_cnots × raw_pairs_per_distilled
    // Account for success probability: expected attempts = 1 / p_s
    double expected_pairs = static_cast<double>(num_remote_cnots * distill.raw_pairs_consumed) /
                           distill.success_probability;

    // Time = pairs / rate
    return expected_pairs / config_.entanglement_rate;
}

// Check if distillation fits within QEC cycle timing constraint (paper Section V.A)
bool SurfaceCodeSimulator::check_timing_constraint(double distillation_time) const {
    // Constraint: T_dist ≤ t_meas (Equation 9)
    return distillation_time <= config_.measurement_delay;
}

void SurfaceCodeSimulator::inject_entanglement_idling_noise() {
    if (!config_.distributed || config_.entanglement_rate <= 0) {
        return;
    }

    // Count remote CNOTs per cycle
    uint32_t N_remote = count_remote_cnots_in_cycle();

    if (N_remote == 0) {
        return;  // No remote gates, no idling
    }

    // Compute distillation parameters
    DistillationResult distill = compute_distillation(config_.distillation_rounds);
    double t_dist = compute_distillation_time(N_remote, distill);

    // Paper Section III.B.2: tidle = max(tmeas, EPR waiting time)
    // For monolithic: tidle = tmeas
    // For distributed: tidle = max(tmeas, t_dist)
    double t_idle = std::max(config_.measurement_delay, t_dist);

    // Check timing constraint (paper Equation 9)
    bool timing_safe = check_timing_constraint(t_dist);

    // Paper Equation 2: Idling error using Pauli twirl approximation
    // pX = pY = (1 - e^(-tidle/T1)) / 4
    // pZ = (1 - e^(-tidle/T2)) / 2 - pX
    double p_X = 0.0, p_Y = 0.0, p_Z = 0.0;

    if (config_.T1_coherence_time > 0) {
        p_X = (1.0 - std::exp(-t_idle / config_.T1_coherence_time)) / 4.0;
        p_Y = p_X;
    }
    if (config_.T2_coherence_time > 0) {
        p_Z = (1.0 - std::exp(-t_idle / config_.T2_coherence_time)) / 2.0 - p_X;
        p_Z = std::max(0.0, p_Z);  // Ensure non-negative
    }

    // Total idle error for depolarizing approximation
    double p_idle_total = p_X + p_Y + p_Z;

    if (mpi_rank_ == 0) {
        std::cout << "Entanglement-limited idling (Eq. 2):" << std::endl;
        std::cout << "  Remote CNOTs per cycle: " << N_remote << std::endl;
        std::cout << "  Raw pairs per distilled: " << distill.raw_pairs_consumed << std::endl;
        std::cout << "  Distillation time: " << (t_dist * 1e9) << " ns" << std::endl;
        std::cout << "  Measurement delay: " << (config_.measurement_delay * 1e9) << " ns" << std::endl;
        std::cout << "  Idling time (max): " << (t_idle * 1e6) << " μs" << std::endl;
        std::cout << "  Timing constraint satisfied: " << (timing_safe ? "YES" : "NO - cycle extended!") << std::endl;
        std::cout << "  Idling errors (pX, pY, pZ): (" << p_X << ", " << p_Y << ", " << p_Z << ")" << std::endl;
        std::cout << "  Total idling error: " << p_idle_total << std::endl;
    }

    if (p_idle_total <= 0) {
        return;  // No idling error to inject
    }

    // Get all qubits
    size_t num_qubits = circuit_.count_qubits();
    std::vector<uint32_t> all_qubits;
    for (uint32_t i = 0; i < num_qubits; i++) {
        all_qubits.push_back(i);
    }

    // Modify circuit to add idling noise to REPEAT blocks
    stim::Circuit modified;

    for (const auto& op : circuit_.operations) {
        if (op.gate_type == stim::GateType::REPEAT) {
            uint64_t repeat_count = op.repeat_block_rep_count();
            stim::Circuit cycle_body = op.repeat_block_body(circuit_);

            // Prepend idling noise at start of cycle
            // (models waiting for entanglement before executing gates)
            stim::Circuit new_body;

            // Use PAULI_CHANNEL_1 for accurate T1/T2 noise (pX, pY, pZ)
            // If not available, fall back to DEPOLARIZE1
            std::vector<double> pauli_probs = {p_X, p_Y, p_Z};
            new_body.safe_append_u("PAULI_CHANNEL_1", all_qubits, pauli_probs);

            // Then add original cycle operations
            for (const auto& cycle_op : cycle_body.operations) {
                new_body.safe_append(cycle_op);
            }

            // Add modified REPEAT block
            modified.append_repeat_block(repeat_count, new_body, "");
        } else {
            // Copy non-REPEAT operations unchanged
            modified.safe_append(op);
        }
    }

    circuit_ = modified;
}

void SurfaceCodeSimulator::run_monte_carlo() {
    uint64_t remaining_shots = shots_this_rank_;
    size_t num_detectors = circuit_.count_detectors();
    size_t num_observables = circuit_.count_observables();

    // Random number generator
    std::mt19937_64 rng(std::random_device{}() + mpi_rank_);

    while (remaining_shots > 0) {
        uint64_t batch_shots = std::min(remaining_shots, BATCH_SIZE);

        // Sample detection events and observables from circuit
        auto sample_result = stim::sample_batch_detection_events<stim::MAX_BITWORD_WIDTH>(
            circuit_,
            batch_shots,
            rng
        );

        auto& detection_events = sample_result.first;
        auto& observable_flips = sample_result.second;

        // Decode the batch
        uint64_t batch_errors = decode_batch(
            mwpm_decoder_,
            detection_events,
            observable_flips,
            batch_shots,
            num_detectors
        );

        local_errors_ += batch_errors;
        remaining_shots -= batch_shots;

        // Progress reporting (only rank 0, occasionally)
        if (mpi_rank_ == 0 && (shots_this_rank_ - remaining_shots) % (10 * BATCH_SIZE) == 0) {
            double progress = 100.0 * (shots_this_rank_ - remaining_shots) / shots_this_rank_;
            std::cout << "Rank 0 progress: " << progress << "%" << std::endl;
        }
    }
}

void SurfaceCodeSimulator::run_bucket() {
    // Analyze DEM to extract error mechanisms
    auto dem = stim::ErrorAnalyzer::circuit_to_detector_error_model(
        circuit_,
        false,  // decompose_errors
        true,   // fold_loops
        false,  // allow_gauge_detectors
        0,      // approximate_disjoint_errors_threshold
        false,  // ignore_decomposition_failures
        false   // block_decomposition_from_introducing_remnant_edges
    );

    DEMAnalyzer analyzer(dem);
    auto errors = analyzer.get_errors();

    if (mpi_rank_ == 0) {
        std::cout << "Bucket mode: " << errors.size() << " error mechanisms" << std::endl;
        std::cout << "Total error rate: " << analyzer.get_total_error_rate() << std::endl;
    }

    // Compute bucket probabilities
    uint32_t max_bucket = config_.max_bucket;
    if (max_bucket == 0) {
        // Auto-detect: use 3*lambda + 5 as cutoff for good coverage
        max_bucket = static_cast<uint32_t>(
            3.0 * analyzer.get_total_error_rate() + 5
        );
    }

    if (mpi_rank_ == 0) {
        std::cout << "Max bucket: " << max_bucket << std::endl;
    }

    auto analysis = analyzer.compute_bucket_probabilities(
        max_bucket, shots_this_rank_, config_.num_sampled_buckets, config_.max_bias_bound
    );
    auto& buckets = analysis.buckets;

    if (mpi_rank_ == 0) {
        std::cout << "Sampled buckets: " << buckets.size() << std::endl;
        std::cout << "Sampled probability mass: " << std::fixed << std::setprecision(6)
                  << analysis.sampled_probability_mass << std::endl;
        std::cout << "Bias bound: " << std::scientific << std::setprecision(2)
                  << analysis.bias_bound << std::endl;
    }

    // Create bucket sampler
    BucketSampler sampler(
        circuit_.count_detectors(),
        circuit_.count_observables(),
        mwpm_decoder_,
        errors,
        std::random_device{}() + mpi_rank_
    );

    // Sample each bucket
    double weighted_error_count = 0.0;
    uint64_t total_samples = 0;

    for (auto& bucket : buckets) {
        if (bucket.target_samples == 0) continue;

        if (mpi_rank_ == 0 && bucket.error_count <= 5) {
            std::cout << "Bucket " << bucket.error_count << ": "
                     << "P=" << bucket.probability << ", "
                     << "samples=" << bucket.target_samples << std::endl;
        }

        for (uint64_t shot = 0; shot < bucket.target_samples; shot++) {
            uint64_t is_error = sampler.sample_shot_with_k_errors(bucket.error_count);
            bucket.logical_errors += is_error;
        }

        bucket.actual_samples = bucket.target_samples;
        total_samples += bucket.actual_samples;

        // Weighted contribution to total error count
        double bucket_ler = static_cast<double>(bucket.logical_errors) / bucket.actual_samples;
        weighted_error_count += bucket.probability * bucket_ler * shots_this_rank_;
    }

    // Compute statistical error bars
    // SE = sqrt[(1/N) × Σ P(k) × LER_k × (1 - LER_k)]
    double variance_sum = 0.0;
    for (const auto& bucket : buckets) {
        if (bucket.actual_samples == 0) continue;
        double bucket_ler = static_cast<double>(bucket.logical_errors) / bucket.actual_samples;
        variance_sum += bucket.probability * bucket_ler * (1.0 - bucket_ler);
    }
    double standard_error = std::sqrt(variance_sum / shots_this_rank_);

    local_errors_ = static_cast<uint64_t>(weighted_error_count);

    // Store bucket statistics
    num_sampled_buckets_ = buckets.size();
    sampled_probability_mass_ = analysis.sampled_probability_mass;
    bias_bound_ = analysis.bias_bound;
    statistical_error_ = standard_error;
    has_bucket_stats_ = true;

    if (mpi_rank_ == 0) {
        std::cout << "Total samples taken: " << total_samples << std::endl;
        double estimated_ler = weighted_error_count / shots_this_rank_;
        std::cout << "Estimated LER: " << std::fixed << std::setprecision(6)
                  << estimated_ler << " ± " << std::setprecision(6)
                  << (2.0 * standard_error) << " (95% CI)" << std::endl;
    }
}

void SurfaceCodeSimulator::run() {
    start_time_ = std::chrono::high_resolution_clock::now();

    if (config_.mode == SimulationMode::MONTE_CARLO) {
        run_monte_carlo();
    } else if (config_.mode == SimulationMode::BUCKET) {
        run_bucket();
    }

    end_time_ = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        end_time_ - start_time_
    );
    local_runtime_ = duration.count() / 1000.0; // Convert to seconds
}

RankStats SurfaceCodeSimulator::get_rank_stats() const {
    RankStats stats;
    stats.rank = mpi_rank_;
    stats.shots = shots_this_rank_;
    stats.runtime_seconds = local_runtime_;
    stats.shots_per_second = local_runtime_ > 0 ? shots_this_rank_ / local_runtime_ : 0.0;
    return stats;
}

} // namespace bucket_sim
