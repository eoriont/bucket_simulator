#pragma once

#include <string>
#include <cstdint>
#include <vector>

namespace bucket_sim {

enum class SimulationMode {
    MONTE_CARLO,
    BUCKET
};

// Lattice surgery merge type
enum class MergeType {
    NONE,                // No lattice surgery (standard memory experiment)
    XX_MERGE,            // X-basis merge (measures XX between patches)
    ZZ_MERGE,            // Z-basis merge (measures ZZ between patches)
    XX_MERGE_DISTRIBUTED // Distributed XX merge: no merge data qubits, remote CNOTs only, flipped Patch B
};

enum class DistillationProtocol {
    NONE,
    PUMPING_2TO1,      // 2→1 Pumping: linear resource cost O(b*k)
    PUMPING_3TO1,      // 3→1 Pumping: linear resource cost O(b*k)
    RECURRENCE_2TO1,   // 2→1 Recurrence: exponential resource cost O(2^k)
    RECURRENCE_3TO1    // 3→1 Recurrence: exponential resource cost O(3^k)
};

struct Config {
    uint32_t code_distance;
    uint32_t rounds;
    double physical_error;
    double measurement_error;
    double reset_error;
    uint64_t total_shots;
    std::string code_type;

    // Bucket mode parameters
    SimulationMode mode;
    uint32_t max_bucket;
    uint64_t min_shots_per_bucket;
    uint32_t num_sampled_buckets;  // 0 = all buckets
    double max_bias_bound;  // 0.0 = disabled, >0 = auto-select buckets to keep bias <= this

    // Distributed QEC parameters
    bool distributed;  // false = monolithic (default), true = distributed
    double interconnect_error;  // Additional error rate for boundary CNOTs

    // Entanglement-limited idling parameters
    double entanglement_rate;  // EPR pairs generated per second (Hz)
    double T1_coherence_time;  // T1 relaxation time in seconds
    double T2_coherence_time;  // T2 dephasing time in seconds

    // Timing parameters (from paper Section III.B)
    double measurement_delay;  // Measurement delay in seconds (default: 660ns)

    // Entanglement distillation parameters
    double raw_epr_fidelity;          // Raw EPR pair fidelity (0 to 1)
    DistillationProtocol distillation_protocol;
    uint32_t distillation_rounds;     // Number of distillation rounds (k)

    // Lattice surgery parameters
    MergeType merge_type;             // Type of merge operation
    uint32_t merge_rounds;            // Number of stabilizer rounds during merge
    std::vector<uint32_t> superstab_ys; // y-coords in merge column to superstabilize (XX_MERGE_DISTRIBUTED only)
    bool split_after_merge;           // Whether to split patches after merge

    // Default constructor
    Config()
        : code_distance(5),
          rounds(5),
          physical_error(0.001),
          measurement_error(0.001),
          reset_error(0.001),
          total_shots(1000000),
          code_type("rotated_memory_x"),
          mode(SimulationMode::MONTE_CARLO),
          max_bucket(0),
          min_shots_per_bucket(1000),
          num_sampled_buckets(0),
          max_bias_bound(0.0),
          distributed(false),
          interconnect_error(0.0),
          entanglement_rate(100e6),    // Default: 100 MHz
          T1_coherence_time(125e-6),   // Default: 125 μs (paper Section III.B)
          T2_coherence_time(200e-6),   // Default: 200 μs (paper Section III.B)
          measurement_delay(660e-9),   // Default: 660 ns (paper Section III.B)
          raw_epr_fidelity(0.99),      // Default: 99% raw fidelity (paper Section VI)
          distillation_protocol(DistillationProtocol::NONE),
          distillation_rounds(1),      // Default: 1 round
          merge_type(MergeType::NONE),
          merge_rounds(0),             // 0 = use code_distance rounds
          split_after_merge(false) {}
};

// Parse configuration from file
Config parse_config(const std::string& filename);

// Helper function to parse magnitude suffixes (K, M, B, G)
uint64_t parse_magnitude(const std::string& input);

} // namespace bucket_sim
