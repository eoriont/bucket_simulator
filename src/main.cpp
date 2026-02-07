#include "config.hpp"
#include "simulator.hpp"
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <sstream>
#include <sys/stat.h>
#include <vector>
#include <optional>

struct BucketStats {
    uint32_t num_buckets;
    double sampled_probability_mass;
    double bias_bound;
    double statistical_error;
};

void print_usage(const char* program_name) {
    std::cout << "Usage: " << program_name << " -config <config_file> [-output <dir>] [-dump-circuit]" << std::endl;
    std::cout << "  -config <file>     Path to configuration file (required)" << std::endl;
    std::cout << "  -output <dir>      Output directory (default: <project_root>/output)" << std::endl;
    std::cout << "  -dump-circuit      Output the Stim circuit and exit (for visualization)" << std::endl;
}

std::string get_timestamp() {
    auto now = std::time(nullptr);
    auto tm = *std::localtime(&now);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y%m%d_%H%M%S");
    return oss.str();
}

void ensure_output_directory(const std::string& output_dir) {
    struct stat info;
    if (stat(output_dir.c_str(), &info) != 0) {
        // Directory doesn't exist, create it
        mkdir(output_dir.c_str(), 0755);
    }
}

void write_output_file(
    const std::string& output_dir,
    const bucket_sim::Config& config,
    uint64_t total_errors,
    double max_runtime,
    const std::vector<bucket_sim::RankStats>& all_rank_stats,
    const std::optional<BucketStats>& bucket_stats = std::nullopt
) {
    ensure_output_directory(output_dir);

    std::string filename = output_dir + "/results_" + get_timestamp() + ".txt";
    std::ofstream outfile(filename);

    if (!outfile.is_open()) {
        std::cerr << "Warning: Could not open output file " << filename << std::endl;
        return;
    }

    // Configuration section
    outfile << "Configuration:" << std::endl;
    outfile << "  Code Distance: " << config.code_distance << std::endl;
    outfile << "  Rounds: " << config.rounds << std::endl;
    outfile << "  Physical Error Rate: " << config.physical_error << std::endl;
    outfile << "  Measurement Error Rate: " << config.measurement_error << std::endl;
    outfile << "  Reset Error Rate: " << config.reset_error << std::endl;
    outfile << "  Total Shots: " << config.total_shots << std::endl;
    outfile << "  Code Type: " << config.code_type << std::endl;
    outfile << "  Simulation Mode: " << (config.mode == bucket_sim::SimulationMode::MONTE_CARLO ? "Monte Carlo" : "Bucket") << std::endl;
    if (config.distributed) {
        outfile << "  Distributed QEC: Yes" << std::endl;

        // Distillation info
        if (config.distillation_protocol != bucket_sim::DistillationProtocol::NONE) {
            outfile << "  Distillation Protocol: ";
            switch (config.distillation_protocol) {
                case bucket_sim::DistillationProtocol::PUMPING_2TO1:
                    outfile << "2→1 Pumping"; break;
                case bucket_sim::DistillationProtocol::PUMPING_3TO1:
                    outfile << "3→1 Pumping"; break;
                case bucket_sim::DistillationProtocol::RECURRENCE_2TO1:
                    outfile << "2→1 Recurrence"; break;
                case bucket_sim::DistillationProtocol::RECURRENCE_3TO1:
                    outfile << "3→1 Recurrence"; break;
                default:
                    outfile << "None"; break;
            }
            outfile << std::endl;
            outfile << "  Distillation Rounds: " << config.distillation_rounds << std::endl;
            outfile << "  Raw EPR Fidelity: " << config.raw_epr_fidelity << std::endl;
        }

        outfile << "  Interconnect Error Rate: " << config.interconnect_error << std::endl;

        if (config.entanglement_rate > 0) {
            outfile << "  Entanglement Rate: " << (config.entanglement_rate / 1e6) << " MHz" << std::endl;
            outfile << "  T1 Coherence Time: " << (config.T1_coherence_time * 1e6) << " μs" << std::endl;
            outfile << "  T2 Coherence Time: " << (config.T2_coherence_time * 1e6) << " μs" << std::endl;
            outfile << "  Measurement Delay: " << (config.measurement_delay * 1e9) << " ns" << std::endl;
        }
    }
    outfile << std::endl;

    // Results section
    double logical_error_rate = static_cast<double>(total_errors) / static_cast<double>(config.total_shots);
    outfile << "Results:" << std::endl;
    outfile << "  Logical Error Rate: " << std::fixed << std::setprecision(6) << logical_error_rate;
    outfile << " (" << total_errors << "/" << config.total_shots << ")" << std::endl;

    // Bucket mode statistics
    if (bucket_stats.has_value()) {
        outfile << "  95% Confidence Interval: ±" << std::fixed << std::setprecision(6)
                << (2.0 * bucket_stats->statistical_error) << std::endl;
        outfile << std::endl;
        outfile << "Bucket Sampling Statistics:" << std::endl;
        outfile << "  Sampled Buckets: " << bucket_stats->num_buckets << std::endl;
        outfile << "  Sampled Probability Mass: " << std::fixed << std::setprecision(6)
                << bucket_stats->sampled_probability_mass << std::endl;
        outfile << "  Truncation Bias Bound: " << std::scientific << std::setprecision(2)
                << bucket_stats->bias_bound << std::endl;
    }
    outfile << std::endl;

    // Runtime statistics
    int hours = static_cast<int>(max_runtime) / 3600;
    int minutes = (static_cast<int>(max_runtime) % 3600) / 60;
    int seconds = static_cast<int>(max_runtime) % 60;

    double total_shots_per_second = config.total_shots / max_runtime;

    outfile << "Runtime Statistics:" << std::endl;
    outfile << "  Total Runtime: " << hours << " hours, " << minutes << " minutes, " << seconds << " seconds" << std::endl;
    outfile << "  Shots per Second: " << std::fixed << std::setprecision(0) << total_shots_per_second << std::endl;
    outfile << std::endl;

    // Per-rank statistics
    outfile << "Per-Rank Statistics:" << std::endl;
    for (const auto& stats : all_rank_stats) {
        outfile << "  Rank " << stats.rank << ": "
                << stats.shots << " shots, "
                << std::fixed << std::setprecision(1) << stats.runtime_seconds << "s, "
                << std::fixed << std::setprecision(0) << stats.shots_per_second << " shots/s"
                << std::endl;
    }

    outfile.close();
    std::cout << "Results written to: " << filename << std::endl;
}

int main(int argc, char** argv) {
    // Initialize MPI
    MPI_Init(&argc, &argv);

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Parse command-line arguments
    std::string config_file;
    std::string output_dir = std::string(PROJECT_SOURCE_DIR) + "/output";  // Default to project root
    bool dump_circuit = false;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-config" && i + 1 < argc) {
            config_file = argv[++i];
        } else if (arg == "-output" && i + 1 < argc) {
            output_dir = argv[++i];
        } else if (arg == "-dump-circuit") {
            dump_circuit = true;
        } else if (arg == "-h" || arg == "--help") {
            if (world_rank == 0) {
                print_usage(argv[0]);
            }
            MPI_Finalize();
            return 0;
        }
    }

    if (config_file.empty()) {
        if (world_rank == 0) {
            std::cerr << "Error: -config argument is required" << std::endl;
            print_usage(argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    try {
        // Parse configuration
        bucket_sim::Config config = bucket_sim::parse_config(config_file);

        // Create simulator (this generates the circuit)
        bucket_sim::SurfaceCodeSimulator simulator(config, world_rank, world_size);

        // If dump-circuit flag is set, output circuit and exit
        if (dump_circuit) {
            if (world_rank == 0) {
                std::cout << simulator.get_circuit().str() << std::endl;
            }
            MPI_Finalize();
            return 0;
        }

        // Run simulation
        simulator.run();

        // Gather results
        uint64_t local_errors = simulator.get_local_errors();
        double local_runtime = simulator.get_local_runtime();

        // Reduce errors across all ranks
        uint64_t total_errors = 0;
        MPI_Reduce(&local_errors, &total_errors, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);

        // Get max runtime across all ranks
        double max_runtime = 0.0;
        MPI_Reduce(&local_runtime, &max_runtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        // Gather per-rank statistics
        bucket_sim::RankStats local_stats = simulator.get_rank_stats();
        std::vector<bucket_sim::RankStats> all_rank_stats;

        if (world_rank == 0) {
            all_rank_stats.resize(world_size);
        }

        MPI_Gather(&local_stats, sizeof(bucket_sim::RankStats), MPI_BYTE,
                   all_rank_stats.data(), sizeof(bucket_sim::RankStats), MPI_BYTE,
                   0, MPI_COMM_WORLD);

        // Rank 0 writes output
        if (world_rank == 0) {
            std::cout << std::endl;
            std::cout << "===============================================" << std::endl;
            std::cout << "Simulation Complete" << std::endl;
            std::cout << "===============================================" << std::endl;

            double logical_error_rate = static_cast<double>(total_errors) / static_cast<double>(config.total_shots);
            std::cout << "Logical Error Rate: " << std::fixed << std::setprecision(6) << logical_error_rate << std::endl;
            std::cout << "(" << total_errors << "/" << config.total_shots << ")" << std::endl;

            // Get bucket statistics if available
            std::optional<BucketStats> bucket_stats;
            if (simulator.has_bucket_stats()) {
                bucket_stats = BucketStats{
                    simulator.get_num_sampled_buckets(),
                    simulator.get_sampled_probability_mass(),
                    simulator.get_bias_bound(),
                    simulator.get_statistical_error()
                };
            }

            write_output_file(output_dir, config, total_errors, max_runtime, all_rank_stats, bucket_stats);
        }

    } catch (const std::exception& e) {
        if (world_rank == 0) {
            std::cerr << "Error: " << e.what() << std::endl;
        }
        MPI_Finalize();
        return 1;
    }

    MPI_Finalize();
    return 0;
}
