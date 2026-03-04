#include "config.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <iostream>

namespace bucket_sim {

uint64_t parse_magnitude(const std::string& input) {
    std::unordered_map<char, uint64_t> magnitude_map = {
        {'K', 1000ULL},
        {'M', 1000000ULL},
        {'B', 1000000000ULL},
        {'G', 1000000000ULL}
    };

    if (input.empty()) {
        throw std::invalid_argument("Empty input string");
    }

    char last_char = input.back();
    auto it = magnitude_map.find(last_char);

    if (it != magnitude_map.end()) {
        // Last character is a magnitude suffix
        std::string number_part = input.substr(0, input.length() - 1);
        uint64_t value = std::stoull(number_part);
        return value * it->second;
    } else {
        // No magnitude suffix, parse as regular number
        return std::stoull(input);
    }
}

Config parse_config(const std::string& filename) {
    Config config;
    std::ifstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("Failed to open config file: " + filename);
    }

    std::string line;
    int line_number = 0;

    while (std::getline(file, line)) {
        line_number++;

        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') {
            continue;
        }

        std::istringstream ss(line);
        std::string key, value;
        ss >> key;

        // Get the rest of the line as value
        std::getline(ss, value);

        // Trim leading whitespace from value
        value.erase(0, value.find_first_not_of(" \t"));

        // Parse key-value pairs
        try {
            if (key == "code_distance") {
                config.code_distance = std::stoul(value);
            } else if (key == "rounds") {
                config.rounds = std::stoul(value);
            } else if (key == "physical_error") {
                config.physical_error = std::stod(value);
            } else if (key == "measurement_error") {
                config.measurement_error = std::stod(value);
            } else if (key == "reset_error") {
                config.reset_error = std::stod(value);
            } else if (key == "total_shots") {
                config.total_shots = parse_magnitude(value);
            } else if (key == "code_type") {
                config.code_type = value;
            } else if (key == "mode") {
                if (value == "monte_carlo" || value == "mc") {
                    config.mode = SimulationMode::MONTE_CARLO;
                } else if (value == "bucket") {
                    config.mode = SimulationMode::BUCKET;
                } else {
                    throw std::invalid_argument("Unknown simulation mode: " + value);
                }
            } else if (key == "max_bucket") {
                config.max_bucket = std::stoul(value);
            } else if (key == "min_shots_per_bucket") {
                config.min_shots_per_bucket = parse_magnitude(value);
            } else if (key == "num_sampled_buckets") {
                config.num_sampled_buckets = std::stoul(value);
            } else if (key == "max_bias_bound") {
                config.max_bias_bound = std::stod(value);
            } else if (key == "distributed") {
                config.distributed = (value == "true" || value == "1");
            } else if (key == "interconnect_error") {
                config.interconnect_error = std::stod(value);
            } else if (key == "entanglement_rate") {
                config.entanglement_rate = std::stod(value);
            } else if (key == "T1") {
                config.T1_coherence_time = std::stod(value);
            } else if (key == "T2") {
                config.T2_coherence_time = std::stod(value);
            } else if (key == "measurement_delay") {
                config.measurement_delay = std::stod(value);
            } else if (key == "raw_epr_fidelity") {
                config.raw_epr_fidelity = std::stod(value);
            } else if (key == "distillation_protocol") {
                if (value == "none") {
                    config.distillation_protocol = DistillationProtocol::NONE;
                } else if (value == "pumping_2to1" || value == "2to1_pumping") {
                    config.distillation_protocol = DistillationProtocol::PUMPING_2TO1;
                } else if (value == "pumping_3to1" || value == "3to1_pumping") {
                    config.distillation_protocol = DistillationProtocol::PUMPING_3TO1;
                } else if (value == "recurrence_2to1" || value == "2to1_recurrence") {
                    config.distillation_protocol = DistillationProtocol::RECURRENCE_2TO1;
                } else if (value == "recurrence_3to1" || value == "3to1_recurrence") {
                    config.distillation_protocol = DistillationProtocol::RECURRENCE_3TO1;
                } else {
                    throw std::invalid_argument("Unknown distillation protocol: " + value +
                        ". Valid options: none, pumping_2to1, pumping_3to1, recurrence_2to1, recurrence_3to1");
                }
            } else if (key == "distillation_rounds") {
                config.distillation_rounds = std::stoul(value);
            } else if (key == "merge_type") {
                if (value == "none") {
                    config.merge_type = MergeType::NONE;
                } else if (value == "xx" || value == "XX" || value == "xx_merge") {
                    config.merge_type = MergeType::XX_MERGE;
                } else if (value == "zz" || value == "ZZ" || value == "zz_merge") {
                    config.merge_type = MergeType::ZZ_MERGE;
                } else if (value == "distributed_xx" || value == "xx_distributed" || value == "DISTRIBUTED_XX") {
                    config.merge_type = MergeType::XX_MERGE_DISTRIBUTED;
                } else {
                    throw std::invalid_argument("Unknown merge type: " + value +
                        ". Valid options: none, xx, zz, distributed_xx");
                }
            } else if (key == "merge_rounds") {
                config.merge_rounds = std::stoul(value);
            } else if (key == "split_after_merge") {
                config.split_after_merge = (value == "true" || value == "1");
            } else if (key == "superstab_ys") {
                config.superstab_ys.clear();
                std::istringstream vs(value);
                uint32_t y;
                while (vs >> y) config.superstab_ys.push_back(y);
            } else {
                std::cerr << "Warning: Unknown config key '" << key
                         << "' at line " << line_number << std::endl;
            }
        } catch (const std::exception& e) {
            throw std::runtime_error("Error parsing config at line "
                                   + std::to_string(line_number)
                                   + ": " + e.what());
        }
    }

    // Validate configuration
    if (config.code_distance == 0) {
        throw std::invalid_argument("code_distance must be greater than 0");
    }
    if (config.rounds == 0) {
        throw std::invalid_argument("rounds must be greater than 0");
    }
    if (config.physical_error < 0.0 || config.physical_error > 1.0) {
        throw std::invalid_argument("physical_error must be in [0, 1]");
    }
    if (config.measurement_error < 0.0 || config.measurement_error > 1.0) {
        throw std::invalid_argument("measurement_error must be in [0, 1]");
    }
    if (config.reset_error < 0.0 || config.reset_error > 1.0) {
        throw std::invalid_argument("reset_error must be in [0, 1]");
    }
    if (config.total_shots == 0) {
        throw std::invalid_argument("total_shots must be greater than 0");
    }

    return config;
}

} // namespace bucket_sim
