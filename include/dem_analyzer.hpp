#pragma once

#include <stim.h>
#include <vector>
#include <cstdint>

namespace bucket_sim {

struct ErrorMechanism {
    uint64_t index;
    double probability;
    std::vector<uint64_t> detectors;
    std::vector<uint64_t> observables;
};

struct BucketInfo {
    uint32_t error_count;
    double probability;
    uint64_t target_samples;
    uint64_t actual_samples = 0;
    uint64_t logical_errors = 0;
};

struct BucketAnalysis {
    std::vector<BucketInfo> buckets;
    double sampled_probability_mass;  // S_sampled
    double bias_bound;  // 1 - S_sampled
};

class DEMAnalyzer {
public:
    explicit DEMAnalyzer(const stim::DetectorErrorModel& dem);

    const std::vector<ErrorMechanism>& get_errors() const { return errors_; }
    double get_total_error_rate() const { return total_error_rate_; }

    BucketAnalysis compute_bucket_probabilities(
        uint32_t max_bucket,
        uint64_t total_shots,
        uint32_t num_sampled_buckets = 0,  // 0 = all buckets
        double max_bias_bound = 0.0  // 0 = disabled, >0 = auto-select to keep bias <= this
    ) const;

private:
    std::vector<ErrorMechanism> errors_;
    double total_error_rate_;
};

} // namespace bucket_sim
