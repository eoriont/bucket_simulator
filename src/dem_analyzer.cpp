#include "dem_analyzer.hpp"
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace bucket_sim {

DEMAnalyzer::DEMAnalyzer(const stim::DetectorErrorModel& dem) {
    total_error_rate_ = 0.0;

    // Flatten to expand repeat blocks
    auto flattened = dem.flattened();

    // Iterate error instructions
    flattened.iter_flatten_error_instructions([&](const stim::DemInstruction& inst) {
        if (inst.arg_data.size() == 0) return;

        ErrorMechanism mech;
        mech.index = errors_.size();
        mech.probability = inst.arg_data[0];

        // Extract detectors and observables
        for (const auto& target : inst.target_data) {
            if (target.is_relative_detector_id()) {
                mech.detectors.push_back(target.raw_id());
            } else if (target.is_observable_id()) {
                mech.observables.push_back(target.raw_id());
            }
        }

        errors_.push_back(mech);
        total_error_rate_ += mech.probability;
    });
}

BucketAnalysis DEMAnalyzer::compute_bucket_probabilities(
    uint32_t max_bucket,
    uint64_t total_shots,
    uint32_t num_sampled_buckets,
    double max_bias_bound
) const {
    std::vector<BucketInfo> all_buckets;

    // Poisson approximation: P(k) = λ^k * exp(-λ) / k!
    // Use log-space to avoid numerical overflow for large k
    double lambda = total_error_rate_;

    // Compute all bucket probabilities
    for (uint32_t k = 0; k <= max_bucket; k++) {
        BucketInfo bucket;
        bucket.error_count = k;

        // Compute log P(k) = k*log(λ) - λ - log(k!)
        // log(k!) = lgamma(k+1)
        double log_prob = k * std::log(lambda) - lambda - std::lgamma(k + 1);
        bucket.probability = std::exp(log_prob);

        // Skip buckets with negligible probability
        if (bucket.probability < 1e-15) {
            continue;
        }

        all_buckets.push_back(bucket);
    }

    // Sort buckets by probability (descending)
    std::sort(all_buckets.begin(), all_buckets.end(),
              [](const BucketInfo& a, const BucketInfo& b) {
                  return a.probability > b.probability;
              });

    // Select buckets based on either num_sampled_buckets or max_bias_bound
    size_t num_to_sample = all_buckets.size();

    if (max_bias_bound > 0.0) {
        // Auto-select buckets to keep bias_bound <= max_bias_bound
        double cumulative_prob = 0.0;
        for (size_t i = 0; i < all_buckets.size(); i++) {
            cumulative_prob += all_buckets[i].probability;
            double bias_bound = 1.0 - cumulative_prob;
            if (bias_bound <= max_bias_bound) {
                num_to_sample = i + 1;
                break;
            }
        }
    } else if (num_sampled_buckets > 0 && num_sampled_buckets < all_buckets.size()) {
        // Use explicit num_sampled_buckets
        num_to_sample = num_sampled_buckets;
    }

    std::vector<BucketInfo> sampled_buckets;
    double sampled_probability_mass = 0.0;

    for (size_t i = 0; i < num_to_sample; i++) {
        auto& bucket = all_buckets[i];
        sampled_probability_mass += bucket.probability;

        // Allocate samples proportional to probability
        bucket.target_samples = static_cast<uint64_t>(
            bucket.probability * total_shots + 0.5
        );

        // Ensure minimum samples if probability is non-negligible
        if (bucket.probability > 1e-10 && bucket.target_samples == 0) {
            bucket.target_samples = 100;  // Minimum for statistics
        }

        sampled_buckets.push_back(bucket);
    }

    BucketAnalysis analysis;
    analysis.buckets = sampled_buckets;
    analysis.sampled_probability_mass = sampled_probability_mass;
    analysis.bias_bound = 1.0 - sampled_probability_mass;

    return analysis;
}

} // namespace bucket_sim
