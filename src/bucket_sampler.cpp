#include "bucket_sampler.hpp"
#include <pymatching/sparse_blossom/driver/mwpm_decoding.h>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <random>

namespace bucket_sim {

BucketSampler::BucketSampler(
    size_t num_detectors,
    size_t num_observables,
    pm::Mwpm& decoder,
    const std::vector<ErrorMechanism>& errors,
    int seed
) : num_detectors_(num_detectors),
    num_observables_(num_observables),
    decoder_(decoder),
    errors_(errors),
    rng_(seed) {}

std::vector<uint64_t> BucketSampler::select_k_errors(uint32_t k) {
    if (k > errors_.size()) {
        throw std::invalid_argument("k exceeds number of available errors");
    }

    if (k == 0) {
        return std::vector<uint64_t>();
    }

    // Weighted sampling without replacement using Exponential/Gumbel method
    // For each error i, compute key_i = U(0,1)^(1/p_i), then take k largest
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    std::vector<std::pair<double, uint64_t>> keys;
    keys.reserve(errors_.size());

    for (uint64_t i = 0; i < errors_.size(); i++) {
        double u = uniform(rng_);
        // Avoid division by zero or log of zero
        if (errors_[i].probability > 0 && u > 0) {
            double key = std::pow(u, 1.0 / errors_[i].probability);
            keys.push_back({key, i});
        }
    }

    // Take k errors with largest keys
    std::partial_sort(keys.begin(), keys.begin() + k, keys.end(),
                     [](const auto& a, const auto& b) { return a.first > b.first; });

    std::vector<uint64_t> indices;
    indices.reserve(k);
    for (uint32_t i = 0; i < k && i < keys.size(); i++) {
        indices.push_back(keys[i].second);
    }

    return indices;
}

void BucketSampler::apply_errors(
    const std::vector<uint64_t>& error_indices,
    std::vector<bool>& detectors,
    std::vector<bool>& observables
) {
    // Initialize to false
    std::fill(detectors.begin(), detectors.end(), false);
    std::fill(observables.begin(), observables.end(), false);

    // XOR each error's effects
    for (uint64_t idx : error_indices) {
        const auto& error = errors_[idx];

        for (uint64_t det : error.detectors) {
            detectors[det] = !detectors[det];
        }

        for (uint64_t obs : error.observables) {
            observables[obs] = !observables[obs];
        }
    }
}

uint64_t BucketSampler::sample_shot_with_k_errors(uint32_t k) {
    // Select k random errors
    auto error_indices = select_k_errors(k);

    // Apply to get syndrome
    std::vector<bool> detectors(num_detectors_);
    std::vector<bool> observables(num_observables_);
    apply_errors(error_indices, detectors, observables);

    // Convert to sparse format (list of triggered detectors)
    std::vector<uint64_t> hits;
    for (size_t i = 0; i < detectors.size(); i++) {
        if (detectors[i]) {
            hits.push_back(i);
        }
    }

    // Decode
    auto result = pm::decode_detection_events_for_up_to_64_observables(
        decoder_, hits, false
    );

    // Check if decoding matches actual observable
    uint64_t actual_obs = observables[0] ? 1ULL : 0ULL;

    return (result.obs_mask != actual_obs) ? 1 : 0;
}

} // namespace bucket_sim
