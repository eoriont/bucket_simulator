#pragma once

#include "dem_analyzer.hpp"
#include <stim.h>
#include <pymatching/sparse_blossom/matcher/mwpm.h>
#include <random>

namespace bucket_sim {

class BucketSampler {
public:
    BucketSampler(
        size_t num_detectors,
        size_t num_observables,
        pm::Mwpm& decoder,
        const std::vector<ErrorMechanism>& errors,
        int seed
    );

    // Sample a single shot with exactly k errors
    // Returns 1 if logical error, 0 otherwise
    uint64_t sample_shot_with_k_errors(uint32_t k);

private:
    // Generate k random error indices
    std::vector<uint64_t> select_k_errors(uint32_t k);

    // Apply errors to get detector/observable flips
    void apply_errors(
        const std::vector<uint64_t>& error_indices,
        std::vector<bool>& detectors,
        std::vector<bool>& observables
    );

    size_t num_detectors_;
    size_t num_observables_;
    pm::Mwpm& decoder_;
    const std::vector<ErrorMechanism>& errors_;
    std::mt19937_64 rng_;
};

} // namespace bucket_sim
