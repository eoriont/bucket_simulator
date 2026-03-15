#pragma once

#include <stim.h>
#include <pymatching/sparse_blossom/driver/mwpm_decoding.h>
#include <pymatching/sparse_blossom/matcher/mwpm.h>
#include <vector>
#include <cstdint>

namespace bucket_sim {

// Decode a batch of detection events using MWPM
// Returns the number of logical errors (mismatches between decoded and actual observables)
uint64_t decode_batch(
    pm::Mwpm& mwpm,
    const stim::simd_bit_table<stim::MAX_BITWORD_WIDTH>& detection_events,
    const stim::simd_bit_table<stim::MAX_BITWORD_WIDTH>& observable_flips,
    size_t num_shots,
    size_t num_detectors,
    size_t num_observables
);

} // namespace bucket_sim
