#include "decoder.hpp"

namespace bucket_sim {

uint64_t decode_batch(
    pm::Mwpm& mwpm,
    const stim::simd_bit_table<stim::MAX_BITWORD_WIDTH>& detection_events,
    const stim::simd_bit_table<stim::MAX_BITWORD_WIDTH>& observable_flips,
    size_t num_shots,
    size_t num_detectors
) {
    uint64_t num_mistakes = 0;

    // Convert detection events to sparse format and decode
    for (size_t shot_idx = 0; shot_idx < num_shots; shot_idx++) {
        // Build sparse shot (list of triggered detector indices)
        std::vector<uint64_t> hits;
        for (size_t det_idx = 0; det_idx < num_detectors; det_idx++) {
            if (detection_events[det_idx][shot_idx]) {
                hits.push_back(det_idx);
            }
        }

        // Decode the detection events to predict observable flips
        auto decoded_result = pm::decode_detection_events_for_up_to_64_observables(mwpm, hits, false);

        // Get actual observable flip (assuming single observable at index 0)
        // For multiple observables, need to iterate and build the mask
        uint64_t actual_obs = observable_flips[0][shot_idx] ? 1ULL : 0ULL;

        // Compare decoded vs actual observable
        if (decoded_result.obs_mask != actual_obs) {
            num_mistakes++;
        }
    }

    return num_mistakes;
}

} // namespace bucket_sim
