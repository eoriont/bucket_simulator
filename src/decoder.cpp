#include "decoder.hpp"
#include <exception>

namespace bucket_sim {

uint64_t decode_batch(
    pm::Mwpm& mwpm,
    const stim::simd_bit_table<stim::MAX_BITWORD_WIDTH>& detection_events,
    const stim::simd_bit_table<stim::MAX_BITWORD_WIDTH>& observable_flips,
    size_t num_shots,
    size_t num_detectors,
    size_t num_observables
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

        // Decode the detection events to predict observable flips.
        // Some distributed circuits can yield odd-parity components that are
        // undecodable by strict MWPM. Treat these shots as decoding failures
        // (counted as mistakes) instead of aborting the full run.
        uint64_t decoded_obs_mask = 0;
        bool decode_ok = true;
        try {
            auto decoded_result = pm::decode_detection_events_for_up_to_64_observables(mwpm, hits, false);
            decoded_obs_mask = decoded_result.obs_mask;
        } catch (const std::exception&) {
            decode_ok = false;
        }

        // Build actual observable mask from all observables
        uint64_t actual_obs = 0;
        for (size_t i = 0; i < num_observables; i++) {
            if (observable_flips[i][shot_idx]) actual_obs |= (1ULL << i);
        }

        // Compare decoded vs actual observable
        if (!decode_ok || decoded_obs_mask != actual_obs) {
            num_mistakes++;
        }
    }

    return num_mistakes;
}

} // namespace bucket_sim
