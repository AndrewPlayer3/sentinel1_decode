#pragma once

#include <math.h>
#include <algorithm>
#include <numeric>
#include <execution>
#include <ctime>

#include "signal_processing.h"
#include "aux_decoding.h"
#include "packet.h"
#include "state_vectors.h"

CF_VEC_1D pulse_compression(
    const CF_VEC_1D& signal,
    const CF_VEC_1D& reference
);

CF_VEC_1D get_reference_function(
    const CF_VEC_1D& replica_chirp
);

std::pair<PACKET_VEC_2D, int> get_azimuth_blocks(
    PACKET_VEC_1D& packets
);

CF_VEC_2D azimuth_frequency_ufr(
    CF_VEC_2D& range_compressed,
    F_VEC_1D&  dc_estimates,
    L0Packet&  initial_packet,
    const double& dc_rate,
    const double& burst_duration,
    const double& prf,
    const double& processing_bandwidth  // B_d
);

CF_VEC_2D azimuth_time_ufr(
    CF_VEC_2D& range_compressed,
    F_VEC_1D&  dc_estimates,
    F_VEC_2D&  az_fm_rate,
    L0Packet&  initial_packet,
    const double& dc_rate,
    const double& burst_duration,
    const double& prf,
    const double& processing_bandwidth,  // B_d
    const int&    swath_number
);

F_VEC_1D get_effective_velocities(
    const F_VEC_1D& position,
    const double& velocity,
    const F_VEC_1D& slant_ranges
);

F_VEC_1D apply_src_and_rcmc(
    CF_VEC_1D& range_line,
    const F_VEC_1D& effective_velocities,
    const F_VEC_1D& slant_ranges,
    const F_VEC_1D& range_freqs,
    const F_VEC_1D& doppler_centroids,
    const double& az_freq
);