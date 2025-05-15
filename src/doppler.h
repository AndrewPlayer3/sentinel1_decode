#pragma once

#include <set>

#include "packet.h"
#include "signal_processing.h"
#include "image_formation.h"
#include "misc_types.h"
#include "state_vectors.h"
#include "structs.h"

F_VEC_1D get_doppler_centroid(CF_VEC_2D& range_compressed, const double& doppler_centroid_rate, const double& burst_duration, L0Packet& first_packet);

double get_doppler_centroid_rate(PACKET_VEC_1D& burst_packets, const double& velocity);

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
    F_VEC_2D&  ka,
    L0Packet&  initial_packet,
    const double& dc_rate,
    const double& burst_duration,
    const double& prf,
    const double& processing_bandwidth  // B_d
);