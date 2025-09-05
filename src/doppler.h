#pragma once

#include <numeric>
#include <set>

#include "packet.h"
#include "signal_processing.h"
#include "image_formation.h"
#include "misc_types.h"
#include "state_vectors.h"
#include "structs.h"

D_VEC_1D get_doppler_centroid(
    CF_VEC_2D& range_compressed,
    const double& doppler_centroid_rate,
    const double& burst_duration,
    L0Packet& first_packet,
    const int& num_rng_blocks = 15
);

double get_doppler_centroid_rate(
    PACKET_VEC_1D& burst_packets,
    const double& velocity
);
