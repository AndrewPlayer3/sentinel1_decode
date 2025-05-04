#pragma once

#include <set>

#include "packet.h"
#include "signal_processing.h"
#include "image_formation.h"
#include "misc_types.h"
#include "state_vectors.h"
#include "structs.h"

F_VEC_1D get_doppler_centroid(CF_VEC_2D& range_compressed, const double& doppler_centroid_rate, const double& burst_duration, L0Packet& first_packet);