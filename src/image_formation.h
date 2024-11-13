#pragma once

#include <math.h>
#include <algorithm>
#include <numeric>
#include <execution>
#include <ctime>

#include "signal_processing.h"
#include "aux_decoding.h"
#include "packet.h"
#include "burst.h"
#include "swath.h"
#include "state_vectors.h"

CF_VEC_1D pulse_compression(
    const CF_VEC_1D& signal,
    const CF_VEC_1D& reference
);

CF_VEC_1D get_reference_function(
    const CF_VEC_1D& replica_chirp
);

SIGNAL_PAIR get_signal_pairs_from_swath(
    const std::string& filename,
    const std::string& swath_name
);

SIGNAL_PAIR get_signal_pairs_from_swath(
    std::ifstream&     data,
    const std::string& swath_name
);

// CF_VEC_2D range_compress_swath(
//     const std::string& filename,
//     const std::string& swath_name
// );

// CF_VEC_2D range_compress_swath(
//     std::ifstream&     data,
//     const std::string& swath_name
// );

// CF_VEC_2D range_compress_swath(
//     Swath& swath
// );

// CF_VEC_2D range_compress_burst(
//     const std::string& filename,
//     const std::string& swath,
//     const int&         burst_num
// );

// CF_VEC_2D range_compress_burst(
//     std::ifstream&     data,
//     const std::string& swath,
//     const int&         burst_num
// );

// CF_VEC_2D range_compress_burst(
//     Burst& burst
// );

// CF_VEC_2D range_compress_burst(
//     std::ifstream&     data,
//     const std::string& swath,
//     const int&         burst_num
// );

// CF_VEC_2D range_doppler_swath(
//     const std::string& filename,
//     const std::string& swath_name
// );

// CF_VEC_2D range_doppler_swath(
//     std::ifstream&     data,
//     const std::string& swath_name
// );

void fftshift(std::vector<std::complex<double>>& data);

// CF_VEC_2D azimuth_compress_wip(
//     PACKET_VEC_1D& packets,
//     CF_VEC_2D& signals
// );

// CF_VEC_2D azimuth_compress_swath(
//     const std::string& filename,
//     const std::string& swath_name
// );

// CF_VEC_2D azimuth_compress_swath(
//    Swath& swath,
//    STATE_VECTORS& state_vectors
// );

// CF_VEC_2D azimuth_compress_burst_mvp(
//     std::ifstream&     data,
//     const std::string& swath,
//     const int&         burst_num,
//     STATE_VECTORS&     state_vectors
// );

// CF_VEC_2D azimuth_compress_burst_mvp(
//     Burst& burst,
//     STATE_VECTORS& state_vectors
// );

// CF_VEC_2D range_doppler_burst(
//     const std::string& filename,
//     const std::string& swath_name,
//     const int&         burst_num
// );

// CF_VEC_2D range_doppler_burst(
//     std::ifstream&     data,
//     const std::string& swath_name,
//     const int&         burst_num
// );
