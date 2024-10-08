#pragma once

#include "signal_processing.h"
#include "packet.h"
#include "burst.h"
#include "swath.h"

CF_VEC_1D pulse_compression(
    const CF_VEC_1D& signal,
    const CF_VEC_1D& replica_chirp
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

CF_VEC_2D range_compress_swath(
    const std::string& filename,
    const std::string& swath_name
);

CF_VEC_2D range_compress_swath(
    std::ifstream&     data,
    const std::string& swath_name
);

CF_VEC_2D range_compress_burst(
    const std::string& filename,
    const std::string& swath,
    const int&         burst_num
);

CF_VEC_2D range_compress_burst(
    std::ifstream&     data,
    const std::string& swath,
    const int&         burst_num
);


CF_VEC_2D range_doppler_swath(
    const std::string& filename,
    const std::string& swath_name
);


CF_VEC_2D range_doppler_swath(
    std::ifstream&     data,
    const std::string& swath_name
);


CF_VEC_2D range_doppler_burst(
    const std::string& filename,
    const std::string& swath_name,
    const int&         burst_num
);


CF_VEC_2D range_doppler_burst(
    std::ifstream&     data,
    const std::string& swath_name,
    const int&         burst_num
);