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

void fftshift(std::vector<std::complex<double>>& data);
