#pragma once

#include <iostream>
#include <vector>
#include <complex>

#include "structs.h"
#include "misc_types.h"
#include "omp.h"
#include "fftw3.h"


F_VEC_1D flatten(const F_VEC_2D& values);

CF_VEC_1D flatten(const CF_VEC_2D& values);

CF_VEC_1D pulse_compression(
    const CF_VEC_1D& signal,
    const CF_VEC_1D& replica_chirp
);

CF_VEC_1D get_reference_function(
    const CF_VEC_1D& replica_chirp
);

CF_VEC_1D conjugate(
    const CF_VEC_1D& complex_samples
);

F_VEC_1D hanning_window(const int& num_samples);

void apply_hanning_window_in_place(
    CF_VEC_1D& complex_samples
);

CF_VEC_1D apply_hanning_window(
    const CF_VEC_1D& complex_samples
);

F_VEC_1D norm_1d(
    const CF_VEC_1D& complex_values,
    const bool& log_scale
);

F_VEC_2D norm_2d(
    const CF_VEC_2D& complex_values,
    const bool& log_scale
);

F_VEC_1D magnitude_1d(
    const CF_VEC_1D& complex_values
);

F_VEC_2D magnitude_2d(
    const CF_VEC_2D& complex_values
);

CF_VEC_1D compute_1d_dft(
    const CF_VEC_1D& complex_signal,
          int   fft_size,
    const bool& inverse
);

CF_VEC_2D compute_axis_dft(
    const CF_VEC_2D& signals,
          int   fft_size,
    const int&  axis,
    const bool& inverse
);

CF_VEC_2D _compute_axis_dft(
    const CF_VEC_2D& signals,
          int   fft_size,
    const bool& inverse
);

CF_VEC_2D compute_2d_dft(
    const CF_VEC_2D& complex_samples,
    const bool& inverse,
    int fft_rows,
    int fft_cols
);