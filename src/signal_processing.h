#pragma once

#include <algorithm>
#include <execution>
#include <iostream>
#include <vector>
#include <complex>
#include <numeric>

#include "structs.h"
#include "misc_types.h"
#include "omp.h"
#include "fftw3.h"


CF_VEC_1D linspace(
    const std::complex<double>& start,
    const std::complex<double>& end,
    const int& size
);

F_VEC_1D linspace(
    const double& start, 
    const double& end, 
    const int&   size
);

std::vector<float> flatten(const std::vector<std::vector<float>>& values);

F_VEC_1D flatten(
    const F_VEC_2D& values
);

CF_VEC_1D flatten(
    const CF_VEC_2D& values
);

void conjugate_in_place(
    CF_VEC_1D& complex_samples
);

void apply_hanning_window_in_place(
    CF_VEC_1D& complex_samples
);


F_VEC_1D norm_1d(
    const CF_VEC_1D& complex_values,
    const bool& log_scale
);

std::vector<std::vector<float>> norm_2d(
    const CF_VEC_2D& complex_values,
    const bool& log_scale
);

F_VEC_1D magnitude_1d(
    const CF_VEC_1D& complex_values
);

std::vector<std::vector<float>> magnitude_2d(
    const CF_VEC_2D& complex_values
);

CF_VEC_2D transpose(const CF_VEC_2D& arr);

F_VEC_1D fftfreq(
    int n,
    double d
);

CF_VEC_1D compute_1d_dft(
    const CF_VEC_1D& complex_signal,
          int   fft_size,
    const bool& inverse
);

void compute_1d_dft_in_place(
          CF_VEC_1D& signal,
          int   fft_size,
    const bool& inverse
);

void _compute_axis_dft(
    CF_VEC_2D&  signals,
          int&  fft_size,
    const bool& inverse
);

CF_VEC_2D compute_axis_dft(
    CF_VEC_2D&  signals,
          int   fft_size,
    const int&  axis,
    const bool& inverse
);

void compute_axis_dft_in_place(
    CF_VEC_2D&  signals,
          int   fft_size,
    const int&  axis,
    const bool& inverse
);

CF_VEC_2D compute_2d_dft(
    const CF_VEC_2D& complex_samples,
    const bool& inverse,
    int fft_rows,
    int fft_cols
);

F_VEC_1D scale(
    const CF_VEC_1D& signal,
    const std::string& scaling_mode
);

std::vector<float> scale(
    const CF_VEC_2D& signal,
    const std::string& scaling_mode
);

std::vector<fftw_plan> get_fftw_plans(CF_VEC_2D& signals);

void destroy_fftw_plans(std::vector<fftw_plan>& plans);

void fftshift_in_place(CF_VEC_1D& signal);
void fftshift_in_place(CF_VEC_2D& signals);