#pragma once

#include <iostream>
#include <vector>
#include <complex>

#include "omp.h"
#include "fftw3.h"

using namespace std;


vector<float> flatten(const vector<vector<float>>& values);

vector<float> norm_1d(const vector<complex<float>>& complex_values, const bool& log_scale);

vector<vector<float>> norm_2d(const vector<vector<complex<float>>>& complex_values, const bool& log_scale);

vector<float> magnitude_1d(const vector<complex<float>>& complex_values);

vector<vector<float>> magnitude_2d(const vector<vector<complex<float>>>& complex_values);

vector<complex<float>> compute_1d_dft(const vector<complex<float>>& complex_signal);

vector<vector<complex<float>>> compute_2d_dft(const vector<vector<complex<float>>>& complex_samples);
