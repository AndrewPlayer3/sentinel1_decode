#pragma once

#include <iostream>
#include <vector>
#include <complex>

#include "omp.h"
#include "fftw3.h"

using namespace std;



template <typename T>
vector<T> flatten(const vector<vector<T>>& values);

vector<float> norm_1d(const vector<complex<float>>& complex_values, const bool& log_scale);

vector<vector<float>> norm_2d(const vector<vector<complex<float>>>& complex_values, const bool& log_scale);

vector<float> magnitude_1d(const vector<complex<float>>& complex_values);

vector<vector<float>> magnitude_2d(const vector<vector<complex<float>>>& complex_values);

vector<complex<float>> compute_1d_dft(const vector<complex<float>>& complex_signal);

vector<vector<complex<float>>> compute_2d_dft(
    const vector<vector<complex<float>>>& complex_samples,
    int fft_rows,
    int fft_cols,
    const bool& inverse
);


template <typename T>
vector<T> flatten(const vector<vector<T>>& values)
{
    int rows = values.size();
    int cols = values[0].size();

    vector<T> flat(rows * cols);

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            flat[i * cols + j] = values[i][j];
        }
    }

    return flat;
}