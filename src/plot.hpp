#include <set>

#include "packet.hpp"
#include "signal_processing.hpp"

#include "../include/matplotlibcpp.h"

using namespace std;


vector<vector<complex<float>>> decode_swath(
    const string& filename,
    const string& swath
);

void plot_pulse(
    const string& filename,
    const int&    packet_index,
    const string& scaling_mode
);

void plot_pulse_compression(
    const string& filename,
    const int&    packet_index,
    const bool&   do_fft,
    const string& scaling_mode
);

void plot_pulse_image(
    const string& filename,
    const string& swath,
    const int&    burst_num,
    const string& scaling_mode
);

void plot_pulse_compressed_image(
    const string& filename,
    const string& swath,
    const int&    burst_num,
    const string& scaling_mode
);

void plot_burst(
    const string& filename,
    const string& swath,
    const int&    burst_num,
    const string& scaling_mode
);

void plot_swath(
    const string& filename,
    const string& swath,
    const string& scaling_mode
);

void plot_fft(
    const string& filename,
    const int&    packet_index,
    const int&    fft_size,
    const bool&   inverse,
    const string& scaling_mode
);

void plot_complex_image(
    const vector<vector<complex<float>>>& complex_samples,
    const string& scaling_mode
);

void plot_fft2d(
    const string& filename,
    const string& swath,
    const int&    burst_num,
          int     fft_rows,
          int     fft_cols,
    const bool&   inverse,
    const string& scaling_mode
);

void plot_fft_axis(
    const string& filename,
    const string& swath,
    const int&    burst_num,
    const int&    axis,
          int     fft_size,
    const bool&   inverse,
    const string& scaling_mode
);

void plot_complex_samples(
    const vector<complex<float>>& complex_samples,
    const string& scaling_mode
);

void plot_complex_samples(
    const string& filename,
    const int&    packet_index,
    const string& scaling_mode
);