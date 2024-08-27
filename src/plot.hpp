#include <set>

#include "packet.hpp"
#include "aux_decoding.hpp"
#include "burst.hpp"
#include "swath.hpp"
#include "signal_processing.hpp"

#include "../include/matplotlibcpp.h"

using namespace std;


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

void plot_range_compressed_burst(
    const string& filename,
    const string& swath,
    const int&    burst_num,
    const string& scaling_mode
);

void plot_range_compressed_swath(
    const string& filename,
    const string& swath,
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
    const vector<vector<complex<float>>>& signals,
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

void plot_signal(
    const vector<complex<float>>& signal,
    const string& scaling_mode
);

void plot_signal(
    const string& filename,
    const int&    packet_index,
    const string& scaling_mode
);