#include <set>

// #include <matplot/matplot.h>
// #include <opencv2/opencv.hpp>
#include <opencv2/highgui.hpp>

#include "packet.h"
#include "aux_decoding.h"
#include "misc_types.h"
#include "signal_processing.h"
#include "image_formation.h"
#include "cli.h"
#include "s1_decoder.h"

// #include "../include/matplotlibcpp.h"


void plot_pulse(
    const std::string& filename,
    const int&         packet_index,
    const std::string& scaling_mode
);

void plot_pulse_compression(
    const std::string& filename,
    const int&         packet_index,
    const std::string& scaling_mode
);

void plot_pulse_image(
    const std::string& filename,
    const std::string& swath,
    const int&         burst_num,
    const std::string& scaling_mode
);

void plot_range_compressed_burst(
    const std::string& filename,
    const std::string& swath,
    const int&         burst_num,
    const std::string& scaling_mode
);

void plot_range_compressed_swath(
    const std::string& filename,
    const std::string& swath,
    const std::string& scaling_mode
);

void plot_azimuth_compressed_burst(
    const std::string& filename,
    const std::string& swath_name,
    const int&         burst_num,
    const std::string& scaling_mode
);

void plot_azimuth_compressed_swath(
    const std::string& filename,
    const std::string& swath_name,
    const std::string& scaling_mode
);

void plot_burst(
    const std::string& filename,
    const std::string& swath,
    const int&         burst_num,
    const std::string& scaling_mode
);

void plot_swath(
    const std::string& filename,
    const std::string& swath,
    const std::string& scaling_mode
);

void plot_complex_image(
    const CF_VEC_2D&   signals,
    const std::string& scaling_mode
);

void plot_fft2d(
    const std::string& filename,
    const std::string& swath,
    const int&         burst_num,
          int          fft_rows,
          int          fft_cols,
    const bool&        inverse,
    const std::string& scaling_mode
);

void plot_fft_axis(
    const std::string& filename,
    const std::string& swath,
    const int&         burst_num,
    const int&         axis,
          int          fft_size,
    const bool&        inverse,
    const std::string& scaling_mode
);

void plot_signal(
    const CF_VEC_1D&   signal,
    const std::string& scaling_mode
);

void plot_signal(
    const std::string& filename,
    const int&         packet_index,
    const std::string& scaling_mode
);
