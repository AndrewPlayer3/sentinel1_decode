#pragma once

#include "tiff.h"
#include "tiffio.h"
#include "signal_processing.h"
#include "image_formation.h"
#include "packet.h"
#include "burst.h"
#include "swath.h"

void _write_tif(
    std::vector<double>& img_data,
    const int& rows,
    const int& cols,
    const std::string out_filename
);

void write_tif(
    CF_VEC_2D& img_data,
    const std::string& out_filename,
    const std::string& scaling_mode
);

void write_burst(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const int&         burst_num,
    const std::string& scaling_mode
);

void write_swath(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const std::string& scaling_mode
);

void write_burst_replica_chirps(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const int&         burst_num,
    const std::string& scaling_mode
);

void write_swath_replica_chirps(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const std::string& scaling_mode
);

void write_range_compressed_swath(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const std::string& scaling_mode
);

void write_range_compressed_burst(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const int&         burst_num,
    const std::string& scaling_mode
);

void write_range_doppler_swath(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const std::string& scaling_mode
);

void write_range_doppler_burst(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const int&         burst_num,
    const std::string& scaling_mode
);

void write_azimuth_compressed_burst(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const int&         burst_num,
    const std::string& scaling_mode
);

void write_azimuth_compressed_swath(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const std::string& scaling_mode
);