#include "image_write.h"

void write_tif(
    std::vector<float>& img_data,
    const int& rows,
    const int& cols,
    const std::string out_filename
) {
    TIFF* tif = TIFFOpen(out_filename.c_str(), "w");
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH, rows);
    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, cols);
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 32);
    TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
    TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
    
    int img_row = 0; 
    for (int data_row = rows-1; data_row >= 0; data_row--)
    {
        TIFFWriteScanline(tif, &(img_data[data_row * cols]), img_row, 0);
        img_row++;
    }
    TIFFClose(tif);
}


void write_burst(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const int&         burst_num,
    const std::string& scaling_mode
) {
    Burst burst(in_filename, swath_name, burst_num);

    CF_VEC_2D signals = burst.get_signals();

    int rows = signals.size();
    int cols = signals[0].size();

    F_VEC_1D scaled = scale(signals, scaling_mode);

    write_tif(scaled, rows, cols, out_filename);
}


void write_swath(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const std::string& scaling_mode
) {
    Swath swath(in_filename, swath_name);

    CF_VEC_2D signals = swath.get_all_signals();

    int rows = signals.size();
    int cols = signals[0].size();

    F_VEC_1D scaled = scale(signals, scaling_mode);

    write_tif(scaled, rows, cols, out_filename);
}


void write_burst_replica_chirps(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const int&         burst_num,
    const std::string& scaling_mode
) {
    Burst burst(in_filename, swath_name, burst_num);

    CF_VEC_2D chirps = burst.get_replica_chirps();

    int rows = chirps.size();
    int cols = chirps[0].size();

    F_VEC_1D scaled = scale(chirps, scaling_mode);

    write_tif(scaled, rows, cols, out_filename);
}


void write_swath_replica_chirps(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const std::string& scaling_mode
) {
    Swath swath(in_filename, swath_name);

    CF_VEC_2D chirps = swath.get_all_replica_chirps();

    int rows = chirps.size();
    int cols = chirps[0].size();

    F_VEC_1D scaled = scale(chirps, scaling_mode);

    write_tif(scaled, rows, cols, out_filename);
}


void write_range_compressed_swath(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const std::string& scaling_mode
) {
    
    CF_VEC_2D compressed_swath = range_compress_swath(in_filename, swath_name);

    int rows = compressed_swath.size();
    int cols = compressed_swath[0].size();

    F_VEC_1D scaled = scale(compressed_swath, scaling_mode);

    write_tif(scaled, rows, cols, out_filename);
}


void write_range_compressed_burst(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const int&         burst_num,
    const std::string& scaling_mode
) {
    CF_VEC_2D compressed_burst = range_compress_burst(in_filename, swath_name, burst_num);

    int rows = compressed_burst.size();
    int cols = compressed_burst[0].size();

    F_VEC_1D scaled = scale(compressed_burst, scaling_mode);

    write_tif(scaled, rows, cols, out_filename);
}
