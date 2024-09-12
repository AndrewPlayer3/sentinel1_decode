#include "image_write.h"

void _write_tif(
    std::vector<float>& img_data,
    const int& rows,
    const int& cols,
    const std::string out_filename
) {
    TIFF* tif = TIFFOpen(out_filename.c_str(), "w8");
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


void write_tif(
    CF_VEC_2D& img_data,
    const std::string& out_filename,
    const std::string& scaling_mode
) {
    int rows = img_data.size();
    int cols = img_data[0].size();

    F_VEC_1D scaled = scale(img_data, scaling_mode);

    _write_tif(scaled, rows, cols, out_filename);
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
    write_tif(signals, out_filename, scaling_mode);
}


void write_swath(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const std::string& scaling_mode
) {
    Swath swath(in_filename, swath_name);
    CF_VEC_2D signals = swath.get_all_signals();
    write_tif(signals, out_filename, scaling_mode);
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
    write_tif(chirps, out_filename, scaling_mode);
}


void write_swath_replica_chirps(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const std::string& scaling_mode
) {
    Swath swath(in_filename, swath_name);
    CF_VEC_2D chirps = swath.get_all_replica_chirps();
    write_tif(chirps, out_filename, scaling_mode);
}


void write_range_compressed_swath(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const std::string& scaling_mode
) {
    CF_VEC_2D compressed_swath = range_compress_swath(in_filename, swath_name);
    write_tif(compressed_swath, out_filename, scaling_mode);
}


void write_range_compressed_burst(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const int&         burst_num,
    const std::string& scaling_mode
) {
    CF_VEC_2D compressed_burst = range_compress_burst(in_filename, swath_name, burst_num);
    write_tif(compressed_burst, out_filename, scaling_mode);
}


void write_range_doppler_swath(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const std::string& scaling_mode
) {
    CF_VEC_2D range_doppler = range_doppler_swath(in_filename, swath_name);
    write_tif(range_doppler, out_filename, scaling_mode);
}


void write_range_doppler_burst(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const int&         burst_num,
    const std::string& scaling_mode
) {
    
    CF_VEC_2D range_doppler = range_doppler_burst(in_filename, swath_name, burst_num);
    write_tif(range_doppler, out_filename, scaling_mode);
}


void write_azimuth_compressed_burst(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const int&         burst_num,
    const std::string& scaling_mode
) {
    std::cout << "Decoding Burst" << std::endl;
    Burst burst(in_filename, swath_name, burst_num);

    std::cout << "Range Compressing Burst" << std::endl;
    CF_VEC_2D range_compressed = range_compress_burst(burst);

    PACKET_VEC_1D packets = burst.get_packets();

    std::cout << "Azimuth Compressing Burst" << std::endl;
    CF_VEC_2D azimuth_compressed = azimuth_compress(packets, range_compressed);

    std::cout << "Calling write_tif" << std::endl;
    write_tif(azimuth_compressed, out_filename, scaling_mode);
}


void write_azimuth_compressed_swath(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const std::string& scaling_mode
) {
    std::cout << "Azimuth Compressing Swath" << std::endl;
    CF_VEC_2D azimuth_compressed = azimuth_compress_swath(in_filename, swath_name);

    std::cout << "Calling write_tif" <<  std::endl;
    write_tif(azimuth_compressed, out_filename, scaling_mode);
}