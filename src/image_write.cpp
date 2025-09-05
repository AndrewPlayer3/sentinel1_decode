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

    std::vector<float> scaled = scale(img_data, scaling_mode);

    _write_tif(scaled, rows, cols, out_filename);
}


void write_spectrogram(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const int&         burst_num,
    const int&         azimuth_line,
    const std::string& scaling_mode
) {
    S1_Decoder s1(in_filename);

    CF_VEC_2D burst = s1.get_burst(swath_name, burst_num);

    if (azimuth_line < 0 || azimuth_line > burst.size())
    {
        std::cout << "Azimuth Line " << azimuth_line 
                  << " is outside of the valid range of 0 to " << burst.size() << std::endl;
    }

    CF_VEC_2D specgram = spectrogram(burst[azimuth_line], 64, 32);

    write_tif(specgram, out_filename, scaling_mode);
}


void write_burst(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const int&         burst_num,
    const std::string& scaling_mode
) {
    S1_Decoder s1(in_filename);

    CF_VEC_2D burst = s1.get_burst(swath_name, burst_num);

    write_tif(burst, out_filename, scaling_mode);
}


void write_swath(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const std::string& scaling_mode
) {
    S1_Decoder s1(in_filename);

    CF_VEC_2D swath = s1.get_swath(swath_name);

    write_tif(swath, out_filename, scaling_mode);
}


void write_range_compressed_swath(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const std::string& scaling_mode
) {
    S1_Decoder s1(in_filename);
    CF_VEC_2D range_compressed_swath = s1.get_range_compressed_swath(swath_name);
    write_tif(range_compressed_swath, out_filename, scaling_mode);
}


void write_range_compressed_burst(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const int&         burst_num,
    const std::string& scaling_mode
) {
    S1_Decoder s1(in_filename);
    CF_VEC_2D range_compressed_burst = s1.get_range_compressed_burst(swath_name, burst_num);
    write_tif(range_compressed_burst, out_filename, scaling_mode);
}


void write_range_doppler_burst(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const int&         burst_num,
    const std::string& scaling_mode
) {
    S1_Decoder s1(in_filename);
    CF_VEC_2D range_compressed_burst = s1.get_range_compressed_burst(swath_name, burst_num, true);
    write_tif(range_compressed_burst, out_filename, scaling_mode);
}


void write_range_doppler_swath(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const std::string& scaling_mode
) {
    S1_Decoder s1(in_filename);
    CF_VEC_2D range_compressed_swath = s1.get_range_compressed_swath(swath_name, true);
    write_tif(range_compressed_swath, out_filename, scaling_mode);
}


void write_azimuth_compressed_burst(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const int&         burst_num,
    const std::string& scaling_mode
) {
    S1_Decoder s1(in_filename);
    CF_VEC_2D azimuth_compressed_burst = s1.get_azimuth_compressed_burst(swath_name, burst_num);
    write_tif(azimuth_compressed_burst, out_filename, scaling_mode);
}


void write_azimuth_compressed_swath(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const std::string& scaling_mode
) {
    S1_Decoder s1(in_filename);
    CF_VEC_2D azimuth_compressed_swath = s1.get_azimuth_compressed_swath(swath_name);
    write_tif(azimuth_compressed_swath, out_filename, scaling_mode);
}


void write_azimuth_compressed_burst_eccm(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const int&         burst_num,
    const int&         detection_threshold,
    const int&         mitigation_threshold,
    const std::string& scaling_mode
) {
    S1_Decoder s1(in_filename);
    s1.apply_eccm(detection_threshold, mitigation_threshold);
    CF_VEC_2D azimuth_compressed_burst = s1.get_azimuth_compressed_burst(swath_name, burst_num);
    write_tif(azimuth_compressed_burst, out_filename, scaling_mode);
}


void write_azimuth_compressed_swath_eccm(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const int&         detection_threshold,
    const int&         mitigation_threshold,
    const std::string& scaling_mode
) {
    S1_Decoder s1(in_filename);
    s1.apply_eccm(detection_threshold, mitigation_threshold);
    CF_VEC_2D azimuth_compressed_swath = s1.get_azimuth_compressed_swath(swath_name);
    write_tif(azimuth_compressed_swath, out_filename, scaling_mode);
}
