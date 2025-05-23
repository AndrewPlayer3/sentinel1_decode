#include "plot.h"


void plot_range_compressed_burst(
    const std::string& filename,
    const std::string& swath,
    const int&         burst_num,
    const std::string& scaling_mode
) {
    S1_Decoder s1(filename);

    plot_complex_image(s1.get_range_compressed_burst(swath, burst_num), scaling_mode);
}


void plot_range_compressed_swath(
    const std::string& filename,
    const std::string& swath_name,
    const std::string& scaling_mode
) {
    S1_Decoder s1(filename);

    plot_complex_image(s1.get_range_compressed_swath(swath_name), scaling_mode);
}


void plot_azimuth_compressed_burst(
    const std::string& filename,
    const std::string& swath_name,
    const int&         burst_num,
    const std::string& scaling_mode
) {
    S1_Decoder s1(filename); 

    plot_complex_image(s1.get_azimuth_compressed_burst(swath_name, burst_num), scaling_mode);
}


void plot_azimuth_compressed_swath(
    const std::string& filename,
    const std::string& swath_name,
    const std::string& scaling_mode
) {
    S1_Decoder s1(filename);

    plot_complex_image(s1.get_azimuth_compressed_swath(swath_name), scaling_mode);
}

cv::Mat g_original_img;
float g_auto_min = 0.0f;
float g_auto_max = 1.0f;

// Trackbar callback
void on_trackbar(int, void*) {
    int min_slider = cv::getTrackbarPos("Min", "Sentinel-1 Decoder");
    int max_slider = cv::getTrackbarPos("Max", "Sentinel-1 Decoder");

    float min_val = g_auto_min + min_slider / 1000.0f * (g_auto_max - g_auto_min);
    float max_val = g_auto_min + max_slider / 1000.0f * (g_auto_max - g_auto_min);

    if (max_val <= min_val) return;

    cv::Mat rescaled;
    g_original_img.convertTo(rescaled, CV_32F, 1.0f / (max_val - min_val), -min_val / (max_val - min_val));
    cv::threshold(rescaled, rescaled, 0.0, 0.0, cv::THRESH_TOZERO);
    cv::threshold(rescaled, rescaled, 1.0, 1.0, cv::THRESH_TRUNC);
    rescaled.convertTo(rescaled, CV_8U, 255.0);

    cv::imshow("Sentinel-1 Decoder", rescaled);
}

void plot_complex_image(
    const CF_VEC_2D& signal,
    const std::string& scaling_mode
) {
    std::cout << "Flattening Vector for Plotting" << std::endl;

    int rows = signal.size();
    int cols = signal[0].size();

    std::vector<float> samples = scale(signal, scaling_mode);
    g_original_img = cv::Mat(rows, cols, CV_32F, samples.data()).clone();

    double min_val, max_val;
    cv::minMaxLoc(g_original_img, &min_val, &max_val);
    g_auto_min = static_cast<float>(min_val);
    g_auto_max = static_cast<float>(max_val);

    std::cout << "Calling Plot" << std::endl;

    cv::namedWindow("Sentinel-1 Decoder", cv::WINDOW_GUI_EXPANDED);
    cv::createTrackbar("Min", "Sentinel-1 Decoder", nullptr, 1000, on_trackbar);
    cv::createTrackbar("Max", "Sentinel-1 Decoder", nullptr, 1000, on_trackbar);
    cv::setTrackbarPos("Min", "Sentinel-1 Decoder", 500);
    cv::setTrackbarPos("Max", "Sentinel-1 Decoder", 1000);

    on_trackbar(0, nullptr);
    cv::waitKey(0);
}


// void plot_complex_image(
//     const CF_VEC_2D&   signal,
//     const std::string& scaling_mode
// ) {
//     std::cout << "Flattening Vector for Plotting" << std::endl;

//     int rows = signal.size();
//     int cols = signal[0].size();

//     std::vector<float> samples = scale(signal, scaling_mode);

//     cv::Mat img(rows, cols, CV_32F, samples.data());

//     std::cout << "Calling Plot" << std::endl;

//     cv::namedWindow("Sentinel-1 Decoder", cv::WINDOW_GUI_EXPANDED);
//     cv::setWindowProperty("Sentinel-1 Decoder", cv::WINDOW_OPENGL, cv::WINDOW_GUI_EXPANDED);
//     cv::imshow("Sentinel-1 Decoder", img);
//     cv::waitKey();
// }


void plot_fft2d(
    const std::string& filename,
    const std::string& swath,
    const int&         burst_num,
          int          fft_rows,
          int          fft_cols,
    const bool&        inverse,
    const std::string& scaling_mode
) {
    S1_Decoder s1(filename);

    CF_VEC_2D signals = s1.get_burst(swath, burst_num);

    if (not fft_rows) fft_rows = signals.size();
    if (not fft_cols) fft_cols = signals[0].size();

    CF_VEC_2D signals_fft = compute_2d_dft(
        signals,
        inverse,
        fft_rows,
        fft_cols
    );

    plot_complex_image(signals_fft, scaling_mode);
}


void plot_fft_axis(
    const std::string& filename,
    const std::string& swath,
    const int&         burst_num,
    const int&         axis,
          int          fft_size,
    const bool&        inverse,
    const std::string& scaling_mode
) {
    S1_Decoder s1(filename);

    CF_VEC_2D signals = s1.get_burst(swath, burst_num);

    int fft_rows = signals.size();
    int fft_cols = signals[0].size();

    CF_VEC_2D signals_fft = compute_axis_dft(
        signals,
        fft_size,
        axis,
        inverse
    );

    int out_rows = signals_fft.size();
    int out_cols = signals_fft[0].size();

    std::cout << "Input Rows: "  << fft_rows  << " Input Cols: "  << fft_cols << std::endl;
    std::cout << "Output Rows: " << out_rows  << " Output Cols: " << out_cols << std::endl;

    plot_complex_image(signals_fft, scaling_mode);
}


void plot_burst(
    const std::string& filename,
    const std::string& swath,
    const int&         burst_num,
    const std::string& scaling_mode
) {
    S1_Decoder s1(filename);

    plot_complex_image(s1.get_burst(swath, burst_num), scaling_mode);
}


void plot_swath(
    const std::string& filename,
    const std::string& swath_name,
    const std::string& scaling_mode
) {
    S1_Decoder s1(filename);

    CF_VEC_2D swath = s1.get_swath(swath_name);

    plot_complex_image(swath, scaling_mode);
}
