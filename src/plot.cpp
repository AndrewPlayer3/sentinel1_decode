#include "plot.h"


void plot_pulse(
    const std::string& filename,
    const int&         packet_index,
    const std::string& scaling_mode
) {
    std::ifstream data   = open_file(filename);
    L0Packet packet = L0Packet::get_packets(data, packet_index + 1)[packet_index];

    CF_VEC_1D replica_chirp = packet.get_replica_chirp();

    plot_signal(replica_chirp, scaling_mode);
}


void plot_pulse_compression(
    const std::string& filename,
    const int&         packet_index,
    const std::string& scaling_mode
) {
    std::ifstream data = open_file(filename);
    L0Packet packet    = L0Packet::get_packets(data, packet_index + 1)[packet_index];

    CF_VEC_1D pulse_compressed = pulse_compression(
        packet.get_signal(),
        packet.get_replica_chirp()
    );

    plot_signal(pulse_compressed, scaling_mode);
}


void plot_pulse_image(
    const std::string& filename,
    const std::string& swath,
    const int&         burst_num,
    const std::string& scaling_mode
) {
    std::cout << "Decoding Burst" << std::endl;

    Burst burst(filename, swath, burst_num);

    std::cout << "Decoded " << burst.get_num_packets() 
              << " signals and replica chirps." << std::endl;

    plot_complex_image(burst.get_replica_chirps(), scaling_mode);
}


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


void plot_complex_image(
    const CF_VEC_2D&   signal,
    const std::string& scaling_mode
) {
    std::cout << "Flattening Vector for Plotting" << std::endl;

    int rows = signal.size();
    int cols = signal[0].size();

    std::vector<float> samples = scale(signal, scaling_mode);

    cv::Mat img(rows, cols, CV_32F, samples.data());

    std::cout << "Calling Plot" << std::endl;

    cv::namedWindow("Sentinel-1 Decoder", cv::WINDOW_GUI_EXPANDED);
    cv::setWindowProperty("Sentinel-1 Decoder", cv::WINDOW_OPENGL, cv::WINDOW_GUI_EXPANDED);
    cv::imshow("Sentinel-1 Decoder", img);
    cv::waitKey();
}


void plot_fft(
    const std::string& filename,
    const int&         packet_index,
    const int&         fft_size,
    const bool&        inverse,
    const std::string& scaling_mode
) {
    std::ifstream data = open_file(filename);

    std::cout << "Parsing Packets" << std::endl;

    PACKET_VEC_1D packets = L0Packet::get_packets(data, packet_index + 1);

    CF_VEC_1D signal     = packets[packet_index].get_signal();
    CF_VEC_1D signal_fft = compute_1d_dft(signal, fft_size, inverse);

    plot_signal(signal_fft, scaling_mode);
}


void plot_fft2d(
    const std::string& filename,
    const std::string& swath,
    const int&         burst_num,
          int          fft_rows,
          int          fft_cols,
    const bool&        inverse,
    const std::string& scaling_mode
) {
    Burst burst = Burst(filename, swath, burst_num);

    CF_VEC_2D signals = burst.get_signals();

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
    Burst burst(filename, swath, burst_num);

    int fft_rows = burst.get_num_packets();
    int fft_cols = burst.get_signal(0).size();

    CF_VEC_2D signals = burst.get_signals();

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


void plot_signal(
    const CF_VEC_1D&   signal,
    const std::string& scaling_mode
) {

    std::cout << "Creating Real-Valued Vector for Plotting" << std::endl;

    F_VEC_1D scaled_signal = scale(signal, scaling_mode);

    std::cout << "Plotting" << std::endl;

    // matplotlibcpp::figure();
    // matplotlibcpp::plot(scaled_signal);
    // matplotlibcpp::show();
}


void plot_signal(
    const std::string& filename,
    const int&    packet_index,
    const std::string& scaling_mode
) {
    PACKET_VEC_1D packets = L0Packet::get_packets(filename, packet_index + 1);
    plot_signal(packets[packet_index].get_signal(), scaling_mode);
}

