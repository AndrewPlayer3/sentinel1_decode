#include "plot.h"


F_VEC_1D scale(const CF_VEC_1D& signal, const std::string& scaling_mode)
{
    F_VEC_1D samples(signal.size());

    if      (scaling_mode == "norm_log") samples = norm_1d(signal, true);
    else if (scaling_mode == "norm"    ) samples = norm_1d(signal, false);
    else if (scaling_mode == "mag"     ) samples = magnitude_1d(signal);      
    else if (scaling_mode == "real" or scaling_mode == "imag")
    {
        bool real = scaling_mode == "real";
        
        for (int i = 0; i < signal.size(); i++)
        {
            samples[i] = real ? signal[i].real() : signal[i].imag();
        }
    }
    else
    {
        throw std::invalid_argument(scaling_mode + " is not a valid scaling mode.");
    }
    return samples;
}


F_VEC_1D scale(const CF_VEC_2D& signal, const std::string& scaling_mode)
{
    int rows = signal.size();
    int cols = signal[0].size();
    
    F_VEC_1D samples(rows*cols);

    if      (scaling_mode == "norm_log") samples = flatten(norm_2d(signal, true));
    else if (scaling_mode == "norm"    ) samples = flatten(norm_2d(signal, false));
    else if (scaling_mode == "mag"     ) samples = flatten(magnitude_2d(signal));    
    else if (scaling_mode == "real" or scaling_mode == "imag")
    {
        bool real = scaling_mode == "real";
        int  size = signal.size() * signal[0].size();
        
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                samples[i*cols+j] = real ? signal[i][j].real() : signal[i][j].imag();
            }
        }
    }
    else
    {
        throw std::invalid_argument(scaling_mode + " is not a valid scaling mode.");
    }
    return samples;
}


SIGNAL_PAIR get_signal_pairs_from_swath(
    const std::string& filename,
    const std::string& swath_name
) {
    Swath swath(filename, swath_name);

    SIGNAL_PAIR signal_pair;
    signal_pair.signals = std::move(swath.get_all_signals());
    signal_pair.replica_chirps = std::move(swath.get_all_replica_chirps());

    return signal_pair;
}


void plot_pulse(
    const std::string& filename,
    const int&         packet_index,
    const std::string& scaling_mode
) {
    std::ifstream data   = open_file(filename);
    L0Packet packet = L0Packet::get_packets(data, packet_index + 1)[packet_index];

    CF_VEC_1D signal        = packet.get_signal();
    CF_VEC_1D replica_chirp = packet.get_replica_chirp();

    plot_signal(replica_chirp, scaling_mode);
}


void plot_pulse_compression(
    const std::string& filename,
    const int&         packet_index,
    const bool&        do_fft,
    const std::string& scaling_mode
) {
    std::ifstream data = open_file(filename);
    L0Packet packet    = L0Packet::get_packets(data, packet_index + 1)[packet_index];

    CF_VEC_1D signal        = packet.get_signal();
    CF_VEC_1D replica_chirp = packet.get_replica_chirp();

    CF_VEC_1D pulse_compressed(signal.size());

    int num_samples = signal.size();

    if (do_fft)
    {
        CF_VEC_1D pulse_compressed = pulse_compression(signal, replica_chirp);
        plot_signal(pulse_compressed, scaling_mode);
        return;
    }
    for (int i = 0; i < num_samples; i++)
    {
        pulse_compressed[i] = signal[i] * replica_chirp[i];
    }
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
    std::cout << "Decoding Burst" << std::endl;

    Burst burst(filename, swath, burst_num);

    int num_packets = burst.get_num_packets();

    std::cout << "Decoded " << num_packets << " signals and replica chirps." << std::endl;

    CF_VEC_2D pulse_compressed(num_packets);

    for (int i = 0; i < num_packets; i++)
    {
        pulse_compressed[i] = pulse_compression(burst.get_signal(i), burst.get_replica_chirp(i));
    }
    std::cout << "Pulse compression completed." << std::endl;

    plot_complex_image(pulse_compressed, scaling_mode);
}


void plot_range_compressed_swath(
    const std::string& filename,
    const std::string& swath_name,
    const std::string& scaling_mode
) {
    std::cout << "Decoding Complex Samples and Replica Chirps" << std::endl;

    SIGNAL_PAIR signal_pair = get_signal_pairs_from_swath(filename, swath_name);

    int num_signals = signal_pair.signals.size();

    std::chrono::time_point compression_start = std::chrono::high_resolution_clock::now();

    CF_VEC_2D pulse_compressed(num_signals);

    for (int i = 0; i < num_signals; i++)
    {
        pulse_compressed[i] = pulse_compression(signal_pair.signals[i], signal_pair.replica_chirps[i]);
    }

    std::chrono::time_point compression_end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<float> compression_time = compression_end - compression_start;

    std::cout << "Pulse compression completed in " << compression_time.count() << "s." << std::endl;

    plot_complex_image(pulse_compressed, scaling_mode);
}


void plot_complex_image(
    const CF_VEC_2D&   signal,
    const std::string& scaling_mode
) {
    std::cout << "Flattening Vector for Plotting" << std::endl;

    int rows = signal.size();
    int cols = signal[0].size();

    F_VEC_1D samples = scale(signal, scaling_mode);

    std::cout << "Calling Plot" << std::endl;

    matplotlibcpp::figure();
    matplotlibcpp::imshow(&samples[0], rows, cols, 1);
    matplotlibcpp::show();
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

    CF_VEC_2D signals_fft = compute_axis_dft(
        burst.get_signals(),
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
    Burst burst(filename, swath, burst_num);
    plot_complex_image(burst.get_signals(), scaling_mode);
}


void plot_swath(
    const std::string& filename,
    const std::string& swath_name,
    const std::string& scaling_mode
) {
    Swath swath(filename, swath_name);
    plot_complex_image(swath.get_all_signals(), scaling_mode);
}


void plot_signal(
    const CF_VEC_1D&   signal,
    const std::string& scaling_mode
) {

    std::cout << "Creating Real-Valued Vector for Plotting" << std::endl;

    F_VEC_1D scaled_signal = scale(signal, scaling_mode);

    std::cout << "Plotting" << std::endl;

    matplotlibcpp::figure();
    matplotlibcpp::plot(scaled_signal);
    matplotlibcpp::show();
}


void plot_signal(
    const std::string& filename,
    const int&    packet_index,
    const std::string& scaling_mode
) {
    PACKET_VEC_1D packets = L0Packet::get_packets(filename, packet_index + 1);
    plot_signal(packets[packet_index].get_signal(), scaling_mode);
}

