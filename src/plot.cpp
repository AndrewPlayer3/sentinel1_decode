#include "plot.hpp"

using namespace std;

vector<float> scale(const vector<complex<float>>& signal, const string& scaling_mode)
{
    vector<float> samples(signal.size());

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
        throw invalid_argument(scaling_mode + " is not a valid scaling mode.");
    }
    return samples;
}


vector<float> scale(const vector<vector<complex<float>>>& signal, const string& scaling_mode)
{
    int rows = signal.size();
    int cols = signal[0].size();
    
    vector<float> samples(rows*cols);

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
        throw invalid_argument(scaling_mode + " is not a valid scaling mode.");
    }
    return samples;
}


struct SIGNAL_PAIR {
    vector<vector<complex<float>>> signals;
    vector<vector<complex<float>>> replica_chirps;
};


SIGNAL_PAIR get_signal_pairs_from_swath(
    const string& filename,
    const string& swath_name
) {
    Swath swath(filename, swath_name);

    SIGNAL_PAIR signal_pair;
    signal_pair.signals = std::move(swath.get_all_signals());
    signal_pair.replica_chirps = std::move(swath.get_all_replica_chirps());

    return signal_pair;
}


void plot_pulse(
    const string& filename,
    const int&    packet_index,
    const string& scaling_mode
) {
    ifstream data   = open_file(filename);
    L0Packet packet = L0Packet::get_packets(data, packet_index + 1)[packet_index];

    vector<complex<float>> signal        = packet.get_signal();
    vector<complex<float>> replica_chirp = packet.get_replica_chirp();

    plot_signal(replica_chirp, scaling_mode);
}


void plot_pulse_compression(
    const string& filename,
    const int&    packet_index,
    const bool&   do_fft,
    const string& scaling_mode
) {
    ifstream data   = open_file(filename);
    L0Packet packet = L0Packet::get_packets(data, packet_index + 1)[packet_index];

    vector<complex<float>> signal        = packet.get_signal();
    vector<complex<float>> replica_chirp = packet.get_replica_chirp();

    vector<complex<float>> pulse_compressed(signal.size());

    int num_samples = signal.size();

    if (do_fft)
    {
        vector<complex<float>> pulse_compressed = pulse_compression(signal, replica_chirp);
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
    const string& filename,
    const string& swath,
    const int&    burst_num,
    const string& scaling_mode
) {
    cout << "Decoding Burst" << endl;

    Burst burst(filename, swath, burst_num);

    cout << "Decoded " << burst.get_num_packets() 
         << " signals and replica chirps." << endl;

    plot_complex_image(burst.get_replica_chirps(), scaling_mode);
}


void plot_range_compressed_burst(
    const string& filename,
    const string& swath,
    const int&    burst_num,
    const string& scaling_mode
) {
    cout << "Decoding Burst" << endl;

    Burst burst(filename, swath, burst_num);

    int num_packets = burst.get_num_packets();

    cout << "Decoded " << num_packets << " signals and replica chirps." << endl;

    vector<vector<complex<float>>> pulse_compressed(num_packets);

    for (int i = 0; i < num_packets; i++)
    {
        pulse_compressed[i] = pulse_compression(burst.get_signal(i), burst.get_replica_chirp(i));
    }
    cout << "Pulse compression completed." << endl;

    plot_complex_image(pulse_compressed, scaling_mode);
}


void plot_range_compressed_swath(
    const string& filename,
    const string& swath_name,
    const string& scaling_mode
) {
    cout << "Decoding Complex Samples and Replica Chirps" << endl;

    SIGNAL_PAIR signal_pair = get_signal_pairs_from_swath(filename, swath_name);

    int num_signals = signal_pair.signals.size();

    auto compression_start = chrono::high_resolution_clock::now();

    vector<vector<complex<float>>> pulse_compressed(num_signals);

    for (int i = 0; i < num_signals; i++)
    {
        pulse_compressed[i] = pulse_compression(signal_pair.signals[i], signal_pair.replica_chirps[i]);
    }

    auto compression_end = chrono::high_resolution_clock::now();

    chrono::duration<float> compression_time = compression_end - compression_start;

    cout << "Pulse compression completed in " << compression_time.count() << "s." << endl;

    plot_complex_image(pulse_compressed, scaling_mode);
}


void plot_complex_image(
    const vector<vector<complex<float>>>& signal,
    const string& scaling_mode
) {
    cout << "Flattening Vector for Plotting" << endl;

    int rows = signal.size();
    int cols = signal[0].size();

    vector<float> samples = scale(signal, scaling_mode);

    cout << "Calling Plot" << endl;

    matplotlibcpp::figure();
    matplotlibcpp::imshow(&samples[0], rows, cols, 1);
    matplotlibcpp::show();
}


void plot_fft(
    const string& filename,
    const int&    packet_index,
    const int&    fft_size,
    const bool&   inverse,
    const string& scaling_mode
) {
    ifstream data = open_file(filename);

    cout << "Parsing Packets" << endl;

    vector<L0Packet> packets = L0Packet::get_packets(data, packet_index + 1);

    vector<complex<float>> signal     = packets[packet_index].get_signal();
    vector<complex<float>> signal_fft = compute_1d_dft(signal, fft_size, inverse);

    plot_signal(signal_fft, scaling_mode);
}


void plot_fft2d(
    const string& filename,
    const string& swath,
    const int&    burst_num,
          int     fft_rows,
          int     fft_cols,
    const bool&   inverse,
    const string& scaling_mode
) {
    Burst burst = Burst(filename, swath, burst_num);

    vector<vector<complex<float>>> signals = burst.get_signals();

    if (not fft_rows) fft_rows = signals.size();
    if (not fft_cols) fft_cols = signals[0].size();

    vector<vector<complex<float>>> signals_fft = compute_2d_dft(
        signals,
        inverse,
        fft_rows,
        fft_cols
    );

    plot_complex_image(signals_fft, scaling_mode);
}


void plot_fft_axis(
    const string& filename,
    const string& swath,
    const int&    burst_num,
    const int&    axis,
          int     fft_size,
    const bool&   inverse,
    const string& scaling_mode
) {
    Burst burst(filename, swath, burst_num);

    int fft_rows = burst.get_num_packets();
    int fft_cols = burst.get_signal(0).size();

    vector<vector<complex<float>>> signals_fft = compute_axis_dft(
        burst.get_signals(),
        fft_size,
        axis,
        inverse
    );

    int out_rows = signals_fft.size();
    int out_cols = signals_fft[0].size();

    cout << "Input Rows: "  << fft_rows  << " Input Cols: "  << fft_cols << endl;
    cout << "Output Rows: " << out_rows  << " Output Cols: " << out_cols << endl;

    plot_complex_image(signals_fft, scaling_mode);
}


void plot_burst(
    const string& filename,
    const string& swath,
    const int&    burst_num,
    const string& scaling_mode
) {
    Burst burst(filename, swath, burst_num);
    plot_complex_image(burst.get_signals(), scaling_mode);
}


void plot_swath(
    const string& filename,
    const string& swath_name,
    const string& scaling_mode
) {
    Swath swath(filename, swath_name);
    plot_complex_image(swath.get_all_signals(), scaling_mode);
}


void plot_signal(
    const vector<complex<float>>& signal,
    const string& scaling_mode
) {

    cout << "Creating Real-Valued Vector for Plotting" << endl;

    vector<float> scaled_signal = scale(signal, scaling_mode);

    cout << "Plotting" << endl;

    matplotlibcpp::figure();
    matplotlibcpp::plot(scaled_signal);
    matplotlibcpp::show();
}


void plot_signal(
    const string& filename,
    const int&    packet_index,
    const string& scaling_mode
) {
    vector<L0Packet> packets = L0Packet::get_packets(filename, packet_index + 1);
    plot_signal(packets[packet_index].get_signal(), scaling_mode);
}

