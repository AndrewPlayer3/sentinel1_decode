#include "plot.hpp"

using namespace std;



pair<vector<vector<complex<float>>>, vector<L0Packet>> decode_bursts(const string& filename, const string& swath, const int& burst_num)
{
    cout << "Parsing Packets" << endl;
    
    ifstream data = open_file(filename);
    vector<L0Packet> packets = L0Packet::get_packets_in_swath(data, swath);
    int num_packets = packets.size();

    cout << "Found " << num_packets << " packets in " << swath << "." << endl;

    vector<vector<L0Packet>> bursts; 
    vector<L0Packet> burst_packets;
    set<int> burst_nums;
    int previous_az = 0;
    int num_bursts = 0;

    for (int i = 0; i < num_packets; i++)
    {
        L0Packet packet = packets[i];
        if (packet.get_data_format() == 'D')
        {
            int az = packet.secondary_header("azimuth_beam_address");
            if (i == 0) previous_az = az;
            if (az != previous_az and az != previous_az + 1)
            {
                bursts.push_back(burst_packets);
                burst_packets = vector<L0Packet>();
                burst_nums.emplace(num_bursts++);
            }
            burst_packets.push_back(packet);
            previous_az = az;
        }
        if (i == num_packets - 1) 
        {
            bursts.push_back(burst_packets);
        }
    }

    vector<L0Packet> burst = bursts[burst_num];

    num_packets = burst.size();
    if (num_packets == 0)
    {
        cout << "No packets found for burst #" << burst_num << endl;
        cout << "Available burst numbers are: " << endl;
        for (int num : burst_nums) cout << num << endl;
        exit(0);
    }
    cout << "Decoding " << num_packets << " Complex Samples for Burst #" << burst_num << endl;

    vector<vector<complex<float>>> complex_samples(num_packets);

    #pragma omp parallel for
    for (int i = 0; i < num_packets; i++)
    {
        L0Packet packet = burst[i];
        complex_samples[i] = packet.get_complex_samples();
    }
    if (complex_samples.size() < 1)
    {
        throw runtime_error("No samples found for swath " + swath);
    }
    return pair(complex_samples, burst);
}



vector<vector<complex<float>>> decode_swath(const string& filename, const string& swath)
{
    cout << "Parsing Packets" << endl;
    
    ifstream data = open_file(filename);

    vector<L0Packet> packets = L0Packet::get_packets_in_swath(data, swath);

    int num_packets = packets.size();

    cout << "Found " << num_packets << " packets in " << swath << "." << endl;

    vector<vector<complex<float>>> complex_samples(num_packets);

    cout << "Decoding Complex Samples" << endl;

    #pragma omp parallel for
    for (int i = 0; i < num_packets; i++)
    {
        L0Packet packet = packets[i];
        complex_samples[i] = packet.get_complex_samples();
    }

    if (complex_samples.size() < 1)
    {
        throw runtime_error("No samples found for swath " + swath);
    }

    return complex_samples;
}


vector<float> scale(const vector<complex<float>>& complex_samples, const string& scaling_mode)
{
    vector<float> samples(complex_samples.size());

    if      (scaling_mode == "norm_log") samples = norm_1d(complex_samples, true);
    else if (scaling_mode == "norm"    ) samples = norm_1d(complex_samples, false);
    else if (scaling_mode == "mag"     ) samples = magnitude_1d(complex_samples);      
    else if (scaling_mode == "real" or scaling_mode == "imag")
    {
        bool real = scaling_mode == "real";
        
        for (int i = 0; i < complex_samples.size(); i++)
        {
            samples[i] = real ? complex_samples[i].real() : complex_samples[i].imag();
        }
    }
    else
    {
        throw invalid_argument(scaling_mode + " is not a valid scaling mode.");
    }
    return samples;
}


vector<float> scale(const vector<vector<complex<float>>>& complex_samples, const string& scaling_mode)
{
    int rows = complex_samples.size();
    int cols = complex_samples[0].size();
    
    vector<float> samples(rows*cols);

    if      (scaling_mode == "norm_log") samples = flatten(norm_2d(complex_samples, true));
    else if (scaling_mode == "norm"    ) samples = flatten(norm_2d(complex_samples, false));
    else if (scaling_mode == "mag"     ) samples = flatten(magnitude_2d(complex_samples));    
    else if (scaling_mode == "real" or scaling_mode == "imag")
    {
        bool real = scaling_mode == "real";
        int  size = complex_samples.size() * complex_samples[0].size();
        
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                samples[i*cols+j] = real ? complex_samples[i][j].real() : complex_samples[i][j].imag();
            }
        }
    }
    else
    {
        throw invalid_argument(scaling_mode + " is not a valid scaling mode.");
    }
    return samples;
}


void plot_pulse(
    const string& filename,
    const int&    packet_index,
    const string& scaling_mode
) {
    ifstream data   = open_file(filename);
    L0Packet packet = L0Packet::get_packets(data, packet_index + 1)[packet_index];

    vector<complex<float>> complex_samples = packet.get_complex_samples();
    vector<complex<float>> replica_chirp   = packet.get_replica_chirp();

    plot_complex_samples(replica_chirp, scaling_mode);
}


vector<complex<float>> process_replica(const vector<complex<float>> replica_chirp)
{
    int num_samples = replica_chirp.size();
    vector<complex<float>> replica_chirp_fft   = replica_chirp;// compute_1d_dft(replica_chirp, 0, false);
    vector<complex<float>> replica_chirp_conj = conjugate(replica_chirp_fft);
    vector<complex<float>> weighted_replica = replica_chirp_conj;
    vector<float> norm = magnitude_1d(weighted_replica);
    float energy = 0.0;
    for (int i = 0; i < num_samples; i++)
    {
        energy += (norm[i] * norm[i]);
    }
    energy /= norm.size();
    vector<complex<float>> replica_out(num_samples);
    for (int i = 0; i < num_samples; i++)
    {
        replica_out[i] = weighted_replica[i] / energy;
    }
    return replica_out;
}


vector<complex<float>> pulse_compression(
    const vector<complex<float>>& complex_samples,
    const vector<complex<float>>& replica_chirp
) {
    vector<complex<float>> pulse_compressed(complex_samples.size());

    int num_samples = complex_samples.size();

    vector<complex<float>> complex_samples_fft = compute_1d_dft(complex_samples,  0, false);
    vector<complex<float>> reference = process_replica(replica_chirp);

    for (int i = 0; i < num_samples; i++)
    {
        pulse_compressed[i] = complex_samples_fft[i] * reference[i];
    }    

    return compute_1d_dft(pulse_compressed, 0, true);
}


void plot_pulse_compression(
    const string& filename,
    const int&    packet_index,
    const bool&   do_fft,
    const string& scaling_mode
) {
    ifstream data   = open_file(filename);
    L0Packet packet = L0Packet::get_packets(data, packet_index + 1)[packet_index];

    vector<complex<float>> complex_samples = packet.get_complex_samples();
    vector<complex<float>> replica_chirp   = packet.get_replica_chirp();

    vector<complex<float>> pulse_compressed(complex_samples.size());

    int num_samples = complex_samples.size();

    if (do_fft)
    {
        vector<complex<float>> pulse_compressed = pulse_compression(complex_samples, replica_chirp);

        plot_complex_samples(pulse_compressed, scaling_mode);

        return;
    }

    for (int i = 0; i < num_samples; i++)
    {
        pulse_compressed[i] = complex_samples[i] * replica_chirp[i];
    }

    plot_complex_samples(pulse_compressed, scaling_mode);
}


void plot_pulse_image(
    const string& filename,
    const string& swath,
    const int&    burst_num,
    const string& scaling_mode
) {

    auto samples_packets = decode_bursts(filename, swath, burst_num);

    vector<vector<complex<float>>> complex_samples = samples_packets.first;
    vector<L0Packet> packets = samples_packets.second;

    int num_packets = packets.size();

    cout << "Found " << num_packets << " packets in " << swath << "." << endl;

    vector<vector<complex<float>>> replica_chirps(num_packets);

    cout << "Decoding Complex Samples and Replica Chirps" << endl;

    auto decoding_start = chrono::high_resolution_clock::now();

    #pragma omp parallel for
    for (int i = 0; i < num_packets; i++)
    {
        L0Packet packet = packets[i];
        replica_chirps[i]  = packet.get_replica_chirp();
    }

    if (replica_chirps.size() < 1)
    {
        throw runtime_error("No samples found for swath " + swath);
    }

    auto decoding_end = chrono::high_resolution_clock::now();

    chrono::duration<float> decode_time = decoding_end - decoding_start;

    cout << "Decoded " << replica_chirps.size() << " replica chirps in " << decode_time.count() << "s" << endl;

    plot_complex_image(replica_chirps, scaling_mode);
}



void plot_pulse_compressed_image(
    const string& filename,
    const string& swath,
    const int&    burst_num,
    const string& scaling_mode
) {

    auto samples_packets = decode_bursts(filename, swath, burst_num);

    vector<vector<complex<float>>> complex_samples = samples_packets.first;
    vector<L0Packet> packets = samples_packets.second;

    int num_packets = packets.size();

    cout << "Found " << num_packets << " packets in " << swath << "." << endl;

    vector<vector<complex<float>>> replica_chirps(num_packets);

    cout << "Decoding Complex Samples and Replica Chirps" << endl;

    auto decoding_start = chrono::high_resolution_clock::now();

    #pragma omp parallel for
    for (int i = 0; i < num_packets; i++)
    {
        L0Packet packet = packets[i];
        replica_chirps[i]  = packet.get_replica_chirp();
    }

    if (complex_samples.size() < 1)
    {
        throw runtime_error("No samples found for swath " + swath);
    }

    auto decoding_end = chrono::high_resolution_clock::now();

    chrono::duration<float> decode_time = decoding_end - decoding_start;

    cout << "Decoded " << complex_samples.size() << " complex sample vectors, and " 
         << replica_chirps.size() << " replica chirps, in " << decode_time.count() << "s" << endl;

    vector<vector<complex<float>>> pulse_compressed(num_packets);

    auto compression_start = chrono::high_resolution_clock::now();

    for (int i = 0; i < num_packets; i++)
    {
        pulse_compressed[i] = pulse_compression(complex_samples[i], replica_chirps[i]);
    }

    // vector<vector<complex<float>>> azimuth = compute_axis_dft(pulse_compressed, 0, 0, false);

    // Not proper yet, just for interesting visualization.
    // vector<vector<complex<float>>> azimuth = compute_2d_dft(pulse_compressed, true, 0, 0);

    // azimuth = compute_axis_dft(azimuth, 0, 0, false);

    auto compression_end = chrono::high_resolution_clock::now();

    chrono::duration<float> compression_time = compression_end - compression_start;

    cout << "Pulse compression completed in " << compression_time.count() << endl;

    plot_complex_image(pulse_compressed, scaling_mode);
}


void plot_complex_image(
    const vector<vector<complex<float>>>& complex_samples,
    const string& scaling_mode
) {
    cout << "Flattening Vector for Plotting" << endl;

    int rows = complex_samples.size();
    int cols = complex_samples[0].size();

    vector<float> samples = scale(complex_samples, scaling_mode);

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

    vector<complex<float>> complex_samples     = packets[packet_index].get_complex_samples();
    vector<complex<float>> complex_samples_fft = compute_1d_dft(complex_samples, fft_size, inverse);

    plot_complex_samples(complex_samples_fft, scaling_mode);
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
    auto samples_packets = decode_bursts(filename, swath, burst_num);

    vector<vector<complex<float>>> complex_samples = samples_packets.first;
    vector<L0Packet> packets = samples_packets.second;

    if (not fft_rows) fft_rows = complex_samples.size();
    if (not fft_cols) fft_cols = complex_samples[0].size();

    vector<vector<complex<float>>> complex_samples_fft = compute_2d_dft(
        complex_samples,
        inverse,
        fft_rows,
        fft_cols
    );

    plot_complex_image(complex_samples_fft, scaling_mode);
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
    auto samples_packets = decode_bursts(filename, swath, burst_num);

    vector<vector<complex<float>>> complex_samples = samples_packets.first;
    vector<L0Packet> packets = samples_packets.second;

    int fft_rows = complex_samples.size();
    int fft_cols = complex_samples[0].size();

    vector<vector<complex<float>>> complex_samples_fft = compute_axis_dft(
        complex_samples,
        fft_size,
        axis,
        inverse
    );

    int out_rows = complex_samples_fft.size();
    int out_cols = complex_samples_fft[0].size();

    cout << "Input Rows: "  << fft_rows  << " Input Cols: "  << fft_cols << endl;
    cout << "Output Rows: " << out_rows  << " Output Cols: " << out_cols << endl;

    plot_complex_image(complex_samples_fft, scaling_mode);
}


void plot_burst(
    const string& filename,
    const string& swath,
    const int&    burst_num,
    const string& scaling_mode
) {
    auto samples_packets = decode_bursts(filename, swath, burst_num);

    vector<vector<complex<float>>> complex_samples = samples_packets.first;

    plot_complex_image(complex_samples, scaling_mode);
}


void plot_swath(
    const string& filename,
    const string& swath,
    const string& scaling_mode
) {
    vector<vector<complex<float>>> complex_samples = decode_swath(filename, swath);

    plot_complex_image(complex_samples, scaling_mode);
}


void plot_complex_samples(
    const vector<complex<float>>& complex_samples,
    const string& scaling_mode
) {

    cout << "Creating Real-Valued Vector for Plotting" << endl;

    vector<float> samples = scale(complex_samples, scaling_mode);

    cout << "Plotting" << endl;

    matplotlibcpp::figure();
    matplotlibcpp::plot(samples);
    matplotlibcpp::show();
}


void plot_complex_samples(
    const string& filename,
    const int&    packet_index,
    const string& scaling_mode
) {
    ifstream data = open_file(filename);

    vector<L0Packet> packets = L0Packet::get_packets(data, packet_index + 1);

    vector<complex<float>> complex_samples = packets[packet_index].get_complex_samples();

    plot_complex_samples(complex_samples, scaling_mode);
}