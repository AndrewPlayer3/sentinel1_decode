#include "plot.hpp"

using namespace std;


const complex<float> I(0.0f, 1.0f); 

const float PI             = 3.14159235;
const float SPEED_OF_LIGHT = 299792458.0;


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


vector<complex<float>> get_replica_chirp(L0Packet& packet)
{
    vector<complex<float>> complex_samples = packet.get_complex_samples();    

    int num_samples = complex_samples.size();

    float txpsf = packet.get_start_frequency();
    float txprr = packet.get_tx_ramp_rate();
    float txpl  = packet.get_pulse_length();
    float phi_1 = txpsf - (txprr * (-0.5 * txpl));
    float phi_2 = txprr / 2;

    float range_start = -0.5 * txpl;
    float range_end   =  0.5 * txpl;
    float delta       = txpl / num_samples;

    vector<float> time(num_samples);
    vector<complex<float>> chirp(num_samples);

    for (int i = 0; i < num_samples; i++)
    {
        if (i == 0) time[i] = range_start;
        else time[i] = time[i-1] + delta;
    }

    for (int i = 0; i < num_samples; i++)
    {
        float t  = time[i]; 
        chirp[i] = float(1.0 / num_samples) * exp(I * 2.0f * PI * ((phi_1 * t) + phi_2 * (t * t)));
    }

    return chirp;
}


void plot_pulse(
    const string& filename,
    const int&    packet_index,
    const string& scaling_mode
) {
    ifstream data   = open_file(filename);
    L0Packet packet = L0Packet::get_packets(data, packet_index + 1)[packet_index];

    vector<complex<float>> complex_samples = packet.get_complex_samples();
    vector<complex<float>> replica_chirp   = get_replica_chirp(packet);

    plot_complex_samples(replica_chirp, scaling_mode);
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
    vector<complex<float>> replica_chirp   = get_replica_chirp(packet);

    vector<complex<float>> pulse_compressed(complex_samples.size());

    int num_samples = complex_samples.size();

    if (do_fft)
    {
        cout << "Doing FFT" << endl;
        
        vector<complex<float>> complex_samples_fft = compute_1d_dft(complex_samples,  0, false);
        vector<complex<float>> replica_chirp_fft   = compute_1d_dft(replica_chirp,    0, false);

        for (int i = 0; i < num_samples; i++)
        {
            pulse_compressed[i] = complex_samples_fft[i] * replica_chirp_fft[i];
        }

        plot_complex_samples(compute_1d_dft(pulse_compressed, 0, true), scaling_mode);

        return;
    }

    for (int i = 0; i < num_samples; i++)
    {
        pulse_compressed[i] = complex_samples[i] * replica_chirp[i];
    }

    plot_complex_samples(pulse_compressed, scaling_mode);
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
          int     fft_rows,
          int     fft_cols,
    const bool&   inverse,
    const string& scaling_mode
) {
    vector<vector<complex<float>>> complex_samples = decode_swath(filename, swath);

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
    const int&    axis,
          int     fft_size,
    const bool&   inverse,
    const string& scaling_mode
) {
    vector<vector<complex<float>>> complex_samples = decode_swath(filename, swath);

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