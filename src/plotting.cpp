/*
By: Andrew Player
Name: main.cpp
Description: Main function with random command-line arguments for testing purposes. 
             See "SAR Space Packet Protocol Data Unit", for more information on the packet specification:
             https://sentinels.copernicus.eu/documents/247904/2142675/Sentinel-1-SAR-Space-Packet-Protocol-Data-Unit.pdf
             For additional information on Level-0 product decoding, see:
             https://sentinel.esa.int/documents/247904/0/Sentinel-1-Level-0-Data-Decoding-Package.pdf/a8742c59-4914-40c4-8309-c77515649f17
*/

#include "packet.hpp"
#include "signal_processing.hpp"

#include "../include/matplotlibcpp.h"

using namespace std;


void plot_complex_samples(const vector<complex<float>>& complex_samples)
{
    vector<float> norm = norm_1d(complex_samples, true);

    matplotlibcpp::figure();
    matplotlibcpp::plot(norm);
    matplotlibcpp::show();
}


void plot_complex_samples(string filename, const int& packet_index)
{
    ifstream data = open_file(filename);

    vector<L0Packet> packets = L0Packet::get_packets(data, packet_index + 1);

    vector<complex<float>> complex_samples = packets[packet_index].get_complex_samples();

    plot_complex_samples(complex_samples);
}


void plot_complex_image(const vector<vector<complex<float>>>& complex_samples)
{
    cout << "Setting Up Normalized Real-Value Vector for Plotting" << endl;

    float max_value = 0;    

    int rows = complex_samples.size();
    int cols = complex_samples[0].size();

    vector<float> flat = flatten(norm_2d(complex_samples, true));

    cout << "Calling Plot" << endl;

    matplotlibcpp::figure();
    matplotlibcpp::imshow(&flat[0], rows, cols, 1);
    matplotlibcpp::show();
}


void plot_fft(string filename, const int& packet_index)
{
    ifstream data = open_file(filename);

    vector<L0Packet> packets = L0Packet::get_packets(data, packet_index + 1);

    vector<complex<float>> complex_samples = packets[packet_index].get_complex_samples();

    vector<complex<float>> complex_samples_fft = compute_1d_dft(complex_samples);

    plot_complex_samples(complex_samples_fft);
}


void plot_fft2d(string filename, const string& swath, const bool& plot_real = true)
{
    cout << "Parsing Packets" << endl;
    
    vector<L0Packet> packets = L0Packet::get_packets(filename);

    vector<vector<complex<float>>> complex_samples;

    cout << "Decoding Complex Samples" << endl;

    #pragma omp parallel for
    for (L0Packet packet : packets)
    {
        if (packet.get_swath() == swath)
        {
            complex_samples.push_back(packet.get_complex_samples());
        }
    }

    int fft_rows = complex_samples.size();
    int fft_cols = complex_samples[0].size();

    vector<vector<complex<float>>> complex_samples_fft = compute_2d_dft(
        complex_samples,
        fft_rows,
        fft_cols,
        false
    );

    plot_complex_image(complex_samples_fft);
}


void plot_swath(string filename, const string& swath)
{
    cout << "Parsing Packets" << endl;
    
    ifstream data = open_file(filename);

    vector<L0Packet> packets = L0Packet::get_packets(data);

    cout << "Decoding Complex Samples" << endl;

    vector<vector<complex<float>>> complex_samples;

    #pragma omp parallel for
    for (L0Packet packet : packets)
    {
        if (packet.get_swath() == swath)
        {
            complex_samples.push_back(packet.get_complex_samples());
        }
    }

    if (complex_samples.size() < 1)
    {
        cout << "No samples found for swath " << swath << endl;
        return;
    }

    plot_complex_image(complex_samples);
}


int main(int argc, char* argv[]) 
{
    if(argv[1] == __null) 
    {
        cout << "Please enter a command." << endl;
        return 1;
    }

    fftwf_init_threads();

    fftwf_plan_with_nthreads(omp_get_max_threads());

    string command = string(argv[1]);

    if (command == "help" or command == "--help" or command == "-h")
    {
        cout << "plot_complex_samples [packet_index] [path]" << endl;
        cout << "plot_swath [swath: S1, IW1, etc...] [path]" << endl;
        cout << "plot_fft [packet_index] [path]"             << endl;
        cout << "plot_fft2 [swath: S1, IW1, etc...] [path]"  << endl;
    }

    else if (command == "plot_complex_samples")
    {
        if(argv[2] == __null || argv[3] == __null) 
        {
            cout << "Please enter the packet index and filename." << endl;
            return 1;
        }

        plot_complex_samples(string(argv[3]), stoi(argv[2]));
    }

    else if (command == "plot_swath")
    {
        if(argv[2] == __null || argv[3] == __null) 
        {
            cout << "Please enter the swath and filename." << endl;
            return 1;
        }

        plot_swath(string(argv[3]), string(argv[2]));
    }

    else if (command == "plot_fft")
    {
        if(argv[2] == __null || argv[3] == __null) 
        {
            cout << "Please enter the packet index and filename." << endl;
            return 1;
        }

        plot_fft(string(argv[3]), stoi(argv[2]));
    }

    else if (command == "plot_fft2")
    {
        if(argv[2] == __null || argv[3] == __null) 
        {
            cout << "Please enter the swath and the filename." << endl;
            return 1;
        }

        plot_fft2d(string(argv[3]), string(argv[2]));
    }

    else
    {
        cout << command << " is not a valid command." << endl;
    }

    fftwf_cleanup_threads();

    return 0;
}