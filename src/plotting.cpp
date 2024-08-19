/*
By: Andrew Player
Name: main.cpp
Description: Main function with random command-line arguments for testing purposes. 
             See "SAR Space Packet Protocol Data Unit", for more information on the packet specification:
             https://sentinels.copernicus.eu/documents/247904/2142675/Sentinel-1-SAR-Space-Packet-Protocol-Data-Unit.pdf
             For additional information on Level-0 product decoding, see:
             https://sentinel.esa.int/documents/247904/0/Sentinel-1-Level-0-Data-Decoding-Package.pdf/a8742c59-4914-40c4-8309-c77515649f17
*/

#include <thread>

#include "../include/matplotlibcpp.h"

#include "packet_decoding.hpp"

using namespace std;

namespace plt = matplotlibcpp;


void plot_complex_samples(string filename, const int& packet_index, const bool& plot_real = true)
{
    ifstream data = open_file(filename);

    vector<L0Packet> packets = get_n_packets(data, packet_index + 1, false, 0);

    vector<complex<float>> complex_samples = packets[packet_index].get_complex_samples();

    vector<float> real_parts;
    vector<float> imag_parts;

    for (complex<float> sample : complex_samples)
    {
        real_parts.push_back(sample.real());
        imag_parts.push_back(sample.imag());
    }

    plt::figure();
    if (plot_real) plt::plot(real_parts);
    else           plt::plot(imag_parts);
    plt::show();
}


void plot_swath(string filename, const int& swath_number)
{
    ifstream data = open_file(filename);

    vector<L0Packet> packets = get_all_packets(data, false, 0);

    vector<vector<complex<float>>> complex_samples;

    #pragma omp parallel for
    for (L0Packet packet : packets)
    {
        if (packet.secondary_header("swath_number") == swath_number)
        {
            complex_samples.push_back(packet.get_complex_samples());
        }
    }

    vector<float> real_samples;
    
    real_samples.reserve(complex_samples.size() * complex_samples[0].size());

    for (int i = 0; i < complex_samples.size(); i++)
    {
        
        vector<float> real_values = {};
        for (complex<float> complex_value : complex_samples[i])
        {
            real_samples.push_back(complex_value.real());
        }
    }

    const float* real_samples_ptr = &(real_samples[0]);

    int rows = complex_samples.size();
    int cols = complex_samples[0].size();

    plt::figure();
    plt::imshow(real_samples_ptr, rows, cols, 1);
    plt::show();
}


int main(int argc, char* argv[]) 
{
    if(argv[1] == __null) 
    {
        cout << "Please enter a command." << endl;
        return 1;
    }

    string command = string(argv[1]);

    if (command == "plot_complex_samples")
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
            cout << "Please enter the packet index and filename." << endl;
            return 1;
        }

        plot_swath(string(argv[3]), stoi(argv[2]));
    }
    else
    {
        cout << command << " is not a valid command." << endl;
    }

    return 0;
}