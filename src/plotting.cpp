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

        ifstream data = open_file(string(argv[3]));

        vector<L0Packet> packets = get_n_packets(data, stoi(argv[2]) + 1, false, 0);

        vector<complex<double>> complex_samples = packets[stoi(argv[2])].get_complex_samples();

        vector<double> real_parts;
        vector<double> imag_parts;

        for (complex<double> sample : complex_samples)
        {
            real_parts.push_back(sample.real());
            imag_parts.push_back(sample.imag());
        }

        plt::figure();
        plt::plot(real_parts);
        plt::xlabel("Time");
        plt::ylabel("Amplitude");
        plt::show();
    }
    else
    {
        cout << command << " is not a valid command." << endl;
    }

    return 0;
}