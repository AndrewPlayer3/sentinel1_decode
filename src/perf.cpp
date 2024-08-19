/*
By: Andrew Player
Name: perf.cpp
Description: Functions for evaluating the performance of packet decoding. 
             See "SAR Space Packet Protocol Data Unit", for more information on the packet specification:
             https://sentinels.copernicus.eu/documents/247904/2142675/Sentinel-1-SAR-Space-Packet-Protocol-Data-Unit.pdf
             For additional information on Level-0 product decoding, see:
             https://sentinel.esa.int/documents/247904/0/Sentinel-1-Level-0-Data-Decoding-Package.pdf/a8742c59-4914-40c4-8309-c77515649f17
*/

#include "perf.hpp"


/* Returns the time in seconds that it takes to decode the complex samples of num_packets packets */
double time_packet_generation(ifstream& data, const int& num_packets, const bool& log, const int& log_interval)
{
    double total_runtime = 0.0;
    for (int i = 0; i < num_packets; i++)
    {
        try
        {
            auto start = chrono::high_resolution_clock::now();

            L0Packet packet = L0Packet::get_next_packet(data);

            vector<complex<float>> complex_samples = packet.get_complex_samples();

            auto end   = chrono::high_resolution_clock::now();

            chrono::duration<double> difference = end - start;
            total_runtime += difference.count();

            if (log && i % log_interval == 0 && i != 0)
            {
                cout << "Decoded " << i << " packets in " << total_runtime << "s." << endl;
            }

            if (data.eof())
            {
                if (log) cout << "EOF reached after decoding " << i << " packets." << endl;
                break;
            }
        }
        catch(runtime_error)
        {
            cout << "Caught runtime error while decoding packet #" << i << endl;
            continue;
        }
        catch(length_error)
        {
            break;
        }
    }
    if (log) cout << "Decoded " << num_packets << " packets in " << total_runtime << "s." << endl;

    return total_runtime;
}


void omp_test(string filename)
{
    ifstream data = open_file(filename);

    vector<L0Packet> packets = L0Packet::get_packets(data);

    const int num_packets = packets.size();

    chrono::time_point start = chrono::high_resolution_clock::now();

    vector<vector<complex<float>>> complex_samples(num_packets);

    #pragma omp parallel for
    for (int i = 0; i < num_packets; i++)
    {
        complex_samples[i] = packets[i].get_complex_samples();
    }

    chrono::duration<double> difference = chrono::high_resolution_clock::now() - start;

    cout << "Decoded " << num_packets << " packets in " << difference.count() << "s." << endl;
}


void thread_runner(
    vector<vector<complex<float>>>& complex_samples,
    vector<L0Packet>& packets,
    const int start_index,
    const int end_index
)
{
    for (int i = start_index; i < end_index; i++)
    {
        complex_samples[i] = packets[i].get_complex_samples();
    }
}


void thread_test(string filename)
{
    ifstream data = open_file(filename);
    
    vector<L0Packet> packets = L0Packet::get_packets(data);
    vector<thread>   threads;

    const int num_packets = packets.size();
    const int num_threads = thread::hardware_concurrency();
    const int chunk_size  = num_packets / num_threads;

    chrono::time_point start = chrono::high_resolution_clock::now();

    vector<vector<complex<float>>> complex_samples(num_packets);

    for (int i = 0; i < num_threads; i++)
    {
        int start_index =  i * chunk_size;
        int end_index   = (i == num_threads - 1) ? num_packets : start_index + chunk_size;
        
        threads.emplace_back(
            thread_runner,
            ref(complex_samples),
            ref(packets),
            start_index,
            end_index
        );
    }

    for (thread& thread : threads)
    {
        if (thread.joinable()) thread.join();
    }

    chrono::duration<double> difference = chrono::high_resolution_clock::now() - start;

    cout << "Decoded " << num_packets << " packets in " << difference.count() << "s." << endl;
}
