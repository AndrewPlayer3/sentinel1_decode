/*
By: Andrew Player
Name: perf.cpp
Description: Functions for evaluating the performance of packet decoding. 
             See "SAR Space Packet Protocol Data Unit", for more information on the packet specification:
             https://sentinels.copernicus.eu/documents/247904/2142675/Sentinel-1-SAR-Space-Packet-Protocol-Data-Unit.pdf
             For additional information on Level-0 product decoding, see:
             https://sentinel.esa.int/documents/247904/0/Sentinel-1-Level-0-Data-Decoding-Package.pdf/a8742c59-4914-40c4-8309-c77515649f17
*/

#include "perf.h"


/* Returns the time in seconds that it takes to decode the complex samples of num_packets packets */
double time_packet_generation(
    const std::string& filename, 
    const int&  num_packets, 
    const bool& log, 
    const int&  log_interval
) {
    std::ifstream data = open_file(filename);
    double total_runtime = 0.0;
    for (int i = 0; i < num_packets; i++)
    {
        try
        {
            std::chrono::time_point start = std::chrono::high_resolution_clock::now();
            L0Packet::get_next_packet(data).get_signal();
            std::chrono::time_point end   = std::chrono::high_resolution_clock::now();

            std::chrono::duration<double> difference = end - start;
            total_runtime += difference.count();

            if (log && i % log_interval == 0 && i != 0)
            {
                std::cout << "Decoded " << i << " packets in " << total_runtime << "s." << std::endl;
            }
            if (data.eof())
            {
                if (log) std::cout << "EOF reached after decoding " << i << " packets." << std::endl;
                break;
            }
        }
        catch(std::runtime_error)
        {
            std::cout << "Caught runtime error while decoding packet #" << i << std::endl;
            continue;
        }
        catch(std::length_error) { break; }
    }
    if (log) std::cout << "Decoded " << num_packets << " packets in " << total_runtime << "s." << std::endl;

    return total_runtime;
}


void omp_test(const std::string& filename)
{
    PACKET_VEC_1D packets = L0Packet::get_packets(filename);

    std::chrono::time_point start = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for
    for (int i = 0; i < packets.size(); i++)
    {
        packets[i].get_signal();
    }

    std::chrono::duration<double> difference = std::chrono::high_resolution_clock::now() - start;

    std::cout << "Decoded " << packets.size() << " packets in " << difference.count() << "s." << std::endl;
}


void thread_runner(
    PACKET_VEC_1D& packets,
    const int start_index,
    const int end_index
) {
    for (int i = start_index; i < end_index; i++)
    {
        packets[i].get_signal();
    }
}


void thread_test(const std::string& filename)
{ 
    PACKET_VEC_1D packets = L0Packet::get_packets(filename);
    std::vector<std::thread> threads;

    const int num_packets = packets.size();
    const int num_threads = std::thread::hardware_concurrency();
    const int chunk_size  = num_packets / num_threads;

    std::chrono::time_point start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < num_threads; i++)
    {
        int start_index =  i * chunk_size;
        int end_index   = (i == num_threads - 1) ? num_packets : start_index + chunk_size;

        threads.emplace_back(
            thread_runner,
            ref(packets),
            start_index,
            end_index
        );
    }
    for (std::thread& thread : threads)
    {
        if (thread.joinable()) thread.join();
    }
    std::chrono::duration<double> difference = std::chrono::high_resolution_clock::now() - start;

    std::cout << "Decoded " << num_packets << " packets in " << difference.count() << "s." << std::endl;
}


/* Returns the time in seconds that it takes to decode the complex samples of num_packets packets */
double time_range_compression(
    const std::string& filename, 
    const std::string& swath_name,
    const int&  burst_num, 
    const bool& log, 
    const int&  log_interval
) {

    Burst burst(filename, swath_name, burst_num);

    std::chrono::time_point start = std::chrono::high_resolution_clock::now();
    
    CF_VEC_2D compressed_burst = range_compress_burst(burst);

    std::chrono::time_point end   = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> difference = end - start;

    return difference.count();
}


/* Returns the time in seconds that it takes to decode the complex samples of num_packets packets */
double time_azimuth_compression(
    const std::string& filename, 
    const std::string& swath_name,
    const int&  burst_num, 
    const bool& log, 
    const int&  log_interval
) {

    Burst burst(filename, swath_name, burst_num);

    PACKET_VEC_1D packets = burst.get_packets();
    CF_VEC_2D signals = burst.get_signals();

    std::chrono::time_point start = std::chrono::high_resolution_clock::now();

    CF_VEC_2D compressed_burst = azimuth_compress(packets, signals);

    std::chrono::time_point end   = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> difference = end - start;

    return difference.count();
}