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

#include "packet_decoding.hpp"

using namespace std;


void omp_test(string filename)
{
    ifstream data = open_file(filename);

    vector<L0Packet> packets = get_all_packets(data, false, 10);

    const int num_packets = packets.size();

    chrono::time_point start = chrono::high_resolution_clock::now();

    vector<vector<complex<double>>> complex_samples(num_packets);

    #pragma omp parallel for
    for (int i = 0; i < num_packets; i++)
    {
        complex_samples[i] = packets[i].get_complex_samples();
    }

    chrono::duration<double> difference = chrono::high_resolution_clock::now() - start;

    cout << "Decoded " << num_packets << " packets in " << difference.count() << "s." << endl;
}


void thread_runner(
    vector<vector<complex<double>>>& complex_samples,
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
    
    vector<L0Packet> packets = get_all_packets(data, false, 0);
    vector<thread>   threads;

    const int num_packets = packets.size();
    const int num_threads = thread::hardware_concurrency();
    const int chunk_size  = num_packets / num_threads;

    chrono::time_point start = chrono::high_resolution_clock::now();

    vector<vector<complex<double>>> complex_samples(num_packets);

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


void print_packet_at_index(
    string filename,
    int index,
    bool headers = true,
    bool pulse_info = false,
    bool modes = false
) {
    ifstream data = open_file(filename);
    
    vector<L0Packet> packets = get_n_packets(data, index + 1, false, 0);

    if (index >= packets.size())
    {
        throw out_of_range("Index is greater than the number of packets.");
    }
    L0Packet packet = packets[index];

    if (headers)
    {
        packet.print_primary_header();
        packet.print_secondary_header();
    }
    if (pulse_info) packet.print_pulse_info();
    if (modes)      packet.print_modes();
}


int main(int argc, char* argv[]) 
{
    if(argv[1] == __null) 
    {
        cout << "Please enter a command." << endl;
        return 1;
    }

    string command = string(argv[1]);


    if (command == "print_headers")
    {
        
        if(argv[2] == __null || argv[3] == __null) 
        {
            cout << "Please enter a packet index and filename." << endl;
            return 1;
        }
        print_packet_at_index(string(argv[3]), stoi(argv[2]));
    }
    
    else if (command == "time")
    {
        
        if(argv[2] == __null || argv[3] == __null) 
        {
            cout << "Please enter the packet count to time and a filename." << endl;
            return 1;
        }
        ifstream data = open_file(string(argv[3]));

        double runtime = time_packet_generation(data, stoi(argv[2]), false, 0);

        cout << "Decoded " << stoi(argv[2]) << " packets in " << runtime << "s." << endl;
    }

    else if (command == "print_pulse_info")
    {

        if(argv[2] == __null || argv[3] == __null) 
        {
            cout << "Please enter the packet index and filename." << endl;
            return 1;
        }
        print_packet_at_index(string(argv[3]), stoi(argv[2]), false, true);
    }
    
    else if (command == "print_modes")
    {
        
        if(argv[2] == __null || argv[3] == __null) 
        {
            cout << "Please enter the packet index and filename." << endl;
            return 1;
        }
        print_packet_at_index(string(argv[3]), stoi(argv[2]), false, false, true);
    }
    
    else if (command == "thread_test")
    {
        
        if(argv[2] == __null) 
        {
            cout << "Please enter the filename." << endl;
            return 1;
        }
        thread_test(string(argv[2]));
    }
    
    else if (command == "omp_test")
    {
        
        if(argv[2] == __null) 
        {
            cout << "Please enter the filename." << endl;
            return 1;
        }
        omp_test(string(argv[2]));
    }

    else if (command == "print_complex_samples")
    {

        if(argv[2] == __null || argv[3] == __null) 
        {
            cout << "Please enter the packet index and filename." << endl;
            return 1;
        }
        ifstream data = open_file(string(argv[3]));

        vector<L0Packet> packets = get_n_packets(data, stoi(argv[2]) + 1, false, 0);

        vector<complex<double>> complex_samples = packets[stoi(argv[2])].get_complex_samples();

        for (complex<double> sample : complex_samples)
        {
            cout << sample << endl;
        }
    }

    else if (command == "find_packets_of_type")
    {

        if(argv[2] == __null || argv[3] == __null) 
        {
            cout << "Please enter the packet type and filename." << endl;
            return 1;
        }
        ifstream data = open_file(string(argv[3]));

        char type = char(argv[2][0]);

        vector<L0Packet> packets = get_all_packets(data, false, 0);

        for (int i = 0; i < packets.size(); i++)
        {
            L0Packet packet = packets[i];
            char data_format    = packet.get_data_format();
            int  sequence_count = packet.primary_header("packet_sequence_count");
            if (data_format == type)
            {
                cout << "Packet #" << sequence_count << " at index " << i << " is type " << data_format << endl;
            }
        }
    }

    else if (command == "print_packet_type")
    {

        if(argv[2] == __null || argv[3] == __null) 
        {
            cout << "Please enter the packet index and filename." << endl;
            return 1;
        }
        ifstream data = open_file(string(argv[3]));

        vector<L0Packet> packets = get_n_packets(data, stoi(argv[2]) + 1, false, 0);

        cout << "Packet Type is " << packets[stoi(argv[2])].get_data_format() << "." << endl;
    }
    else
    {
        cout << command << " is not a valid command." << endl;
    }

    return 0;
}