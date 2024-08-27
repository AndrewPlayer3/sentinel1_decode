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

#include "perf.hpp"
#include "aux_decoding.hpp"

using namespace std;


void validate_args(
    const string& command,
    const vector<string>& command_args,
    const vector<string>& arg_types,
          char*   args[],
    const int&    start_index = 2,
    const string& help = ""
) {
    int end_index = command_args.size() + start_index;
    for (int i = start_index; i < end_index; i++)
    {
        int index = i - start_index;
        if (args[i] == __null)
        {
            cout << command  << " is missing the following argument: "
                 << command_args[index] << endl;
            exit(1);
        }
        try
        {
            if      (arg_types[index] == "string") string(args[i]);
            else if (arg_types[index] == "char")   char(args[i][0]);
            else if (arg_types[index] == "int")    stoi(args[i]);
            else if (arg_types[index] == "float")  stof(args[i]);
            else if (arg_types[index] == "path")   open_file(string(args[i]));
        }
        catch(...)
        {
            cout << args[i] << " is not a valid " << command_args[index] << "." << endl;
            exit(1);
        }
    }
}


void print_index_records(const string& filename)
{
    vector<unordered_map<string, u_int64_t>> records = index_decoder(filename);

    int index = 0;
    for (unordered_map<string, u_int64_t> record : records)
    {
        cout << "Index Record #"  << index << ":"            << endl;
        cout << "Date/Time: "     << record["date_time"]     << endl;
        cout << "Time Delta: "    << record["time_delta"]    << endl;
        cout << "Data Size: "     << record["data_size"]     << endl;
        cout << "Unit Offset: "   << record["unit_offset"]   << endl;
        cout << "Byte Offset: "   << record["byte_offset"]   << endl;
        cout << "Variable Flag: " << record["variable_flag"] << endl;
        cout << "Spare Data: "    << record["spare_data"]    << endl;
        cout << "---"                                        << endl;

        index += 1;
    }
}


void print_annotation_record(
    const string& filename,
    const int&    index
) {
    unordered_map<string, u_int64_t> record = annotation_decoder(filename)[index];

    cout << "Days Uplink: "           << record["days_ul"]         << endl;
    cout << "Milliseconds Uplink: "   << record["milliseconds_ul"] << endl;
    cout << "Microseconds Uplink: "   << record["microseconds_ul"] << endl;
    cout << "Days Downlink: "         << record["days_dl"]         << endl;
    cout << "Milliseconds Downlink: " << record["milliseconds_dl"] << endl;
    cout << "Microseconds Downlink: " << record["microseconds_dl"] << endl;
    cout << "Packet Length: "         << record["packet_length"]   << endl;
    cout << "Error Flag: "            << record["error_flag"]      << endl;
}


void print_packet_at_index(
    string filename,
    int index,
    bool headers = true,
    bool pulse_info = false,
    bool modes = false
) {
    ifstream data = open_file(filename);

    vector<L0Packet> packets = L0Packet::get_packets(data, index + 1);

    if (index >= packets.size())
    {
        throw out_of_range("Index is greater than the number of packets.");
    }
    L0Packet packet = packets[index];

    if (modes)      packet.print_modes();
    if (pulse_info) packet.print_pulse_info();
    if (headers)
    {
        packet.print_primary_header();
        packet.print_secondary_header();
    }
}


int main(int argc, char* argv[]) 
{
    if(argv[1] == __null) 
    {
        cout << "Please enter a command." << endl;
        return 1;
    }
    string command = string(argv[1]);

    if (command == "help" or command  == "--help" or command == "-h")
    {
        cout << "print [packet_index] [path]"                    << endl;
        cout << "print_headers [packet_index] [path]"            << endl;
        cout << "print_modes [packet_index] [path]"              << endl;
        cout << "print_pulse_info [packet_index] [path]"         << endl;
        cout << "print_complex_samples [packet_index] [path]"    << endl;
        cout << "print_index_records [path]"                     << endl;
        cout << "print_annotation_record [record_index] [path]"  << endl;
        cout << "print_swaths [path]"                            << endl;
        cout << "find_packets_of_type [packet_type] [path]"      << endl;
        cout << "time [num_packets] [path]"                      << endl;
        cout << "thread_test [path]"                             << endl;
        cout << "omp_test [path]"                                << endl;
    }

    else if (command == "print")
    {
        vector<string> args = {"packet_index", "path"};
        vector<string> arg_types = {"int", "path"};
        validate_args(command, args, arg_types, &(argv[0]));
        print_packet_at_index(string(argv[3]), stoi(argv[2]), true, true, true);
    }

    else if (command == "print_headers")
    {
        vector<string> args = {"packet_index", "path"};
        vector<string> arg_types = {"int", "path"};
        validate_args(command, args, arg_types, &(argv[0]));
        print_packet_at_index(string(argv[3]), stoi(argv[2]));
    }

    else if (command == "print_pulse_info")
    {
        vector<string> args = {"packet_index", "path"};
        vector<string> arg_types = {"int", "path"};
        validate_args(command, args, arg_types, &(argv[0]));
        print_packet_at_index(string(argv[3]), stoi(argv[2]), false, true);
    }

    else if (command == "print_modes")
    {
        vector<string> args = {"packet_index", "path"};
        vector<string> arg_types = {"int", "path"};
        validate_args(command, args, arg_types, &(argv[0]));
        print_packet_at_index(string(argv[3]), stoi(argv[2]), false, false, true);
    }

    else if (command == "print_index_records")
    {
        vector<string> args = {"path"};
        vector<string> arg_types = {"path"};
        validate_args(command, args, arg_types, &(argv[0]));
        print_index_records(string(argv[2]));
    }

    else if (command == "print_annotation_record")
    {
        vector<string> args = {"record_index", "path"};
        vector<string> arg_types = {"int", "path"};
        validate_args(command, args, arg_types, &(argv[0]));
        print_annotation_record(string(argv[3]), stoi(argv[2]));
    }

    else if (command == "print_complex_samples")
    {
        vector<string> args = {"packet_index", "path"};
        vector<string> arg_types = {"int", "path"};
        validate_args(command, args, arg_types, &(argv[0]));

        ifstream data = open_file(string(argv[3]));

        vector<L0Packet>       packets = L0Packet::get_packets(data, stoi(argv[2]) + 1);
        vector<complex<float>> signal  = packets[stoi(argv[2])].get_signal();

        for (complex<float> sample : signal)
        {
            cout << sample << endl;
        }
    }

    else if (command == "find_packets_of_type")
    {
        vector<string> args = {"packet_type", "path"};
        vector<string> types = {"char", "path"};
        validate_args(command, args, types, &(argv[0]));

        char type = char(argv[2][0]);

        vector<L0Packet> packets = L0Packet::get_packets(string(argv[3]));

        for (int i = 0; i < packets.size(); i++)
        {
            L0Packet packet = packets[i];
            char data_format    = packet.get_data_format();
            int  sequence_count = packet.primary_header("packet_sequence_count");
            if (data_format == type)
            {
                cout << "Packet #" << sequence_count << " at index "
                     << i << " is type " << data_format << endl;
            }
        }
    }

    else if (command == "print_swaths")
    {
        vector<string> args = {"path"};
        vector<string> types = {"path"};
        validate_args(command, args, types, &(argv[0]));

        vector<L0Packet> packets = L0Packet::get_packets(string(argv[2]));
        
        unordered_map<string, int> swaths;

        for (L0Packet packet : packets) 
        {
            swaths[packet.get_swath()]++;
        }
        for (pair<string, int> swath : swaths)
        {
            cout << swath.first << ": " << swath.second << endl;
        }
    }

    else if (command == "time")
    {
        vector<string> args = {"packet_count", "path"};
        vector<string> types = {"int", "path"};
        validate_args(command, args, types, &(argv[0]));

        double runtime = time_packet_generation(string(argv[3]), stoi(argv[2]), false, 0);

        cout << "Decoded " << stoi(argv[2]) << " packets in " << runtime << "s." << endl;
    }

    else if (command == "thread_test")
    {
        vector<string> args = {"path"};
        vector<string> types = {"path"};
        validate_args(command, args, types, &(argv[0]));

        thread_test(string(argv[2]));
    }

    else if (command == "omp_test")
    {
        vector<string> args = {"path"};
        vector<string> types = {"path"};
        validate_args(command, args, types, &(argv[0]));

        omp_test(string(argv[2]));
    }

    else
    {
        cout << command << " is not a valid command." << endl;
    }
    return 0;
}