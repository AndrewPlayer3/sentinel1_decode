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

#include "perf.h"
#include "aux_decoding.h"
#include "misc_types.h"


void validate_args(
    const std::string&   command,
    const STRING_VEC_1D& command_args,
    const STRING_VEC_1D& arg_types,
          char*          args[],
    const int&           start_index = 2,
    const std::string&   help = ""
) {
    int end_index = command_args.size() + start_index;
    for (int i = start_index; i < end_index; i++)
    {
        int index = i - start_index;
        if (args[i] == __null)
        {
            std::cout << command  << " is missing the following argument: "
                 << command_args[index] << std::endl;
            exit(1);
        }
        try
        {
            if      (arg_types[index] == "string") std::string(args[i]);
            else if (arg_types[index] == "char")   char(args[i][0]);
            else if (arg_types[index] == "int")    std::stoi(args[i]);
            else if (arg_types[index] == "float")  std::stof(args[i]);
            else if (arg_types[index] == "path")   open_file(std::string(args[i]));
        }
        catch(...)
        {
            std::cout << args[i] << " is not a valid " << command_args[index] << "." << std::endl;
            exit(1);
        }
    }
}


void print_index_records(const std::string& filename)
{
    std::vector<std::unordered_map<std::string, u_int64_t>> records = index_decoder(filename);

    int index = 0;
    for (std::unordered_map<std::string, u_int64_t> record : records)
    {
        std::cout << "Index Record #"  << index << ":"            << std::endl;
        std::cout << "Date/Time: "     << record["date_time"]     << std::endl;
        std::cout << "Time Delta: "    << record["time_delta"]    << std::endl;
        std::cout << "Data Size: "     << record["data_size"]     << std::endl;
        std::cout << "Unit Offset: "   << record["unit_offset"]   << std::endl;
        std::cout << "Byte Offset: "   << record["byte_offset"]   << std::endl;
        std::cout << "Variable Flag: " << record["variable_flag"] << std::endl;
        std::cout << "Spare Data: "    << record["spare_data"]    << std::endl;
        std::cout << "---"                                        << std::endl;

        index += 1;
    }
}


void print_annotation_record(
    const std::string& filename,
    const int&         index
) {
    std::unordered_map<std::string, u_int64_t> record = annotation_decoder(filename)[index];

    std::cout << "Days Uplink: "           << record["days_ul"]         << std::endl;
    std::cout << "Milliseconds Uplink: "   << record["milliseconds_ul"] << std::endl;
    std::cout << "Microseconds Uplink: "   << record["microseconds_ul"] << std::endl;
    std::cout << "Days Downlink: "         << record["days_dl"]         << std::endl;
    std::cout << "Milliseconds Downlink: " << record["milliseconds_dl"] << std::endl;
    std::cout << "Microseconds Downlink: " << record["microseconds_dl"] << std::endl;
    std::cout << "Packet Length: "         << record["packet_length"]   << std::endl;
    std::cout << "Error Flag: "            << record["error_flag"]      << std::endl;
}


void print_packet_at_index(
    std::string filename,
    int  index,
    bool headers    = true,
    bool pulse_info = false,
    bool modes      = false
) {
    std::ifstream data = open_file(filename);

    PACKET_VEC_1D packets = L0Packet::get_packets(data, index + 1);

    if (index >= packets.size())
    {
        throw std::out_of_range("Index is greater than the number of packets.");
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
        std::cout << "Please enter a command." << std::endl;
        return 1;
    }
    std::string command = std::string(argv[1]);

    if (command == "help" or command  == "--help" or command == "-h")
    {
        std::cout << "print [packet_index] [path]"                    << std::endl;
        std::cout << "print_headers [packet_index] [path]"            << std::endl;
        std::cout << "print_modes [packet_index] [path]"              << std::endl;
        std::cout << "print_pulse_info [packet_index] [path]"         << std::endl;
        std::cout << "print_complex_samples [packet_index] [path]"    << std::endl;
        std::cout << "print_index_records [path]"                     << std::endl;
        std::cout << "print_annotation_record [record_index] [path]"  << std::endl;
        std::cout << "print_swaths [path]"                            << std::endl;
        std::cout << "find_packets_of_type [packet_type] [path]"      << std::endl;
        std::cout << "time [num_packets] [path]"                      << std::endl;
        std::cout << "thread_test [path]"                             << std::endl;
        std::cout << "omp_test [path]"                                << std::endl;
    }

    else if (command == "print")
    {
        STRING_VEC_1D args  = {"packet_index", "path"};
        STRING_VEC_1D types = {"int", "path"};
        validate_args(command, args, types, &(argv[0]));
        print_packet_at_index(std::string(argv[3]), std::stoi(argv[2]), true, true, true);
    }

    else if (command == "print_headers")
    {
        STRING_VEC_1D args  = {"packet_index", "path"};
        STRING_VEC_1D types = {"int", "path"};
        validate_args(command, args, types, &(argv[0]));
        print_packet_at_index(std::string(argv[3]), std::stoi(argv[2]));
    }

    else if (command == "print_pulse_info")
    {
        STRING_VEC_1D args  = {"packet_index", "path"};
        STRING_VEC_1D types = {"int", "path"};
        validate_args(command, args, types, &(argv[0]));
        print_packet_at_index(std::string(argv[3]), std::stoi(argv[2]), false, true);
    }

    else if (command == "print_modes")
    {
        STRING_VEC_1D args  = {"packet_index", "path"};
        STRING_VEC_1D types = {"int", "path"};
        validate_args(command, args, types, &(argv[0]));
        print_packet_at_index(std::string(argv[3]), std::stoi(argv[2]), false, false, true);
    }

    else if (command == "print_index_records")
    {
        STRING_VEC_1D args  = {"path"};
        STRING_VEC_1D types = {"path"};
        validate_args(command, args, types, &(argv[0]));
        print_index_records(std::string(argv[2]));
    }

    else if (command == "print_annotation_record")
    {
        STRING_VEC_1D args  = {"record_index", "path"};
        STRING_VEC_1D types = {"int", "path"};
        validate_args(command, args, types, &(argv[0]));
        print_annotation_record(std::string(argv[3]), std::stoi(argv[2]));
    }

    else if (command == "print_complex_samples")
    {
        STRING_VEC_1D args  = {"packet_index", "path"};
        STRING_VEC_1D types = {"int", "path"};
        validate_args(command, args, types, &(argv[0]));

        std::ifstream data = open_file(std::string(argv[3]));

        PACKET_VEC_1D       packets = L0Packet::get_packets(data, std::stoi(argv[2]) + 1);
        CF_VEC_1D signal  = packets[std::stoi(argv[2])].get_signal();

        for (std::complex<float> sample : signal)
        {
            std::cout << sample << std::endl;
        }
    }

    else if (command == "find_packets_of_type")
    {
        STRING_VEC_1D args = {"packet_type", "path"};
        STRING_VEC_1D types = {"char", "path"};
        validate_args(command, args, types, &(argv[0]));

        char type = char(argv[2][0]);

        PACKET_VEC_1D packets = L0Packet::get_packets(std::string(argv[3]));

        for (int i = 0; i < packets.size(); i++)
        {
            L0Packet packet     = packets[i];
            char data_format    = packet.get_data_format();
            int  sequence_count = packet.primary_header("packet_sequence_count");
            if (data_format == type)
            {
                std::cout << "Packet #" << sequence_count << " at index "
                     << i << " is type " << data_format << std::endl;
            }
        }
    }

    else if (command == "print_swaths")
    {
        STRING_VEC_1D args  = {"path"};
        STRING_VEC_1D types = {"path"};
        validate_args(command, args, types, &(argv[0]));

        PACKET_VEC_1D packets = L0Packet::get_packets(std::string(argv[2]));
        
        std::unordered_map<std::string, int> swaths;

        for (L0Packet packet : packets) 
        {
            swaths[packet.get_swath()]++;
        }
        for (std::pair<std::string, int> swath : swaths)
        {
            std::cout << swath.first << ": " << swath.second << std::endl;
        }
    }

    else if (command == "time")
    {
        STRING_VEC_1D args  = {"packet_count", "path"};
        STRING_VEC_1D types = {"int", "path"};
        validate_args(command, args, types, &(argv[0]));

        double runtime = time_packet_generation(std::string(argv[3]), std::stoi(argv[2]), false, 0);

        std::cout << "Decoded " << std::stoi(argv[2]) << " packets in " << runtime << "s." << std::endl;
    }

    else if (command == "thread_test")
    {
        STRING_VEC_1D args  = {"path"};
        STRING_VEC_1D types = {"path"};
        validate_args(command, args, types, &(argv[0]));

        thread_test(std::string(argv[2]));
    }

    else if (command == "omp_test")
    {
        STRING_VEC_1D args  = {"path"};
        STRING_VEC_1D types = {"path"};
        validate_args(command, args, types, &(argv[0]));

        omp_test(std::string(argv[2]));
    }

    else
    {
        std::cout << command << " is not a valid command." << std::endl;
    }
    return 0;
}