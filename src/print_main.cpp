/*
By: Andrew Player
Name: main.cpp
Description: Main function with random command-line arguments for testing purposes. 
             See "SAR Space Packet Protocol Data Unit", for more information on the packet specification:
             https://sentinels.copernicus.eu/documents/247904/2142675/Sentinel-1-SAR-Space-Packet-Protocol-Data-Unit.pdf
             For additional information on Level-0 product decoding, see:
             https://sentinel.esa.int/documents/247904/0/Sentinel-1-Level-0-Data-Decoding-Package.pdf/a8742c59-4914-40c4-8309-c77515649f17
*/

#include "cli.h"
#include "aux_decoding.h"
#include "misc_types.h"
#include "state_vectors.h"


void print_state_vectors(
    std::string filename
) {
    PACKET_VEC_1D packets = L0Packet::get_packets(filename, 0);
    STATE_VECTORS state_vectors(packets);
    state_vectors.print();
}


void print_packet_at_index(
    std::string filename,
    int  index
) {
    std::ifstream data = open_file(filename);
    PACKET_VEC_1D packets = L0Packet::get_packets(data, index + 1);

    if (index >= packets.size())
    {
        throw std::out_of_range("Index is greater than the number of packets.");
    }
    L0Packet packet = packets[index];

    std::cout << "Primary Header:" << std::endl;
    packet.print_primary_header();
    std::cout << "\nSecondary Header: " << std::endl;
    packet.print_secondary_header();
    std::cout << "\nOperating Mode Info:" << std::endl;
    packet.print_modes();
    std::cout << "\nPulse Info:" << std::endl;
    packet.print_pulse_info();
}


int main(int argc, char* argv[]) 
{
    STRING_VEC_1D help_strings = {
        "packet_info [packet_index] [path]",
        "complex_samples [packet_index] [path]",
        "swath_names [path]",
        "index_records [path]",
        "annotation_record [record_index] [path]",
        "state_vectors [path]"
    };

    if(argv[1] == __null) 
    {
        std::cout << "Please enter a command:" << std::endl;
        print_help(help_strings);
        return 1;
    }
    std::string command = std::string(argv[1]);

    if (command == "packet_info")
    {
        STRING_VEC_1D args  = {"packet_index", "path"};
        STRING_VEC_1D types = {"int", "path"};
        validate_args(command, args, types, &(argv[0]));
        print_packet_at_index(std::string(argv[3]), std::stoi(argv[2]));
    }

    else if (command == "swath_names")
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

    else if (command == "complex_samples")
    {
        STRING_VEC_1D args  = {"packet_index", "path"};
        STRING_VEC_1D types = {"int", "path"};
        validate_args(command, args, types, &(argv[0]));

        std::ifstream data = open_file(std::string(argv[3]));

        PACKET_VEC_1D       packets = L0Packet::get_packets(data, std::stoi(argv[2]) + 1);
        CF_VEC_1D signal  = packets[std::stoi(argv[2])].get_signal();

        for (std::complex<double> sample : signal)
        {
            std::cout << sample << std::endl;
        }
    }

    else if (command == "index_records")
    {
        STRING_VEC_1D args  = {"path"};
        STRING_VEC_1D types = {"path"};
        validate_args(command, args, types, &(argv[0]));
        print_index_records(std::string(argv[2]));
    }

    else if (command == "annotation_record")
    {
        STRING_VEC_1D args  = {"record_index", "path"};
        STRING_VEC_1D types = {"int", "path"};
        validate_args(command, args, types, &(argv[0]));
        print_annotation_record(std::string(argv[3]), std::stoi(argv[2]));
    }

    else if (command == "state_vectors")
    {
        STRING_VEC_1D args  = {"path"};
        STRING_VEC_1D types = {"path"};
        validate_args(command, args, types, &(argv[0]));
        print_state_vectors(std::string(argv[2]));
    }

    else if (command == "help" or command  == "--help" or command == "-h")
    {
        print_help(help_strings);
    }

    else
    {
        std::cout << command << " is not a valid command." << std::endl;
    }
    return 0;
}
