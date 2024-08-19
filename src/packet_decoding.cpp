/*
By: Andrew Player
Name: packet_decoding.cpp
Description: Funtions for reading packets, index info, and annotation info from the Level-0 dat files. 
             See "SAR Space Packet Protocol Data Unit", for more information on the packet specification:
             https://sentinels.copernicus.eu/documents/247904/2142675/Sentinel-1-SAR-Space-Packet-Protocol-Data-Unit.pdf
             For additional information on Level-0 product decoding, see:
             https://sentinel.esa.int/documents/247904/0/Sentinel-1-Level-0-Data-Decoding-Package.pdf/a8742c59-4914-40c4-8309-c77515649f17
*/

#include "packet_decoding.hpp"

using namespace std;


/* Returns the header dictionary with the values cast to integers with no calculations */
unordered_map<string, int> get_header_dict(
    const vector<u_int8_t>&  bytes,
    const vector<int>&       bit_lengths,
    const vector<string>&    field_names
) {
    int num_fields = field_names.size();
    int bit_index = 0;

    unordered_map<string, int> header;

    for (int i = 0; i < num_fields; i++) 
    {
        header[field_names[i]] = read_n_bits(bytes, bit_index, bit_lengths[i]);

        bit_index += bit_lengths[i];
    }
    return header;
}


/* Decode the next packet within the data stream */
L0Packet decode_next_packet(ifstream& data)
{
    vector<u_int8_t> primary_bytes = read_bytes(data, 6);
    unordered_map<string, int> primary_header = get_header_dict(
        primary_bytes,
        PRIMARY_HEADER,
        PRIMARY_HEADER_FIELDS
    );

    if (data.eof()) return L0Packet();

    vector<u_int8_t> secondary_bytes = read_bytes(data, 62);
    unordered_map<string, int> secondary_header = get_header_dict(
        secondary_bytes,
        SECONDARY_HEADER,
        SECONDARY_HEADER_FIELDS
    );

    u_int32_t packet_length    = primary_header["packet_data_length"];
    u_int32_t user_data_length = packet_length + 1 - SECONDARY_HEADER_SIZE;

    vector<u_int8_t> user_data = read_bytes(data, user_data_length);

    L0Packet packet = L0Packet(
        primary_header,
        secondary_header,
        user_data
    );

    return packet;
}


/* Returns a vector of maps containing the index records */
vector<unordered_map<string, int>> index_decoder(ifstream& data)
{
    int record_size = 36;

    vector<unordered_map<string, int>> index_records;

    while (!data.eof())
    {
        vector<u_int8_t> record_bytes = read_bytes(data, record_size);

        int date_time     = read_n_bits(record_bytes,   0, 64);
        int time_delta    = read_n_bits(record_bytes,  64, 64);
        int data_size     = read_n_bits(record_bytes, 128, 32);
        int unit_offset   = read_n_bits(record_bytes, 160, 32);
        int byte_offset   = read_n_bits(record_bytes, 192, 64);
        int variable_flag = read_n_bits(record_bytes, 256,  8);
        int spare_data    = read_n_bits(record_bytes, 264, 24);

        index_records.push_back(unordered_map<string, int>({
            {"date_time",     date_time},
            {"time_delta",    time_delta},
            {"data_size",     data_size},
            {"unit_offset",   unit_offset},
            {"byte_offset",   byte_offset},
            {"variable_flag", variable_flag},
            {"spare_data",    spare_data}
        }));
    }
    return index_records;
}


/* Returns a vector of maps containing the annotation records */
vector<unordered_map<string, int>> annotation_decoder(ifstream& data)
{
    int record_size = 26;

    vector<unordered_map<string, int>> annotation_records;

    while (!data.eof())
    {
        vector<u_int8_t> record_bytes = read_bytes(data, record_size);

        int days_ul             = read_n_bits(record_bytes,   0, 16);
        int milliseconds_ul     = read_n_bits(record_bytes,  16, 32);
        int microseconds_ul     = read_n_bits(record_bytes,  48, 16);
        int days_dl             = read_n_bits(record_bytes,  64, 16);
        int milliseconds_dl     = read_n_bits(record_bytes,  80, 32);
        int microseconds_dl     = read_n_bits(record_bytes, 112, 16);
        int packet_length       = read_n_bits(record_bytes, 128, 16);
        int num_transfer_frames = read_n_bits(record_bytes, 144, 16);
        int error_flag          = read_n_bits(record_bytes, 160,  8);
        int bit_1_type          = read_n_bits(record_bytes, 168,  1);
        int bit_2_type          = read_n_bits(record_bytes, 169,  2);
        int bit_6_type          = read_n_bits(record_bytes, 171,  6);
        int spare_field         = read_n_bits(record_bytes, 177,  8);

        annotation_records.push_back(unordered_map<string, int>({
            {"days_ul",             days_ul},
            {"milliseconds_ul",     milliseconds_ul},
            {"microseconds_ul",     microseconds_ul},
            {"days_dl",             days_dl},
            {"milliseconds_dl",     milliseconds_dl},
            {"microseconds_dl",     microseconds_dl},
            {"packet_length",       packet_length},
            {"num_transfer_frames", num_transfer_frames},
            {"error_flag",          error_flag},
            {"bit_1_type",          bit_1_type},
            {"bit_2_type",          bit_2_type},
            {"bit_6_type",          bit_6_type},
            {"spare_field",         spare_field}
        }));
    }
    return annotation_records;
}


/* Returns the time in seconds that it takes to decode the complex samples of num_packets packets */
double time_packet_generation(ifstream& data, const int& num_packets, const bool& log, const int& log_interval)
{
    double total_runtime = 0.0;
    for (int i = 0; i < num_packets; i++)
    {
        try
        {
            auto start = chrono::high_resolution_clock::now();

            L0Packet packet = decode_next_packet(data);
            
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


/* Returns a vector containing all of the packets left in a data stream */
vector<L0Packet> get_all_packets(ifstream& data, const bool& log, const int& log_interval)
{
    vector<L0Packet> packets;

    int index = 0;

    while (!data.eof())
    {
        try
        {
            L0Packet packet = decode_next_packet(data);

            if (!packet.is_empty())
            {
                packets.push_back(packet);
            }

            if (log && index % log_interval == 0 && index != 0)
            {
                cout << "Read " << index << " packets." << endl;
            }

            index += 1;
        }
        catch(runtime_error)
        {
            cout << "Runtime error encountered while decoding packet #" << index << endl;
            continue;
        }
        catch(length_error)
        {
            break;
        }
    }
    if (log) cout << "EOF reached after decoding " << index << " packets." << endl;
    return packets;
}


/* Returns a vector containing the n first packets left in a data stream */
vector<L0Packet> get_n_packets(ifstream& data, const int& n, const bool& log, const int& log_interval)
{
    vector<L0Packet> packets;

    int index = 0;

    while (index < n)
    {
        try
        {
            L0Packet packet = decode_next_packet(data);
            
            packets.push_back(packet);

            if (log && index % log_interval == 0 && index != 0)
            {
                cout << "Read " << index << " packets." << endl;
            }

            index += 1;
        }
        catch(runtime_error)
        {
            cout << "Caught a runtime error while decoding packet #" << index << endl;
            continue;
        }
        catch(length_error)
        {
            throw out_of_range("n is greater than the available number of packets.");
        }
    }
    return packets;
}