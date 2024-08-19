/*
By: Andrew Player
Name: packet_decoding.cpp
Description: Funtions for reading packets, index info, and annotation info from the Level-0 dat files. 
             See "SAR Space Packet Protocol Data Unit", for more information on the packet specification:
             https://sentinels.copernicus.eu/documents/247904/2142675/Sentinel-1-SAR-Space-Packet-Protocol-Data-Unit.pdf
             For additional information on Level-0 product decoding, see:
             https://sentinel.esa.int/documents/247904/0/Sentinel-1-Level-0-Data-Decoding-Package.pdf/a8742c59-4914-40c4-8309-c77515649f17
*/

#include "aux_decoding.hpp"

using namespace std;


vector<unordered_map<string, u_int64_t>> index_decoder(const string& filename)
{
    ifstream data = open_file(filename);
    return index_decoder(data);
}


/* Returns a vector of maps containing the index records */
vector<unordered_map<string, u_int64_t>> index_decoder(ifstream& data)
{
    int record_size = 36;

    vector<unordered_map<string, u_int64_t>> index_records;

    while (!data.eof())
    {
        vector<u_int8_t> record_bytes = read_bytes(data, record_size);

        index_records.push_back(unordered_map<string, u_int64_t>({
            {"date_time",     read_n_bits(record_bytes,   0, 64)},
            {"time_delta",    read_n_bits(record_bytes,  64, 64)},
            {"data_size",     read_n_bits(record_bytes, 128, 32)},
            {"unit_offset",   read_n_bits(record_bytes, 160, 32)},
            {"byte_offset",   read_n_bits(record_bytes, 192, 64)},
            {"variable_flag", read_n_bits(record_bytes, 256,  8)},
            {"spare_data",    read_n_bits(record_bytes, 264, 24)}
        }));
    }
    return index_records;
}


vector<unordered_map<string, u_int64_t>> annotation_decoder(const string& filename)
{
    ifstream data = open_file(filename);
    return annotation_decoder(data);
}


/* Returns a vector of maps containing the annotation records */
vector<unordered_map<string, u_int64_t>> annotation_decoder(ifstream& data)
{
    int record_size = 26;

    vector<unordered_map<string, u_int64_t>> annotation_records;

    while (!data.eof())
    {
        vector<u_int8_t> record_bytes = read_bytes(data, record_size);

        annotation_records.push_back(unordered_map<string, u_int64_t>({
            {"days_ul",             read_n_bits(record_bytes,   0, 16)},
            {"milliseconds_ul",     read_n_bits(record_bytes,  16, 32)},
            {"microseconds_ul",     read_n_bits(record_bytes,  48, 16)},
            {"days_dl",             read_n_bits(record_bytes,  64, 16)},
            {"milliseconds_dl",     read_n_bits(record_bytes,  80, 32)},
            {"microseconds_dl",     read_n_bits(record_bytes, 112, 16)},
            {"packet_length",       read_n_bits(record_bytes, 128, 16)},
            {"num_transfer_frames", read_n_bits(record_bytes, 144, 16)},
            {"error_flag",          read_n_bits(record_bytes, 160,  8)},
            {"bit_1_type",          read_n_bits(record_bytes, 168,  1)},
            {"bit_2_type",          read_n_bits(record_bytes, 169,  2)},
            {"bit_6_type",          read_n_bits(record_bytes, 171,  6)},
            {"spare_field",         read_n_bits(record_bytes, 177,  8)}
        }));
    }
    return annotation_records;
}