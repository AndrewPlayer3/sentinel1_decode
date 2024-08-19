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