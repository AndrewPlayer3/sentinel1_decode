/*
By: Andrew Player
Name: decoding_utils.hpp
Description: Functions to assist reading and decompressing binary data. 
             See "SAR Space Packet Protocol Data Unit", for more information on the packet specification:
             https://sentinels.copernicus.eu/documents/247904/2142675/Sentinel-1-SAR-Space-Packet-Protocol-Data-Unit.pdf
             For additional information on Level-0 product decoding, see:
             https://sentinel.esa.int/documents/247904/0/Sentinel-1-Level-0-Data-Decoding-Package.pdf/a8742c59-4914-40c4-8309-c77515649f17
*/

#pragma once

#include <chrono>
#include <fstream>
#include <unordered_map>
#include <complex>

#include "structs.hpp"


u_int16_t huffman_decode(
    const vector<u_int8_t>& data,
    const int& brc,
          int& bit_index
);

u_int16_t huffman_decode_with_length(
    const vector<u_int8_t>& data,
    const int& brc,
          int& bit_index
);

vector<u_int8_t> read_bytes(
    ifstream&  data,
    const int& num_bytes
);

u_int64_t read_n_bits(
    const std::vector<u_int8_t>& data,
    const int& start_bit,
    const int& n
);

ifstream open_file(const string& filename);