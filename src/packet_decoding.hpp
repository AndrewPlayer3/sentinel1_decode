/*
By: Andrew Player
Name: packet_decoding.hpp
Description: Funtions for reading packets, index info, and annotation info from the Level-0 dat files. 
             See "SAR Space Packet Protocol Data Unit", for more information on the packet specification:
             https://sentinels.copernicus.eu/documents/247904/2142675/Sentinel-1-SAR-Space-Packet-Protocol-Data-Unit.pdf
             For additional information on Level-0 product decoding, see:
             https://sentinel.esa.int/documents/247904/0/Sentinel-1-Level-0-Data-Decoding-Package.pdf/a8742c59-4914-40c4-8309-c77515649f17
*/

#pragma once

#include "packet.hpp"

unordered_map<string, int> get_header_dict(
    const vector<u_int8_t>& bytes,
    const vector<int>&      bit_lengths,
    const vector<string>&   field_names
);

L0Packet decode_next_packet(ifstream& data);

vector<L0Packet> get_all_packets(ifstream& data, const bool& log, const int& log_interval);

vector<L0Packet> get_n_packets(ifstream& data, const int& n, const bool& log, const int& log_interval);

vector<unordered_map<string, int>> annotation_decoder(ifstream& data);

vector<unordered_map<string, int>> index_decoder(ifstream& data);

double time_packet_generation(ifstream& data, const int& num_packets, const bool& log, const int& log_interval);