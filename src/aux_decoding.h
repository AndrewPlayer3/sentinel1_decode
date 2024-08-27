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

#include "packet.h"
#include "misc_types.h"

std::vector<std::unordered_map<std::string, u_int64_t>> annotation_decoder(const std::string& filename);
std::vector<std::unordered_map<std::string, u_int64_t>> annotation_decoder(std::ifstream& data);
std::vector<std::unordered_map<std::string, u_int64_t>> index_decoder(const std::string& filename);
std::vector<std::unordered_map<std::string, u_int64_t>> index_decoder(std::ifstream& data);