/*
By: Andrew Player
Name: perf.hpp
Description: Functions for evaluating the performance of packet decoding. 
             See "SAR Space Packet Protocol Data Unit", for more information on the packet specification:
             https://sentinels.copernicus.eu/documents/247904/2142675/Sentinel-1-SAR-Space-Packet-Protocol-Data-Unit.pdf
             For additional information on Level-0 product decoding, see:
             https://sentinel.esa.int/documents/247904/0/Sentinel-1-Level-0-Data-Decoding-Package.pdf/a8742c59-4914-40c4-8309-c77515649f17
*/

#pragma once

#include <thread>

#include "packet.h"
#include "misc_types.h"
#include "image_formation.h"
#include "s1_decoder.h"

double time_packet_generation(
    const std::string& filename,
    const int&  num_packets,
    const bool& log,
    const int&  log_interval
);

void omp_test(
    const std::string& filename
);

void thread_test(
    const std::string& filename
);

/* Returns the time in seconds that it takes to range compress the given burst */
double time_range_compression(
    const std::string& filename, 
    const std::string& swath_name,
    const int&  burst_num, 
    const bool& log, 
    const int&  log_interval
);
