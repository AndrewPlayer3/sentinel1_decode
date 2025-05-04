#pragma once

#include "packet.h"
#include "image_formation.h"
#include "misc_types.h"

/* Class representing a burst as a vector of packets, complex_samples, and chirps. */
class Burst {
private:
    int _burst_num;
    int _num_packets;

    double _pri;
    double _pulse_len;
    double _start_freq;
    double _swl;
    double _swst;

    double _start_time;
    double _stop_time;

    PACKET_VEC_1D _packets;

    CF_VEC_1D _reference_function;

    F_VEC_1D _slant_ranges;
    F_VEC_1D _slant_range_times;

    std::string _swath;

public:
    Burst(
        PACKET_VEC_1D& packets,
        const std::string& swath,
        const int& burst_index
    ) {
        _packets = packets;
        _burst_num = burst_index;
        _swath = swath;

        _packets = L0Packet::get_packets_in_bursts(packets, swath, false)[burst_index];
        _num_packets = _packets.size();

        L0Packet packet = _packets[0];
        _pri = packet.get_pri() * 1e-6;
        _swst = packet.get_swst() * 1e-6;
        _swl = packet.get_swl() * 1e-6;
        _pulse_len = packet.get_pulse_length() * 1e-6;
        _start_freq = packet.get_start_frequency();

        _reference_function = get_reference_function(packet.get_replica_chirp());

        _slant_ranges = packet.get_slant_ranges();
        _slant_range_times = packet.get_slant_range_times();

        _start_time = packet.get_time();
        _stop_time = _packets.back().get_time();

    }
};


typedef struct std::vector<Burst>              BURST_VEC_1D;
typedef struct std::vector<std::vector<Burst>> BURST_VEC_2D;
