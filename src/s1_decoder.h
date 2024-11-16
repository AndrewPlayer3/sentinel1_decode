#pragma once

#include <set>

#include "packet.h"
#include "image_formation.h"
#include "misc_types.h"
#include "state_vectors.h"
#include "structs.h"


class S1_Decoder {
private:

    std::string _filename;
    PACKET_VEC_1D _flat_packets;
    std::unordered_map<std::string, int> _swath_counts;
    std::unordered_map<std::string, PACKET_VEC_2D> _echo_packets;
    std::unordered_map<std::string, PACKET_VEC_2D> _cal_packets;
    STATE_VECTORS _state_vectors;

public:

    S1_Decoder() {}

    S1_Decoder(const std::string& filename) 
    {
        _filename = filename;
        _set_packets();
        _set_state_vectors();
    }

    void _set_state_vectors()
    {
        _state_vectors = STATE_VECTORS(_flat_packets);
    }

    void _set_packets()
    {
        _flat_packets = L0Packet::get_packets(_filename, 0);

        for (L0Packet packet : _flat_packets) 
        {
            _swath_counts[packet.get_swath()]++;
        }

        for (std::pair<std::string, int> swath_count : _swath_counts)
        {
            std::string name = swath_count.first;
            if (ECHO_SWATHS.contains(name))  
            {
                _echo_packets[name] = L0Packet::get_packets_in_bursts(_flat_packets, name);
            }
            else if (CAL_SWATHS.contains(name))
            {
                _cal_packets[name] = L0Packet::get_packets_in_bursts(_flat_packets, name, true);
            }
        }
    }

    STATE_VECTORS get_state_vectors();

    std::pair<PACKET_VEC_2D, int> get_azimuth_blocks(PACKET_VEC_1D& packets);

    CF_VEC_2D get_burst(const std::string& swath, const int& burst);
    CF_VEC_2D get_swath(const std::string& swath);

    CF_VEC_2D range_compress(PACKET_VEC_1D& packets);
    CF_VEC_2D azimuth_compress(PACKET_VEC_1D& packets);

    CF_VEC_2D get_range_compressed_burst(const std::string& swath, const int& burst);
    CF_VEC_2D get_range_compressed_swath_sm(const std::string& swath);
    CF_VEC_2D get_range_compressed_swath_iw(const std::string& swath);
    CF_VEC_2D get_range_compressed_swath(const std::string& swath);

    CF_VEC_2D get_azimuth_compressed_burst(const std::string& swath, const int& burst);
    CF_VEC_2D get_azimuth_compressed_swath_sm(const std::string& swath);
    CF_VEC_2D get_azimuth_compressed_swath_iw(const std::string& swath);
    CF_VEC_2D get_azimuth_compressed_swath(const std::string& swath);
};
