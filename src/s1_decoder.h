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

    CF_VEC_2D _range_compress(PACKET_VEC_1D& packets);
    CF_VEC_2D _azimuth_compress(PACKET_VEC_1D& packets);

    CF_VEC_2D _get_range_compressed_swath_sm(const std::string& swath);
    CF_VEC_2D _get_range_compressed_swath_iw(const std::string& swath);

    CF_VEC_2D _get_azimuth_compressed_swath_sm(const std::string& swath);
    CF_VEC_2D _get_azimuth_compressed_swath_iw(const std::string& swath);

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

    void _set_packets();

    STATE_VECTORS get_state_vectors();

    std::pair<PACKET_VEC_2D, int> get_azimuth_blocks(PACKET_VEC_1D& packets);

    CF_VEC_2D get_burst(const std::string& swath, const int& burst);
    CF_VEC_2D get_swath(const std::string& swath);

    CF_VEC_2D get_range_compressed_burst(const std::string& swath, const int& burst);
    CF_VEC_2D get_range_compressed_swath(const std::string& swath);

    CF_VEC_2D get_azimuth_compressed_burst(const std::string& swath, const int& burst);
    CF_VEC_2D get_azimuth_compressed_swath(const std::string& swath);
};
