#pragma once

#include <set>

#include "packet.h"
#include "image_formation.h"
#include "misc_types.h"
#include "state_vectors.h"
#include "structs.h"
#include "doppler.h"


class S1_Decoder {
private:

    std::string _filename;

    bool _use_eccm = false;
    int _eccm_detection_threshold;
    int _eccm_mitigation_threshold;

    PACKET_VEC_1D _flat_packets;

    D_VEC_1D _times;

    std::unordered_map<std::string, int> _swath_counts;

    std::unordered_map<std::string, PACKET_VEC_2D> _echo_packets;
    std::unordered_map<std::string, PACKET_VEC_2D> _cal_packets;

    STATE_VECTORS _state_vectors;

    void _validate_request(const std::string& swath, const int& burst=0);

    CF_VEC_2D _get_cal_swath(const std::string& swath);

    CF_VEC_2D _range_compress(PACKET_VEC_1D& packets, const bool& do_ifft=true, const bool& do_azimuth_fft=false);
    CF_VEC_2D _azimuth_compress(PACKET_VEC_1D& packets, const bool& tops_mode=false);

    CF_VEC_2D _get_range_compressed_swath_sm(const std::string& swath, const bool& range_doppler=false);
    CF_VEC_2D _get_range_compressed_swath_iw(const std::string& swath, const bool& range_doppler=false);

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

    void apply_eccm() {_use_eccm = true;}
    void apply_eccm(const int& detection_threshold, const int& mitigation_threshold);

    CF_VEC_2D get_burst(const std::string& swath, const int& burst);
    CF_VEC_2D get_swath(const std::string& swath);

    CF_VEC_2D get_range_compressed_burst(const std::string& swath, const int& burst, const bool& range_doppler=false);
    CF_VEC_2D get_range_compressed_swath(const std::string& swath, const bool& range_doppler=false);

    CF_VEC_2D get_azimuth_compressed_burst(const std::string& swath, const int& burst);
    CF_VEC_2D get_azimuth_compressed_swath(const std::string& swath);
};
