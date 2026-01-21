/*
By: Andrew Player
Name: packet.hpp
Description: L0Packet class for storing and decoding Level-0 Packets in a convinient and easy to use manner. 
             See "SAR Space Packet Protocol Data Unit", for more information on the packet specification:
             https://sentinels.copernicus.eu/documents/247904/2142675/Sentinel-1-SAR-Space-Packet-Protocol-Data-Unit.pdf
             For additional information on Level-0 product decoding, see:
             https://sentinel.esa.int/documents/247904/0/Sentinel-1-Level-0-Data-Decoding-Package.pdf/a8742c59-4914-40c4-8309-c77515649f17
*/

#pragma once

#include <numeric>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <complex>
#include <iomanip>

#include "decoding_utils.h"
#include "misc_types.h"
#include "signal_processing.h"


/* H_CODE / S_CODE struct representing one element of a quadrature */
struct H_CODE {
    UINT8_VEC_1D  signs;
    UINT16_VEC_1D m_codes;

    H_CODE() {
        signs.resize(128);
        m_codes.resize(128);
    }

    H_CODE(const int& num_codes)
    {
        signs.resize(num_codes);
        m_codes.resize(num_codes);
    }
};


/* Quadrature struct representing one of the IE, IO, QE, or QO quads */
struct QUAD {
    std::vector<H_CODE> blocks;

    std::string key;

    QUAD(const std::string& component_key) {
        key = component_key;
    }

    QUAD(const std::string& component_key, const int& num_blocks)
    {
        key = component_key;
        blocks.resize(num_blocks);
    }
};


/* Class for Sentinel-1 Level-0 Packets with decoding functionality for all packet types. */
class L0Packet 
{
private:
    std::unordered_map<std::string, int> _primary_header;
    std::unordered_map<std::string, int> _secondary_header;

    UINT8_VEC_1D _raw_user_data;

    int _packet_index;
    int _num_quads;
    int _test_mode;
    int _baq_mode;
    int _num_baq_blocks;
    int _user_data_length;

    bool _is_empty = true;

    char _data_format;

    void _set_data_format();

    /**********************************/
    /* DECODING COMPLEX SAMPLES       */
    /* PAGES 61 -> 85                 */
    /**********************************/

    INT_VEC_1D _brc = {};
    INT_VEC_1D _thresholds = {};
    CF_VEC_1D _signal;

    bool _signal_set_flag = false;

    void _set_signal();
    void _decode();

    [[nodiscard]] int  _get_next_word_boundary(const int& bit_index);

    [[nodiscard]] H_CODE _get_h_code_type_c(
              int&  bit_index, 
        const bool& is_last_block
    );
    [[nodiscard]] H_CODE _get_h_code_type_d(
        const u_int8_t& brc,
              int&      bit_index,
        const bool&     is_last_block
    );
    [[nodiscard]] double _get_s_values_type_c(
        const u_int16_t& threshold_id,
        const int&       sign,
        const int&       m_code
    );
    [[nodiscard]] double _get_s_values_type_d(
        const u_int8_t&  brc,
        const u_int16_t& threshold_id,
        const int&       sign,
        const int&       m_code
    );

    void _set_quad_type_d(QUAD& component, int& bit_index);
    void _set_quad_type_c(QUAD& component, int& bit_index);
    void _set_quad_types_a_and_b(H_CODE& component, int& bit_index);

    [[nodiscard]] CF_VEC_1D _get_signal_type_d(
        QUAD& IE, QUAD& IO, 
        QUAD& QE, QUAD& QO
    );
    [[nodiscard]] CF_VEC_1D _get_signal_type_c(
        QUAD& IE, QUAD& IO,
        QUAD& QE, QUAD& QO
    );
    [[nodiscard]] CF_VEC_1D _get_signal_types_a_and_b(
        H_CODE& IE, H_CODE& IO,
        H_CODE& QE, H_CODE& QO
    );


public:
    L0Packet(){}    

    L0Packet(
        std::unordered_map<std::string, int> primary_header,
        std::unordered_map<std::string, int> secondary_header,
        UINT8_VEC_1D raw_user_data,
        int packet_index
    ) {
        _primary_header   = primary_header;
        _secondary_header = secondary_header;
        _raw_user_data    = raw_user_data;
        _packet_index     = packet_index;
        _num_quads        = secondary_header["num_quadratures"];
        _test_mode        = secondary_header["test_mode"];
        _baq_mode         = secondary_header["baq_mode"];
        _user_data_length = primary_header["packet_data_length"] + 1 - SECONDARY_HEADER_SIZE;
        _num_baq_blocks   = ceil((2.0 * double(_num_quads)) / 256.0);

        if (_user_data_length != int(_raw_user_data.size())) 
        {
            std::cout << _user_data_length << " != " << _raw_user_data.size() << std::endl;
            throw std::runtime_error("The lenght of the user data field is invalid.");
        }

        _set_data_format();

        _is_empty = false;
    }

    [[nodiscard]] int primary_header(const std::string& key) const {return _primary_header.at(key);}
    [[nodiscard]] int secondary_header(const std::string& key) const {return _secondary_header.at(key);}

    [[nodiscard]] bool is_empty() const {return _is_empty;}
    [[nodiscard]] int  get_num_quads() const {return _num_quads;}
    [[nodiscard]] int  get_num_baq_blocks() const {return _num_baq_blocks;}
    [[nodiscard]] int  get_user_data_length() const {return _user_data_length;}
    [[nodiscard]] char get_data_format() const {return _data_format;}

    [[nodiscard]] int get_packet_index() const;
    [[nodiscard]] int get_baq_block_length() const;
    [[nodiscard]] int get_replica_chirp_length() const;
    [[nodiscard]] double get_time() const;
    [[nodiscard]] double get_pulse_length() const;
    [[nodiscard]] double get_tx_ramp_rate() const;
    [[nodiscard]] double get_start_frequency() const;
    [[nodiscard]] double get_pri() const;
    [[nodiscard]] double get_swl() const;
    [[nodiscard]] double get_swst() const;
    [[nodiscard]] double get_rx_gain() const;
    [[nodiscard]] double get_azimuth_beam_angle() const; 
    [[nodiscard]] double get_range_sample_rate() const;
    [[nodiscard]] char get_rx_polarization() const;
    [[nodiscard]] char get_tx_polarization() const;

    [[nodiscard]] std::string get_baq_mode() const;
    [[nodiscard]] std::string get_test_mode() const;
    [[nodiscard]] std::string get_sensor_mode() const;
    [[nodiscard]] std::string get_signal_type() const;
    [[nodiscard]] std::string get_error_status() const;
    [[nodiscard]] std::string get_swath() const;

    [[nodiscard]] D_VEC_1D get_slant_ranges(int num_ranges=0) const;
    [[nodiscard]] D_VEC_1D get_slant_range_times(int num_times=0) const;
    [[nodiscard]] D_VEC_1D get_timing_corrections() const;
    [[nodiscard]] CF_VEC_1D get_replica_chirp() const;
    [[nodiscard]] CF_VEC_1D get_signal();


    void print_primary_header();
    void print_secondary_header();
    void print_modes();
    void print_pulse_info();

    typedef struct std::vector<L0Packet>              PACKET_VEC_1D;
    typedef struct std::vector<std::vector<L0Packet>> PACKET_VEC_2D;
};

typedef L0Packet::PACKET_VEC_1D PACKET_VEC_1D;
typedef L0Packet::PACKET_VEC_2D PACKET_VEC_2D;

std::unordered_map<std::string, int> parse_packet_header(
    const UINT8_VEC_1D&  bytes,
    const INT_VEC_1D&    bit_lengths,
    const STRING_VEC_1D& field_names
);

L0Packet get_next_packet(std::ifstream& data, int packet_index);

void decode_packets_in_place(PACKET_VEC_1D& packets);

PACKET_VEC_1D get_packets(std::ifstream& data, int num_packets = 0);
PACKET_VEC_1D get_packets(const std::string& filename, int num_packets = 0);
PACKET_VEC_1D get_packets_in_swath(const std::string& filename, const std::string& swath);
PACKET_VEC_1D get_packets_in_swath(std::ifstream& data, const std::string& swath);

PACKET_VEC_2D get_packets_in_bursts(const std::string& filename, const std::string& swath);
PACKET_VEC_2D get_packets_in_bursts(std::ifstream& data, const std::string& swath, bool get_cal_packets = false);
PACKET_VEC_2D get_packets_in_bursts(PACKET_VEC_1D& packets, const std::string& swath, bool get_cal_packets = false);
PACKET_VEC_1D decode_packets(const PACKET_VEC_1D& packets);
