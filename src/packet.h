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

#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <complex>

#include "decoding_utils.h"
#include "misc_types.h"


/* H_CODE / S_CODE struct representing one element of a quadrature */
struct H_CODE {
    UINT8_VEC_1D  signs;
    UINT16_VEC_1D m_codes;

    H_CODE() {
        signs.reserve(128);
        m_codes.reserve(128);
    }

    H_CODE(const int& num_codes)
    {
        signs.reserve(num_codes);
        m_codes.reserve(num_codes);
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
        blocks.reserve(num_blocks);
    }
};


/* Class for Sentinel-1 Level-0 Packets with decoding functionality for all packet types. */
class L0Packet 
{
private:
    std::unordered_map<std::string, int> _primary_header;
    std::unordered_map<std::string, int> _secondary_header;

    UINT8_VEC_1D _raw_user_data;

    int _num_quads;
    int _test_mode;
    int _baq_mode;
    int _num_baq_blocks;
    int _user_data_length;
    int _bit_index;

    bool _is_empty = true;

    char _data_format;

    void _set_data_format();

    static std::unordered_map<std::string, int> _parse_header(
        const UINT8_VEC_1D&  bytes,
        const INT_VEC_1D&    bit_lengths,
        const STRING_VEC_1D& field_names
    );


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
    int  _get_next_word_boundary(const int& bit_index);

    H_CODE _get_h_code_type_c(
              int&  bit_index, 
        const bool& is_last_block
    );
    H_CODE _get_h_code_type_d(
        const u_int8_t& brc,
              int&      bit_index,
        const bool&     is_last_block
    );
    double _get_s_values_type_c(
        const u_int16_t& threshold_id,
        const int&       sign,
        const int&       m_code
    );
    double _get_s_values_type_d(
        const u_int8_t&  brc,
        const u_int16_t& threshold_id,
        const int&       sign,
        const int&       m_code
    );

    void _set_quad_type_d(QUAD& component, int& bit_index);
    void _set_quad_type_c(QUAD& component, int& bit_index);
    void _set_quad_types_a_and_b(H_CODE& component, int& bit_index);

    CF_VEC_1D _get_signal_type_d(
        QUAD& IE, QUAD& IO, 
        QUAD& QE, QUAD& QO
    );
    CF_VEC_1D _get_signal_type_c(
        QUAD& IE, QUAD& IO,
        QUAD& QE, QUAD& QO
    );
    CF_VEC_1D _get_signal_types_a_and_b(
        H_CODE& IE, H_CODE& IO,
        H_CODE& QE, H_CODE& QO
    );


public:
    L0Packet(){}    

    L0Packet(
        std::unordered_map<std::string, int> primary_header,
        std::unordered_map<std::string, int> secondary_header,
        UINT8_VEC_1D raw_user_data
    ) {
        _primary_header   = primary_header;
        _secondary_header = secondary_header;
        _raw_user_data    = raw_user_data;

        _num_quads        = secondary_header["num_quadratures"];
        _test_mode        = secondary_header["test_mode"];
        _baq_mode         = secondary_header["baq_mode"];
        _user_data_length = primary_header["packet_data_length"] + 1 - SECONDARY_HEADER_SIZE;
        _num_baq_blocks   = ceil((2.0 * double(_num_quads)) / 256.0);

        if (_user_data_length != _raw_user_data.size()) 
        {
            std::cout << _user_data_length << " != " << _raw_user_data.size() << std::endl;
            throw std::runtime_error("The lenght of the user data field is invalid.");
        }

        _set_data_format();

        _is_empty = false;
    }

    bool is_empty() {return _is_empty;}

    int  get_num_quads() {return _num_quads;}
    int  get_num_baq_blocks  () {return _num_baq_blocks;}
    int  get_user_data_length() {return _user_data_length;}
    char get_data_format() {return _data_format;}

    const int primary_header(const std::string& key) {return _primary_header.at(key);}
    const int secondary_header(const std::string& key) {return _secondary_header.at(key);}

    int get_baq_block_length();

    double get_pulse_length();
    double get_tx_ramp_rate();
    double get_start_frequency();
    double get_pri();
    double get_swl();
    double get_swst();
    double get_rx_gain();

    char get_rx_polarization();
    char get_tx_polarization();

    std::string get_baq_mode();
    std::string get_test_mode();
    std::string get_sensor_mode();
    std::string get_signal_type();
    std::string get_error_status();
    std::string get_swath();

    void print_primary_header();
    void print_secondary_header();
    void print_modes();
    void print_pulse_info();

    CF_VEC_1D get_signal();
    CF_VEC_1D get_replica_chirp();

    static L0Packet get_next_packet(std::ifstream& data);
    static std::vector<L0Packet> get_packets(std::ifstream& data, const int& num_packets = 0);
    static std::vector<L0Packet> get_packets(const std::string& filename, const int& num_packets = 0);
    static std::vector<L0Packet> get_packets_in_swath(const std::string& filename, const std::string& swath);
    static std::vector<L0Packet> get_packets_in_swath(std::ifstream& data, const std::string& swath);
    static std::vector<std::vector<L0Packet>> get_packets_in_bursts(const std::string& filename, const std::string& swath);
    static std::vector<std::vector<L0Packet>> get_packets_in_bursts(std::ifstream& data, const std::string& swath);
    static std::vector<L0Packet> decode_packets(const std::vector<L0Packet>& packets);
    static void decode_packets_in_place(std::vector<L0Packet>& packets);
};


typedef struct std::vector<L0Packet>              PACKET_VEC_1D;
typedef struct std::vector<std::vector<L0Packet>> PACKET_VEC_2D;
