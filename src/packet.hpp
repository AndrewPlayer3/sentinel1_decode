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
#include <complex>

#include "decoding_utils.hpp"


using namespace std;


/* H_CODE or S_CODE struct representing one element of a quadrature */
struct H_CODE {
    vector<u_int8_t> signs;
    vector<u_int16_t> m_codes;

    int bits_read;
};


/* Quadrature struct representing one of the IE, IO, QE, or QO quads */
struct QUAD {
    vector<vector<u_int8_t>> signs;
    vector<vector<u_int16_t>> m_codes;

    string key;

    int bits_read;

    QUAD(const string& component_key) {
        key = component_key;
    }
};


/* Class for Sentinel-1 Level-0 Packets with decoding functionality for all packet types. */
class L0Packet 
{
private:
    unordered_map<string, int> _primary_header;
    unordered_map<string, int> _secondary_header;

    vector<u_int8_t> _raw_user_data;

    int _num_quads;
    int _test_mode;
    int _baq_mode;
    int _num_baq_blocks;
    int _user_data_length;
    int _bit_index;

    bool _is_empty = true;

    char _data_format;

    void _set_data_format();


    /**********************************/
    /* DECODING COMPLEX SAMPLES       */
    /* PAGES 61 -> 85                 */
    /**********************************/

    vector<int> _brc = {};
    vector<int> _thresholds = {};
    vector<complex<double>> _complex_samples;

    bool _complex_samples_set_flag = false;

    void _set_complex_samples();
    void _decode();
    int  _get_next_word_boundary(const int& bit_index);

    H_CODE _get_h_code_type_c(int& bit_index, const bool& is_last_block);
    H_CODE _get_h_code_type_d(const u_int8_t& brc, int& bit_index, const bool& is_last_block);

    double _get_s_values_type_c(const u_int16_t& threshold_id, const int& sign, const int& m_code);
    double _get_s_values_type_d(const u_int8_t& brc, const u_int16_t& threshold_id, const int& sign, const int& m_code);

    void _set_quad_type_d(QUAD& component, int& bit_index);
    void _set_quad_type_c(QUAD& component, int& bit_index);
    void _set_quad_types_a_and_b(QUAD& component, int& bit_index);

    vector<complex<double>> _get_complex_samples_type_d(QUAD& IE, QUAD& IO, QUAD& QE, QUAD& QO);
    vector<complex<double>> _get_complex_samples_type_c(QUAD& IE, QUAD& IO, QUAD& QE, QUAD& QO);
    vector<complex<double>> _get_complex_samples_types_a_and_b(QUAD& IE, QUAD& IO, QUAD& QE, QUAD& QO);


public:
    L0Packet(){}    

    L0Packet(
        unordered_map<string, int> primary_header,
        unordered_map<string, int> secondary_header,
        vector<u_int8_t> raw_user_data
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
            cout << _user_data_length << " != " << _raw_user_data.size() << endl;
            throw runtime_error("The lenght of the user data field is invalid.");
        }

        _set_data_format();

        _is_empty = false;
    }

    bool is_empty() {return _is_empty;}

    int  get_num_quads() {return _num_quads;}
    int  get_num_baq_blocks  () {return _num_baq_blocks;}
    int  get_user_data_length() {return _user_data_length;}
    char get_data_format() {return _data_format;}

    int primary_header(const string& key) {return _primary_header[key];}
    int secondary_header(const string& key) {return _secondary_header[key];}

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

    string get_baq_mode();
    string get_test_mode();
    string get_sensor_mode();
    string get_signal_type();
    string get_error_status();

    void print_primary_header();
    void print_secondary_header();
    void print_modes();
    void print_pulse_info();

    vector<complex<double>> get_complex_samples();
};