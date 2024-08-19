/*
By: Andrew Player
Name: packet.cpp
Description: L0Packet class for storing and decoding Level-0 Packets in a convinient and easy to use manner. 
             See "SAR Space Packet Protocol Data Unit", for more information on the packet specification:
             https://sentinels.copernicus.eu/documents/247904/2142675/Sentinel-1-SAR-Space-Packet-Protocol-Data-Unit.pdf
             For additional information on Level-0 product decoding, see:
             https://sentinel.esa.int/documents/247904/0/Sentinel-1-Level-0-Data-Decoding-Package.pdf/a8742c59-4914-40c4-8309-c77515649f17
*/

#include "packet.hpp"


/* Returns the length of the baq blocks in bytes */
int L0Packet::get_baq_block_length()
{   
    return 8 * (_secondary_header["baq_block_length"] + 1);
}


/* Returns the PRI in microseconds */
double L0Packet::get_pri()
{
    return _secondary_header["pri"] / F_REF;
}


/* Returns the pulse length in microseconds */
double L0Packet::get_pulse_length()
{
    return _secondary_header["pulse_length"] / F_REF;
}


/* Returns the sampling window length in microseconds */
double L0Packet::get_swl()
{
    return _secondary_header["swl"] / F_REF;
}


/* Returns the start time of the sampling window in the PRI in microseconds */
double L0Packet::get_swst()
{
    return _secondary_header["swst"] / F_REF;
}


/* Returns the RX gain in dB */
double L0Packet::get_rx_gain()
{
    return -0.5 * _secondary_header["rx_gain"];
}


/* Returns the tx pulse start frequency in MHz */
double L0Packet::get_start_frequency()
{
    int sign = _secondary_header["pulse_start_frequency_sign"] == 0 ? -1 : 1;
    int mag  = _secondary_header["pulse_start_frequency_mag"];

    double txprr = get_tx_ramp_rate();

    return (sign * mag * (F_REF / 16384)) + (txprr / (4 * F_REF));
}


/* Returns the linear FM rate at which the chirp frequency changes in MHz/microsecond */
double L0Packet::get_tx_ramp_rate()
{
    int sign = _secondary_header["tx_ramp_rate_sign"] == 0 ? -1 : 1;
    int mag  = _secondary_header["tx_ramp_rate_mag"];

    return sign * mag * (pow(F_REF, 2) / 2097152);
}


/* Returns the polarization of the tx pulse */
char L0Packet::get_tx_polarization()
{
    int pol_code = _secondary_header["polarisation"];
    
    if (pol_code >= 0 && pol_code <= 3)
    {
        return 'H';
    }
    else if (pol_code >= 4 && pol_code <= 7)
    {
        return 'V';
    }
    else
    {
        throw runtime_error("The polarization code is invalid.");
    }
}


/* Returns the rx polarization */
char L0Packet::get_rx_polarization()
{
    switch (_secondary_header["rx_channel_id"])
    {
        case 0:
            return 'V';
        case 1:
            return 'H';
        default:
            throw runtime_error("The rx_channel_id is invalid.");
    }
}


/* Returns the parity error status of the SES SSB message */
string L0Packet::get_error_status()
{
    return _secondary_header["error_flag"] == 0 ? "nominal" : "ssb_corrupt";
}


/* Returns the measurement, test, or rf characterization mode */
string L0Packet::get_sensor_mode()
{
    return ECC_CODE_TO_SENSOR_MODE[_secondary_header["ecc_number"]];
}


/* Returns the type of baq compression */
string L0Packet::get_baq_mode()
{
    string baq_mode = unordered_map<int, string>({
        {0,  "bypass_mode"},
        {3,  "3_bit_mode"},
        {4,  "4_bit_mode"},
        {5,  "5_bit_mode"},
        {12, "fdbaq_mode_0"},
        {13, "fdbaq_mode_1"},
        {14, "fdbaq_mode_2"}
    })[_baq_mode];
    
    return baq_mode == "" ? "n_a" : baq_mode;
}


/* Returns the test mode, or measurement mode if not testing */
string L0Packet::get_test_mode()
{
    string test_mode = unordered_map<int, string>({
        {0, "measurement_mode"},
        {1, "n_a"},
        {3, "n_a"},
        {4, "contingency"},
        {5, "contingency"},
        {6, "test_mode_baq"},
        {7, "test_mode_bypass"}
    })[_test_mode];

    if (test_mode == "")
    {
        throw out_of_range("The test mode number is not in the range of valid values.");
    }
    return test_mode;
}


/* Returns the type of signal that is being received */
string L0Packet::get_signal_type()
{
    string signal_type = unordered_map<int, string>({
        {0,  "echo"},
        {1,  "noise"},
        {8,  "tx_cal"},
        {9,  "rx_cal"},
        {10, "epdn_cal"},
        {11, "ta_cal"},
        {12, "apdn_cal"},
        {15, "txh_cal_iso"}
    })[_secondary_header["signal_type"]];
    
    return signal_type == "" ? "n_a" : signal_type;
}


/* Prints the keys and values of the primary header */
void L0Packet::print_primary_header() 
{
    for (string key : PRIMARY_HEADER_FIELDS) 
    {
        cout << key << ": " << _primary_header[key] << endl;
    }
}


/* Prints the keys and values of the secondary header with no calculations applied */
void L0Packet::print_secondary_header() 
{
    for (string key : SECONDARY_HEADER_FIELDS) 
    {
        if (key.substr(0, 3) == "na_")
        {
            continue;
        }
        cout << key << ": " << _secondary_header[key] << endl;
    }
}


/* Prints information related to the operating modes */
void L0Packet::print_modes()
{
    cout << "Data Format: "      << get_data_format()      << endl;
    cout << "BAQ Mode: "         << get_baq_mode()         << endl;
    cout << "BAQ Block Length: " << get_baq_block_length() << endl;
    cout << "Test Mode: "        << get_test_mode()        << endl;
    cout << "Sensor Mode: "      << get_sensor_mode()      << endl;
    cout << "Signal Type: "      << get_signal_type()      << endl;
    cout << "Error Status: "     << get_error_status()     << endl;
}


/* Prints the pulse descriptors with the relevant calculations applied */
void L0Packet::print_pulse_info() 
{
    int range_decimation = _secondary_header["range_decimation"];
    int tx_pulse_number  = _secondary_header["tx_pulse_number"];
    
    cout << "RX Polarization: "         << get_rx_polarization() << endl;
    cout << "TX Polarization: "         << get_tx_polarization() << endl;
    cout << "Pulse Length: "            << get_pulse_length()    << endl;
    cout << "TX Ramp Rate (TXPRR): "    << get_tx_ramp_rate()    << endl;
    cout << "Start Frequency (TXPSF): " << get_start_frequency() << endl;
    cout << "PRI: "                     << get_pri()             << endl;
    cout << "SWL: "                     << get_swl()             << endl;
    cout << "SWST: "                    << get_swst()            << endl;
    cout << "RX Gain: "                 << get_rx_gain()         << endl;
    cout << "Range Decimation: "        << range_decimation      << endl;
    cout << "TX Pulse Number: "         << tx_pulse_number       << endl;
}


/* Sets the data format of the packet. The data format determines what type of decoding to do. */
void L0Packet::_set_data_format() 
{
    unordered_set<int> type_c_modes = { 3,  4,  5};
    unordered_set<int> type_d_modes = {12, 13, 14};

    if (_baq_mode == 0) 
    {
        if (_test_mode % 3 == 0) 
        {
            _data_format = 'A';
        } 
        else 
        {
            _data_format = 'B';
        }
    }
    else if (type_c_modes.contains(_baq_mode)) 
    {
        _data_format = 'C';
    }
    else if (type_d_modes.contains(_baq_mode)) 
    {
        _data_format = 'D';
    }
    else 
    {
        throw out_of_range("BAQ Mode is invalid.");
    }
}


/***********************************************************************/
/*                                                                     */
/* DECODING COMPLEX SAMPLES                                            */
/*                                                                     */
/* PAGES 61 -> 85                                                      */
/*                                                                     */
/***********************************************************************/


/* Returns the decoded complex sample data - only does the calculations the first time.  */
vector<complex<float>> L0Packet::get_complex_samples()
{
    if (!_complex_samples_set_flag)
    {
        _set_complex_samples();
    }
    return _complex_samples;
}


/* Decodes the complex data based on the packets data format, and sets _complex_samples */
void L0Packet::_decode() 
{
    int bit_index = 0;

    if (_data_format == 'A' || _data_format == 'B') 
    {
        H_CODE IE;
        H_CODE IO;
        H_CODE QE;
        H_CODE QO;

        _set_quad_types_a_and_b(IE, bit_index);
        _set_quad_types_a_and_b(IO, bit_index);
        _set_quad_types_a_and_b(QE, bit_index);
        _set_quad_types_a_and_b(QO, bit_index);

        _complex_samples = _get_complex_samples_types_a_and_b(IE, IO, QE, QO);

        return;
    }

    QUAD IE = QUAD("IE", _num_baq_blocks);
    QUAD IO = QUAD("IO", _num_baq_blocks);
    QUAD QE = QUAD("QE", _num_baq_blocks);
    QUAD QO = QUAD("QO", _num_baq_blocks);

    _thresholds.reserve(_num_baq_blocks);

    if (_data_format == 'C')
    {
        _set_quad_type_c(IE, bit_index);
        _set_quad_type_c(IO, bit_index);
        _set_quad_type_c(QE, bit_index);
        _set_quad_type_c(QO, bit_index);

        _complex_samples = _get_complex_samples_type_c(IE, IO, QE, QO);
    }
    else if (_data_format == 'D')
    {
        _brc.reserve(_num_baq_blocks);
        
        _set_quad_type_d(IE, bit_index);
        _set_quad_type_d(IO, bit_index);
        _set_quad_type_d(QE, bit_index);
        _set_quad_type_d(QO, bit_index);

        _complex_samples = _get_complex_samples_type_d(IE, IO, QE, QO);
    }
}


/* Setter for _complex_samples - destroys _raw_user_data to free space */
void L0Packet::_set_complex_samples()
{
    _decode();
    _complex_samples_set_flag = true;

    vector<u_int8_t>().swap(_raw_user_data);
}


/* Gets the bit_index at the next word boundary (needed after a quad block is decoded) */
int L0Packet::_get_next_word_boundary(const int& bit_index)
{
    int offset = bit_index % WORD_SIZE;
    if (offset == 0) return bit_index;
    return bit_index + (WORD_SIZE - offset);
}


/***********************************************************************/
/*                                                                     */
/* TYPE A AND B PACKETS                                                */
/*                                                                     */
/* SEE PAGE 61 FOR THE SPEC                                            */
/*                                                                     */
/***********************************************************************/


vector<complex<float>> L0Packet::_get_complex_samples_types_a_and_b(
    H_CODE& IE,
    H_CODE& IO,
    H_CODE& QE,
    H_CODE& QO
) {
    auto get_s_value = [](u_int8_t sign, u_int16_t m_code) {return pow(-1, sign) * m_code;};

    vector<vector<double>> s_values;

    s_values.reserve(_num_quads);

    for(int i = 0; i < _num_quads; i++)
    {
        s_values.push_back({
            get_s_value(IE.signs[i], IE.m_codes[i]),
            get_s_value(IO.signs[i], IO.m_codes[i]),
            get_s_value(QE.signs[i], QE.m_codes[i]),
            get_s_value(QO.signs[i], QO.m_codes[i])
        });
    }

    vector<complex<float>> complex_samples;

    complex_samples.reserve(_num_quads * 4);

    for (int i = 1; i <= _num_quads; i++)
    {
        vector<double> components = s_values[i-1];

        complex_samples.push_back(complex<float>(components[0], components[2]));
        complex_samples.push_back(complex<float>(components[1], components[3]));
    }

    return complex_samples;
}


void L0Packet::_set_quad_types_a_and_b(H_CODE& component, int& bit_index)
{
    int num_bits = 0;
    int sign_bits = 1;
    int m_code_bits = 9;

    for (int i = 0; i < _num_quads; i++)
    {
        u_int8_t sign   = read_n_bits(_raw_user_data, bit_index, sign_bits);
        bit_index += sign_bits;
        
        u_int16_t m_code = read_n_bits(_raw_user_data, bit_index, m_code_bits);
        bit_index += m_code_bits;

        component.signs.push_back(sign);
        component.m_codes.push_back(m_code);
    }

    bit_index = _get_next_word_boundary(bit_index);
}


/***********************************************************************/
/*                                                                     */
/* TYPE C PACKETS                                                      */
/*                                                                     */
/* SEE PAGES 62 -> 66 FOR THE SPEC                                     */
/*                                                                     */
/***********************************************************************/


double L0Packet::_get_s_values_type_c(
    const u_int16_t& threshold_index,
    const int& sign,
    const int& m_code
) {
        if (threshold_index <= BAQ_MODE_TO_THIDX.at(_baq_mode))
        {
            int flag = BAQ_MODE_TO_M_CODE.at(_baq_mode);

            if      (m_code <  flag) return pow(-1, sign) * double(m_code);
            else if (m_code == flag) return pow(-1, sign) * SIMPLE_RECONSTRUCTION_METHOD[0][_baq_mode - 3][threshold_index];
            else throw runtime_error("Invalid m_code in s_value generation.");
        }
        double sigma_factor              = THIDX_TO_SF_ARRAY[threshold_index];
        double norm_reconstruction_level = NORMALIZED_RECONSTRUCTION_LEVELS[0][_baq_mode - 3][m_code];
        return pow(-1.0, sign) * norm_reconstruction_level * sigma_factor;
}


vector<complex<float>> L0Packet::_get_complex_samples_type_c(
    QUAD& IE,
    QUAD& IO,
    QUAD& QE,
    QUAD& QO
) {
    vector<vector<double>> s_values;

    s_values.reserve(_num_quads);

    for (int block_id = 0; block_id < _num_baq_blocks; block_id++)
    {
        bool is_last_block = (block_id == _num_baq_blocks - 1);
        int  block_length  = is_last_block ? _num_quads - (128 * (_num_baq_blocks - 1)) : 128;
        int  threshold_id  = _thresholds[block_id];

        for (int s_id = 0; s_id < block_length; s_id++)
        {
            s_values.push_back({
                _get_s_values_type_c(threshold_id, IE.signs[block_id][s_id], IE.m_codes[block_id][s_id]),
                _get_s_values_type_c(threshold_id, IO.signs[block_id][s_id], IO.m_codes[block_id][s_id]),
                _get_s_values_type_c(threshold_id, QE.signs[block_id][s_id], QE.m_codes[block_id][s_id]),
                _get_s_values_type_c(threshold_id, QO.signs[block_id][s_id], QO.m_codes[block_id][s_id])
            });
        }
    }
    vector<complex<float>> complex_samples;

    complex_samples.reserve(_num_quads * 4);

    for (int i = 1; i <= _num_quads; i++)
    {
        vector<double> components = s_values[i-1];

        complex_samples.push_back(complex<float>(components[0], components[2]));
        complex_samples.push_back(complex<float>(components[1], components[3]));
    }

    return complex_samples;
}


H_CODE L0Packet::_get_h_code_type_c(int& bit_index, const bool& is_last_block)
{
        int num_codes = is_last_block ? _num_quads - (128 * (_num_baq_blocks - 1)) : 128;
        
        H_CODE h_code(num_codes);

        for (int i = 0; i < num_codes; i++)
        {
            int start_bit_index = bit_index;
            int sign   = read_n_bits(_raw_user_data, bit_index, 1);
            bit_index += 1;

            int m_code = read_n_bits(_raw_user_data, bit_index, _baq_mode - 1);
            bit_index += _baq_mode - 1;

            h_code.signs.push_back(sign);
            h_code.m_codes.push_back(m_code);
        }

        return h_code;
}


void L0Packet::_set_quad_type_c(QUAD& component, int& bit_index)
{
    u_int16_t threshold;

    int threshold_bits = 8;       
    int s_code_bits    = _baq_mode;

    for (int i = 0; i < _num_baq_blocks; i++)
    {
        bool is_qe = (component.key == "QE");
        bool is_last_block = (i == _num_baq_blocks - 1);
        if (is_qe)
        {
            threshold = read_n_bits(_raw_user_data, bit_index, threshold_bits);
            if (threshold > 256)
            {
                throw runtime_error("Threshold Index is invalid.");
            }
            _thresholds.push_back(threshold);
            bit_index += threshold_bits;
        }

        H_CODE h_code = _get_h_code_type_c(bit_index, is_last_block);

        component.signs.push_back(h_code.signs);
        component.m_codes.push_back(h_code.m_codes);
        component.bits_read += h_code.bits_read + (is_qe ? threshold_bits : 0);
    }

    bit_index = _get_next_word_boundary(bit_index);
}


/***********************************************************************/
/*                                                                     */
/* TYPE D PACKETS                                                      */
/*                                                                     */
/* SEE PAGES 67 -> 74 FOR THE SPEC                                     */
/*                                                                     */
/***********************************************************************/


double L0Packet::_get_s_values_type_d(
    const u_int8_t& brc,
    const u_int16_t& threshold_index,
    const int& sign,
    const int& m_code
) {
    if (threshold_index <= BRC_TO_THIDX[brc])
    {
        int flag = BRC_TO_M_CODE[brc];

        if      (m_code <  flag) return pow(-1.0, sign) * m_code;
        else if (m_code == flag) return pow(-1.0, sign) * SIMPLE_RECONSTRUCTION_METHOD[1][brc][threshold_index];
        else throw runtime_error("MCode is greater than the comparison flag.");
    }
    return pow(-1.0, sign) * NORMALIZED_RECONSTRUCTION_LEVELS[1][brc][m_code] * THIDX_TO_SF_ARRAY[threshold_index];
}


vector<complex<float>> L0Packet::_get_complex_samples_type_d(
    QUAD& IE,
    QUAD& IO,
    QUAD& QE,
    QUAD& QO
) {
    vector<vector<double>> s_values;

    s_values.reserve(_num_quads);

    for (int block_id = 0; block_id < _num_baq_blocks; block_id++)
    {
        bool is_last_block = (block_id == _num_baq_blocks - 1);
        int  block_length  = is_last_block ? _num_quads - (128 * (_num_baq_blocks - 1)) : 128;
        int  brc           = _brc[block_id];
        int  threshold_id  = _thresholds[block_id];

        for (int s_id = 0; s_id < block_length; s_id++)
        {
            s_values.push_back({
                _get_s_values_type_d(brc, threshold_id, IE.signs[block_id][s_id], IE.m_codes[block_id][s_id]),
                _get_s_values_type_d(brc, threshold_id, IO.signs[block_id][s_id], IO.m_codes[block_id][s_id]),
                _get_s_values_type_d(brc, threshold_id, QE.signs[block_id][s_id], QE.m_codes[block_id][s_id]),
                _get_s_values_type_d(brc, threshold_id, QO.signs[block_id][s_id], QO.m_codes[block_id][s_id]),
            });
        }
    }
    vector<complex<float>> complex_samples;

    complex_samples.reserve(_num_quads * 4);

    for (int i = 1; i <= _num_quads; i++)
    {
        vector<double> components = s_values[i-1];

        complex_samples.push_back(complex<float>(components[0], components[2]));
        complex_samples.push_back(complex<float>(components[1], components[3]));
    }

    return complex_samples;
}


H_CODE L0Packet::_get_h_code_type_d(const u_int8_t& brc, int& bit_index, const bool& is_last_block)
{
        int num_codes = is_last_block ? _num_quads - (128 * (_num_baq_blocks - 1)) : 128;

        H_CODE h_code(num_codes);

        for (int i = 0; i < num_codes; i++)
        {
            int start_bit_index = bit_index;
            int sign   = read_n_bits(_raw_user_data, bit_index, 1);

            bit_index += 1;

            u_int16_t m_code = huffman_decode(_raw_user_data, brc, bit_index);

            h_code.signs.push_back(sign);
            h_code.m_codes.push_back(m_code);
        }

        return h_code;
}


void L0Packet::_set_quad_type_d(QUAD& component, int& bit_index)
{
    u_int8_t brc;
    u_int16_t threshold;

    int brc_bits       = 3;
    int threshold_bits = 8;       

    for (int i = 0; i < _num_baq_blocks; i++)
    {
        bool is_ie = (component.key == "IE");
        bool is_qe = (component.key == "QE");
        bool is_last_block = (i == _num_baq_blocks - 1);
        if (is_ie)
        {
            brc = read_n_bits(_raw_user_data, bit_index, brc_bits);
            if (brc > 4)
            {
                throw runtime_error("BRC value is invalid.");
            }
            _brc.push_back(brc);
            bit_index += brc_bits;
        }
        else if (is_qe)
        {
            threshold = read_n_bits(_raw_user_data, bit_index, threshold_bits);
            if (threshold > 256)
            {
                throw runtime_error("Threshold Index is invalid.");
            }
            _thresholds.push_back(threshold);
            bit_index += threshold_bits;
        }
        brc = _brc[i];

        H_CODE h_code = _get_h_code_type_d(brc, bit_index, is_last_block);

        component.signs.push_back(h_code.signs);
        component.m_codes.push_back(h_code.m_codes);

        int brc_offset       = is_ie ? brc_bits       : 0;
        int threshold_offset = is_qe ? threshold_bits : 0;

        component.bits_read += h_code.bits_read + brc_offset + threshold_offset;
    }

    bit_index = _get_next_word_boundary(bit_index);
}
