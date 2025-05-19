/*
By: Andrew Player
Name: packet.cpp
Description: L0Packet class for storing and decoding Level-0 Packets in a convinient and easy to use manner. 
             See "SAR Space Packet Protocol Data Unit", for more information on the packet specification:
             https://sentinels.copernicus.eu/documents/247904/2142675/Sentinel-1-SAR-Space-Packet-Protocol-Data-Unit.pdf
             For additional information on Level-0 product decoding, see:
             https://sentinel.esa.int/documents/247904/0/Sentinel-1-Level-0-Data-Decoding-Package.pdf/a8742c59-4914-40c4-8309-c77515649f17
*/

#include "packet.h"


/* Returns the length of the baq blocks in bytes */
const int L0Packet::get_baq_block_length()
{   
    return 8 * (_secondary_header["baq_block_length"] + 1);
}


/* Returns the index of the packet within the data */
const int L0Packet::get_packet_index()
{
    return _packet_index;
}


/* Returns the PRI in microseconds */
const double L0Packet::get_pri()
{
    return _secondary_header.at("pri") / F_REF;
}


/* Returns the pulse length in microseconds */
const double L0Packet::get_pulse_length()
{
    return _secondary_header["pulse_length"] / F_REF;
}


/* Returns the sampling window length in microseconds */
const double L0Packet::get_swl()
{
    return _secondary_header["swl"] / F_REF;
}


/* Returns the start time of the sampling window in the PRI in microseconds */
const double L0Packet::get_swst()
{
    return _secondary_header["swst"] / F_REF;
}


const double L0Packet::get_time()
{
    double coarse_time = static_cast<double>(_secondary_header["coarse_time"]);
    double fine_time   = static_cast<double>(_secondary_header["fine_time"]);

    double time = coarse_time + (1.0 / fine_time);
    return time;
}


/* Returns the RX gain in dB */
const double L0Packet::get_rx_gain()
{
    return -0.5 * _secondary_header["rx_gain"];
}


/* Returns the tx pulse start frequency in MHz */
const double L0Packet::get_start_frequency()
{
    int sign = _secondary_header["pulse_start_frequency_sign"] == 0 ? -1 : 1;
    int mag  = _secondary_header["pulse_start_frequency_mag"];
    double txprr = get_tx_ramp_rate();
    return (sign * mag * (F_REF / 16384)) + (txprr / (4 * F_REF));
}


static const F_VEC_1D AZIMUTH_BEAM_ADDRESS_TO_ANGLE = linspace(-0.018, 0.018, 1024);


/* Returns the azimuth beam angle in radians*/
const double L0Packet::get_azimuth_beam_angle()
{
    return AZIMUTH_BEAM_ADDRESS_TO_ANGLE[_secondary_header["azimuth_beam_address"]];
}


/* Returns the linear FM rate at which the chirp frequency changes in MHz/microsecond */
const double L0Packet::get_tx_ramp_rate()
{
    int sign = _secondary_header["tx_ramp_rate_sign"] == 0 ? -1 : 1;
    int mag  = _secondary_header["tx_ramp_rate_mag"];

    return sign * mag * (pow(F_REF, 2) / 2097152);
}


/* Returns the polarization of the tx pulse */
const char L0Packet::get_tx_polarization()
{
    int pol_code = _secondary_header["polarisation"];

    if      (pol_code >= 0 && pol_code <= 3) return 'H';
    else if (pol_code >= 4 && pol_code <= 7) return 'V';
    else throw std::runtime_error("The polarization code is invalid.");
}


/* Returns the rx polarization */
const char L0Packet::get_rx_polarization()
{
    switch (_secondary_header["rx_channel_id"])
    {
        case 0: return 'V';
        case 1: return 'H';
        default: throw std::runtime_error("The rx_channel_id is invalid.");
    }
}


/* Returns the parity error status of the SES SSB message */
std::string L0Packet::get_error_status()
{
    return _secondary_header.at("error_flag") == 0 ? "nominal" : "ssb_corrupt";
}


/* Returns the measurement, test, or rf characterization mode */
std::string L0Packet::get_sensor_mode()
{
    return ECC_CODE_TO_SENSOR_MODE[_secondary_header["ecc_number"]];
}


/* Returns the swath string representation */
std::string L0Packet::get_swath()
{
    return SWATH_NUM_TO_STRING.at(_secondary_header["swath_number"]);
}


/* Returns the type of baq compression */
std::string L0Packet::get_baq_mode()
{
    std::string baq_mode = std::unordered_map<int, std::string>({
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
std::string L0Packet::get_test_mode()
{
    std::string test_mode = std::unordered_map<int, std::string>({
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
        throw std::out_of_range("The test mode number is not in the range of valid values.");
    }
    return test_mode;
}


/* Returns the type of signal that is being received */
std::string L0Packet::get_signal_type()
{
    std::string signal_type = std::unordered_map<int, std::string>({
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
    for (std::string key : PRIMARY_HEADER_FIELDS) 
    {
        std::cout << key << ": " << _primary_header.at(key) << std::endl;
    }
}


/* Prints the keys and values of the secondary header with no calculations applied */
void L0Packet::print_secondary_header() 
{
    for (std::string key : SECONDARY_HEADER_FIELDS) 
    {
        if (key.substr(0, 3) == "na_")
        {
            continue;
        }
        std::cout << key << ": " << _secondary_header.at(key) << std::endl;
    }
}


/* Prints information related to the operating modes */
void L0Packet::print_modes()
{
    std::cout << "Data Format: "      << get_data_format()      << std::endl;
    std::cout << "BAQ Mode: "         << get_baq_mode()         << std::endl;
    std::cout << "BAQ Block Length: " << get_baq_block_length() << std::endl;
    std::cout << "Num BAQ Blocks: "   << _num_baq_blocks        << std::endl;
    std::cout << "Test Mode: "        << get_test_mode()        << std::endl;
    std::cout << "Sensor Mode: "      << get_sensor_mode()      << std::endl;
    std::cout << "Signal Type: "      << get_signal_type()      << std::endl;
    std::cout << "Error Status: "     << get_error_status()     << std::endl;
}


/* Prints the pulse descriptors with the relevant calculations applied */
void L0Packet::print_pulse_info() 
{
    int range_decimation = _secondary_header["range_decimation"];
    int tx_pulse_number  = _secondary_header["tx_pulse_number"];

    std::cout << std::fixed                  << std::setprecision(8);
    std::cout << "Time: "                    << get_time()            << std::endl;
    std::cout << "Swath: "                   << get_swath()           << std::endl;
    std::cout << "RX Polarization: "         << get_rx_polarization() << std::endl;
    std::cout << "TX Polarization: "         << get_tx_polarization() << std::endl;
    std::cout << "Pulse Length: "            << get_pulse_length()    << std::endl;
    std::cout << "TX Ramp Rate (TXPRR): "    << get_tx_ramp_rate()    << std::endl;
    std::cout << "Start Frequency (TXPSF): " << get_start_frequency() << std::endl;
    std::cout << "PRI: "                     << get_pri()             << std::endl;
    std::cout << "SWL: "                     << get_swl()             << std::endl;
    std::cout << "SWST: "                    << get_swst()            << std::endl;
    std::cout << "RX Gain: "                 << get_rx_gain()         << std::endl;
    std::cout << "Range Decimation: "        << range_decimation      << std::endl;
    std::cout << "TX Pulse Number: "         << tx_pulse_number       << std::endl;
}


/* Sets the data format of the packet. The data format determines what type of decoding to do. */
void L0Packet::_set_data_format() 
{
    std::unordered_set<int> type_c_modes = { 3,  4,  5};
    std::unordered_set<int> type_d_modes = {12, 13, 14};

    if (_baq_mode == 0) _test_mode % 3 == 0 ?  _data_format = 'A' : _data_format = 'B';
    else if (type_c_modes.contains(_baq_mode)) _data_format = 'C';
    else if (type_d_modes.contains(_baq_mode)) _data_format = 'D';
    else throw std::out_of_range("BAQ Mode is invalid.");
}


double L0Packet::get_range_sample_rate()
{
    int range_dec = secondary_header("range_decimation");
    return RANGE_DECIMATION[range_dec];
}


D_VEC_1D L0Packet::get_slant_ranges(int num_ranges)
{
    if (num_ranges <= 0) num_ranges = 4 * _num_quads;

    double txpsf = get_start_frequency();
    double txprr = get_tx_ramp_rate();
    double txpl  = get_pulse_length() * 1e-6;

    double start_time = get_swst() * 1e-6;
    double pri = get_pri() *1e-6;
    double rank = _secondary_header.at("rank");

    double delta_t = (320 / (8 * F_REF)) * 1e-6;

    double delay = rank * pri + start_time + delta_t;

    double min_slant_range = delay * SPEED_OF_LIGHT / 2;
    double max_slant_range = (delay + txpl) * SPEED_OF_LIGHT / 2;

    return linspace(min_slant_range, max_slant_range, num_ranges);
}



D_VEC_1D L0Packet::get_slant_range_times(int num_ranges)
{
    if (num_ranges <= 0) num_ranges = 4 * _num_quads;

    double txpsf = get_start_frequency();
    double txprr = get_tx_ramp_rate();
    double txpl  = get_pulse_length() * 1e-6;

    double start_time = get_swst() * 1e-6;
    double pri = get_pri() *1e-6;
    double rank = _secondary_header.at("rank");

    double delta_t = (320 / (8 * F_REF)) * 1e-6;

    double delay = rank * pri + start_time + delta_t;

    return linspace(delay, delay + txpl, num_ranges);
}


CF_VEC_1D L0Packet::get_replica_chirp()
{
    int num_range = 2 * _num_quads;

    double txpsf = get_start_frequency();
    double txprr = get_tx_ramp_rate();
    double txpl  = get_pulse_length();

    double phi_1 = txpsf;
    double phi_2 = txprr * 0.5;

    int range_dec   = secondary_header("range_decimation");
    int num_samples = int(floor(RANGE_DECIMATION[range_dec] * txpl));

    F_VEC_1D  time = linspace(0.0, txpl, num_samples);

    int min_index = int(ceil((num_range - num_samples)/2))-1;
    int max_index = min_index + num_samples;

    CF_VEC_1D chirp(num_range);
    for (int i = min_index; i < max_index; i++)
    {
        double t  = time[i - min_index]; 
        chirp[i] = double(1.0 / num_samples) * exp(I * 2.0 * PI * ((phi_1 * t) + phi_2 * (t * t)));
    }

    return chirp;
}



/***********************************************************************/
/*                                                                     */
/* DECODING COMPLEX SAMPLES                                            */
/*                                                                     */
/* PAGES 61 -> 85                                                      */
/*                                                                     */
/***********************************************************************/


/* Returns the decoded complex sample data - only does the calculations the first time.  */
CF_VEC_1D L0Packet::get_signal()
{
    if (!_signal_set_flag)
    {
        _set_signal();
    }
    return _signal;
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

        _signal = _get_signal_types_a_and_b(IE, IO, QE, QO);

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

        _signal = _get_signal_type_c(IE, IO, QE, QO);
    }
    else if (_data_format == 'D')
    {
        _brc.resize(_num_baq_blocks);
        _thresholds.resize(_num_baq_blocks);

        _set_quad_type_d(IE, bit_index);
        _set_quad_type_d(IO, bit_index);
        _set_quad_type_d(QE, bit_index);
        _set_quad_type_d(QO, bit_index);

        _signal = _get_signal_type_d(IE, IO, QE, QO);
    }
}


/* Setter for _complex_samples - destroys _raw_user_data to free space */
void L0Packet::_set_signal()
{
    _signal.reserve(_num_quads*4);
    _decode();
    _signal_set_flag = true;

    UINT8_VEC_1D().swap(_raw_user_data);
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


CF_VEC_1D L0Packet::_get_signal_types_a_and_b(
    H_CODE& IE,
    H_CODE& IO,
    H_CODE& QE,
    H_CODE& QO
) {
    auto get_s_value = [](u_int8_t sign, u_int16_t m_code) {return pow(-1, sign) * m_code;};

    D_VEC_2D s_values;
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
    CF_VEC_1D complex_samples;
    complex_samples.reserve(_num_quads*4);

    for (int i = 1; i <= _num_quads; i++)
    {
        D_VEC_1D components = s_values[i-1];

        complex_samples.push_back(std::complex<double>(components[0], components[2]));
        complex_samples.push_back(std::complex<double>(components[1], components[3]));
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
            else if (m_code == flag) return pow(-1, sign) * SIMPLE_RECONSTRUCTION[0][_baq_mode-3][threshold_index];
            else throw std::runtime_error("Invalid m_code in s_value generation.");
        }
        double sigma_factor              = THIDX_TO_SF[threshold_index];
        double norm_reconstruction_level = NORMALIZED_RECONSTRUCTION[0][_baq_mode - 3][m_code];
        return pow(-1.0, sign) * norm_reconstruction_level * sigma_factor;
}


CF_VEC_1D L0Packet::_get_signal_type_c(
    QUAD& IE,
    QUAD& IO,
    QUAD& QE,
    QUAD& QO
) {
    D_VEC_2D s_values;
    s_values.reserve(_num_quads);

    for (int block_id = 0; block_id < _num_baq_blocks; block_id++)
    {
        bool is_last_block = (block_id == _num_baq_blocks - 1);
        int  block_length  = is_last_block ? _num_quads - (128 * (_num_baq_blocks - 1)) : 128;
        int  threshold_id  = _thresholds[block_id];

        for (int s_id = 0; s_id < block_length; s_id++)
        {
            s_values.push_back({
                _get_s_values_type_c(
                    threshold_id,
                    IE.blocks[block_id].signs[s_id],
                    IE.blocks[block_id].m_codes[s_id]
                ),
                _get_s_values_type_c(
                    threshold_id,
                    IO.blocks[block_id].signs[s_id],
                    IO.blocks[block_id].m_codes[s_id]
                ),
                _get_s_values_type_c(
                    threshold_id,
                    QE.blocks[block_id].signs[s_id],
                    QE.blocks[block_id].m_codes[s_id]
                ),
                _get_s_values_type_c(
                    threshold_id,
                    QO.blocks[block_id].signs[s_id],
                    QO.blocks[block_id].m_codes[s_id]
                )
            });
        }
    }
    CF_VEC_1D complex_samples;
    complex_samples.reserve(_num_quads*4);

    for (int i = 1; i <= _num_quads; i++)
    {
        D_VEC_1D components = s_values[i-1];

        complex_samples.push_back(std::complex<double>(components[0], components[2]));
        complex_samples.push_back(std::complex<double>(components[1], components[3]));
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
                throw std::runtime_error("Threshold Index is invalid.");
            }
            _thresholds.push_back(threshold);
            bit_index += threshold_bits;
        }
        component.blocks.push_back(_get_h_code_type_c(bit_index, is_last_block));
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
        else if (m_code == flag) return pow(-1.0, sign) * SIMPLE_RECONSTRUCTION[1][brc][threshold_index];
        else throw std::runtime_error("MCode is greater than the comparison flag.");
    }
    return pow(-1.0, sign) * NORMALIZED_RECONSTRUCTION[1][brc][m_code] * THIDX_TO_SF[threshold_index];
}


CF_VEC_1D L0Packet::_get_signal_type_d(
    QUAD& IE,
    QUAD& IO,
    QUAD& QE,
    QUAD& QO
) {
    D_VEC_2D s_values(_num_quads, D_VEC_1D(4));

    for (int block_id = 0; block_id < _num_baq_blocks; block_id++)
    {
        bool is_last_block = (block_id == _num_baq_blocks - 1);
        int  block_length  = is_last_block ? _num_quads - (128 * (_num_baq_blocks - 1)) : 128;
        int  brc           = _brc[block_id];
        int  threshold_id  = _thresholds[block_id];

        for (int s_id = 0; s_id < block_length; s_id++)
        {
            s_values[block_id*block_length+s_id] = D_VEC_1D({
                _get_s_values_type_d(
                    brc, threshold_id, 
                    IE.blocks[block_id].signs[s_id],
                    IE.blocks[block_id].m_codes[s_id]
                ),
                _get_s_values_type_d(
                    brc, threshold_id,
                    IO.blocks[block_id].signs[s_id],
                    IO.blocks[block_id].m_codes[s_id]
                ),
                _get_s_values_type_d(
                    brc, threshold_id,
                    QE.blocks[block_id].signs[s_id],
                    QE.blocks[block_id].m_codes[s_id]
                ),
                _get_s_values_type_d(
                    brc, threshold_id,
                    QO.blocks[block_id].signs[s_id],
                    QO.blocks[block_id].m_codes[s_id]
                )
            });
        }
    }
    CF_VEC_1D complex_samples(_num_quads * 2);

    for (int i = 1; i <= _num_quads*2; i+=2)
    {
        int s_index = ceil((i-1)/2);
        
        D_VEC_1D components = s_values[s_index];

        complex_samples[i-1] = std::complex<double>(components[0], components[2]);
        complex_samples[i] = std::complex<double>(components[1], components[3]);
    }
    return complex_samples;
}


H_CODE L0Packet::_get_h_code_type_d(
    const u_int8_t& brc,
          int&      bit_index,
    const bool&     is_last_block
) {
    int num_codes = is_last_block ? _num_quads - (128 * (_num_baq_blocks - 1)) : 128;

    H_CODE h_code(num_codes);

    for (int i = 0; i < num_codes; i++)
    {
        int start_bit_index = bit_index;
        int sign   = read_n_bits(_raw_user_data, bit_index, 1);

        bit_index += 1;

        u_int16_t m_code = huffman_decode(_raw_user_data, brc, bit_index);

        h_code.signs[i] = sign;
        h_code.m_codes[i] = m_code;
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
                throw std::runtime_error("BRC value is invalid.");
            }
            _brc[i] = brc;
            bit_index += brc_bits;
        }
        else if (is_qe)
        {
            threshold = read_n_bits(_raw_user_data, bit_index, threshold_bits);
    
            if (threshold > 256)
            {
                throw std::runtime_error("Threshold Index is invalid.");
            }
            _thresholds[i] = threshold;
            bit_index += threshold_bits;
        }
        brc = _brc[i];

        component.blocks[i] = _get_h_code_type_d(brc, bit_index, is_last_block);
    }
    bit_index = _get_next_word_boundary(bit_index);
}


/***********************************************************************/
/*                                                                     */
/* STATIC DECODING METHODS                                             */
/*                                                                     */
/***********************************************************************/


/* Returns the header dictionary with the values cast to integers with no calculations */
std::unordered_map<std::string, int> L0Packet::_parse_header(
    const UINT8_VEC_1D&  bytes,
    const INT_VEC_1D&   bit_lengths,
    const STRING_VEC_1D& field_names
) {
    int num_fields = field_names.size();
    int bit_index = 0;

    std::unordered_map<std::string, int> header;

    for (int i = 0; i < num_fields; i++) 
    {   
        header[field_names[i]] = read_n_bits(bytes, bit_index, bit_lengths[i]);

        bit_index += bit_lengths[i];
    }
    return header;
}


/* Decode the next packet within the data stream */
L0Packet L0Packet::get_next_packet(std::ifstream& data, int& packet_index)
{
    UINT8_VEC_1D primary_bytes = read_bytes(data, 6);
    std::unordered_map<std::string, int> primary_header = _parse_header(
        primary_bytes,
        PRIMARY_HEADER,
        PRIMARY_HEADER_FIELDS
    );
    if (data.eof()) return L0Packet();

    UINT8_VEC_1D secondary_bytes = read_bytes(data, 62);
    std::unordered_map<std::string, int> secondary_header = _parse_header(
        secondary_bytes,
        SECONDARY_HEADER,
        SECONDARY_HEADER_FIELDS
    );

    u_int32_t packet_length    = primary_header["packet_data_length"];
    u_int32_t user_data_length = packet_length + 1 - SECONDARY_HEADER_SIZE;

    UINT8_VEC_1D user_data = read_bytes(data, user_data_length);

    L0Packet packet = L0Packet(
        primary_header,
        secondary_header,
        user_data,
        packet_index
    );
    return packet;
}


/* Returns num_packets packets from the data stream or all packets if num_packets is 0  */
PACKET_VEC_1D L0Packet::get_packets(const std::string& filename, const int& num_packets)
{
    std::ifstream data = open_file(filename);
    return get_packets(data, num_packets);
}


/* Returns num_packets packets from the data stream or all packets if num_packets is 0  */
PACKET_VEC_1D L0Packet::get_packets(std::ifstream& data, const int& num_packets)
{
    PACKET_VEC_1D packets;

    int index = 0;

    bool get_all_packets = num_packets == 0;

    while (!data.eof() and (index < num_packets or get_all_packets))
    {
        try
        {
            L0Packet packet = L0Packet::get_next_packet(data, index);

            if (!packet.is_empty()) packets.push_back(packet);
            else break;
            index += 1;
        }
        catch(std::runtime_error)
        {
            std::cout << "Caught a runtime error while decoding packet #"
                 << index << ". Skipping..." << std::endl;
            continue;
        }
    }
    return packets;
}


/* Returns all packets that are inside of the provided swath */
PACKET_VEC_1D L0Packet::get_packets_in_swath(const std::string& filename, const std::string& swath)
{
    std::ifstream data = open_file(filename);
    return get_packets_in_swath(data, swath);
}


/* Returns all packets that are inside of the provided swath */
PACKET_VEC_1D L0Packet::get_packets_in_swath(std::ifstream& data, const std::string& swath)
{
    PACKET_VEC_1D packets;

    int index = 0;

    while (!data.eof())
    {
        try
        {
            L0Packet packet = L0Packet::get_next_packet(data, index);

            if (packet.is_empty()) break;
            if (packet.get_swath() == swath) packets.push_back(packet);
            index += 1;
        }
        catch(std::runtime_error)
        {
            std::cout << "Caught a runtime error while decoding packet #" 
                 << index << ". Skipping..." << std::endl;
            continue;
        }
    }

    if (packets.size() == 0)
    {
        throw std::runtime_error("No packets found for swath: " + swath);
    }

    return packets;
}


/* Returns all packets in the provided swath, with each burst in its own vector. */
PACKET_VEC_2D L0Packet::get_packets_in_bursts(const std::string& filename, const std::string& swath)
{
    std::ifstream data = open_file(filename);
    return get_packets_in_bursts(data, swath);
}


/* Returns all packets in the provided swath, with each burst in its own vector. */
PACKET_VEC_2D L0Packet::get_packets_in_bursts(std::ifstream& data, const std::string& swath)
{
    PACKET_VEC_1D packets = L0Packet::get_packets_in_swath(data, swath);
    int num_packets = packets.size();

    PACKET_VEC_2D bursts; 
    PACKET_VEC_1D burst_packets;
    int previous_az = 0;

    for (int i = 0; i < num_packets; i++)
    {
        L0Packet packet = packets[i];
        if (packet.get_data_format() == 'D')
        {
            int az = packet.secondary_header("azimuth_beam_address");

            if (i == 0) previous_az = az;

            if (az != previous_az and az != previous_az + 1)
            {
                bursts.push_back(burst_packets);
                burst_packets = PACKET_VEC_1D();
            }
            burst_packets.push_back(packet);

            previous_az = az;
        }
        if (i == num_packets - 1) bursts.push_back(burst_packets);
    }

    return bursts;
}


/* Returns all packets in the provided swath, with each burst in its own vector. */
PACKET_VEC_2D L0Packet::get_packets_in_bursts(PACKET_VEC_1D& packets, const std::string& swath, const bool& get_cal_packets)
{
    int num_packets = packets.size();

    PACKET_VEC_2D bursts; 
    PACKET_VEC_1D burst_packets;
    int previous_az = 0;
    int previous_pri_count = 0;

    for (int i = 0; i < num_packets; i++)
    {

        L0Packet packet = packets[i];
        bool type_check = get_cal_packets ? packet.get_data_format() != 'D' : packet.get_data_format() == 'D';

        if (type_check and packet.get_swath() == swath)
        {
            int az = packet.secondary_header("azimuth_beam_address");
            int pri_count = packet.secondary_header("pri_count");
            if (previous_az == 0) 
            {
                previous_az = az;
                previous_pri_count = pri_count;
            }
            if (az < previous_az or pri_count > previous_pri_count + 1)
            {
                if (burst_packets.size() < 1000)
                {
                    burst_packets = PACKET_VEC_1D();
                }
                else
                {
                    bursts.push_back(burst_packets);
                    burst_packets = PACKET_VEC_1D();
                }
            }
            burst_packets.push_back(packet);

            previous_pri_count = pri_count;
            previous_az = az;
        }
        if (i == num_packets - 1 && burst_packets.size() > 0) bursts.push_back(burst_packets);
    }

    return bursts;
}


PACKET_VEC_1D L0Packet::decode_packets(const PACKET_VEC_1D& packets)
{
    int num_packets = packets.size();

    PACKET_VEC_1D packets_out(num_packets);

    #pragma omp parallel for
    for (int i = 0; i < num_packets; i++)
    {
        L0Packet packet = packets[i];
        packet.get_signal();
        packets_out[i] = packet;
    }

    return packets_out;
}


void L0Packet::decode_packets_in_place(PACKET_VEC_1D& packets)
{
    int num_packets = packets.size();
    
    #pragma omp parallel for
    for (int i = 0; i < num_packets; i++)
    {
        packets[i].get_signal();
    }
}
