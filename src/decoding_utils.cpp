/*
By: Andrew Player
Name: decoding_utils.cpp
Description: Functions to assist reading and decompressing binary data. 
             See "SAR Space Packet Protocol Data Unit", for more information on the packet specification:
             https://sentinels.copernicus.eu/documents/247904/2142675/Sentinel-1-SAR-Space-Packet-Protocol-Data-Unit.pdf
             For additional information on Level-0 product decoding, see:
             https://sentinel.esa.int/documents/247904/0/Sentinel-1-Level-0-Data-Decoding-Package.pdf/a8742c59-4914-40c4-8309-c77515649f17
*/

#include <iostream>
#include "decoding_utils.h"


/* Returns the next m_code using the relevant huffman tree for the given brc
   Checks length of bits in tree to handle edge cases but is slower than huffman_decode */
u_int16_t huffman_decode_with_length(
    const UINT8_VEC_1D& data,
    const int& brc, 
          int& bit_index
) {
    std::vector<std::unordered_map<u_int16_t, u_int8_t>> huffman_coding = HUFFMAN_CODINGS_WITH_LENGTH[brc];

    u_int16_t bits     = -1;
    u_int8_t  bit_len  = BRC_TO_HUFFMAN_START_BIT_LEN[brc];
    u_int8_t  max_bits = 10;

    for (int i = 0; i < max_bits; i++)
    {
        bits = read_n_bits(data, bit_index, bit_len);

        if (!huffman_coding[bit_len - 1].contains(bits)) bit_len += 1;
        else break;

        if (bit_len > 10) throw std::out_of_range("Max bit length exceeded in Huffman decoding.");
    }
    bit_index += bit_len;

    return huffman_coding[bit_len - 1][bits];
}


/* Returns the next m_code using the relevant huffman tree for the given brc */
u_int16_t huffman_decode(
    const UINT8_VEC_1D& data,
    const int& brc, 
          int& bit_index
) {
    u_int16_t bits     = -1;
    u_int8_t  bit_len  = BRC_TO_HUFFMAN_START_BIT_LEN[brc];
    u_int8_t  max_bits = 10;

    if (brc == 4)
    {
        if (read_n_bits(data, bit_index, bit_len) == 0)
        {
            bit_index += 2;
            return 0;
        }
        bit_len += 1;
    }
    for (int i = 0; i < max_bits; i++)
    {
        bits = read_n_bits(data, bit_index, bit_len);

        if (!HUFFMAN_CODINGS[brc].contains(bits)) bit_len += 1;
        else break;

        if (bit_len > 10) 
        {
            throw std::out_of_range("Max bit length exceeded in Huffman decoding.");
        }
    }
    bit_index += bit_len;

    return HUFFMAN_CODINGS[brc].at(bits);
}


/* Returns a vector of the next num_bytes bytes from a data stream */
UINT8_VEC_1D read_bytes(
    std::ifstream&  data,
    const int& num_bytes
) {
    UINT8_VEC_1D buffer(num_bytes);
    data.read(
        reinterpret_cast<char*>(buffer.data()),
        num_bytes
    );
    return buffer;
}


/* Returns the integer encoded by the n bits starting at start_bit in a vector of bytes */
u_int64_t read_n_bits(
    const UINT8_VEC_1D& data,
    const int& start_bit,
    const int& n
) {
    if (n < 1 || n > 64) 
    {
        throw std::invalid_argument("Invalid number of bits to read. Must be between 1 and 64.");
    }
    int byte_index = start_bit / 8;
    int bit_offset = start_bit % 8;

    if (byte_index >= data.size()) 
    {
        throw std::out_of_range("Start bit is outside of the data vector.");
    }
    u_int64_t result = 0;
    int bits_read = 0;

    while (bits_read < n) 
    {
        if (byte_index >= data.size()) 
        {
            throw std::out_of_range("Byte index exceeded data vector size.");
        }

        int bits_left_in_byte = 8 - bit_offset;
        int bits_to_read = std::min(bits_left_in_byte, n - bits_read);

        u_int8_t mask = (1 << bits_to_read) - 1;
        u_int8_t shifted_data = (data[byte_index] >> (bits_left_in_byte - bits_to_read)) & mask;

        result = (result << bits_to_read) | shifted_data;

        bits_read += bits_to_read;
        bit_offset = 0;
        byte_index++;
    }
    return result;
}


/* Return ifstream for the given file */
std::ifstream open_file(const std::string& filename)
{
    std::ifstream data(filename, std::ios::binary);
    if (!data.is_open()) 
    {
        throw std::runtime_error("Unable to open: " + filename);
    }
    return data;
}


double int_to_ieee754(unsigned long number, int is_double)
{
    int mantissaShift = is_double ? 52 : 23;
    unsigned long exponentMask = is_double ? 0x7FF0000000000000 : 0x7f800000;
    int bias = is_double ? 1023 : 127;
    int signShift = is_double ? 63 : 31;

    int sign = (number >> signShift) & 0x01;
    int exponent = ((number & exponentMask) >> mantissaShift) - bias;

    int power = -1;
    double total = 0.0;
    for ( int i = 0; i < mantissaShift; i++ )
    {
        int calc = (number >> (mantissaShift-i-1)) & 0x01;
        total += calc * pow(2.0, power);
        power--;
    }
    double value = (sign ? -1 : 1) * pow(2.0, exponent) * (total + 1.0);

    return value;
}
