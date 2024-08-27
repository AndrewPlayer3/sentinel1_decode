/*
By: Andrew Player
Name: structs.hpp
Description: Some data structs/tables and constants containing information necessary for decoding Level-0 Products.
             They are all from "SAR Space Packet Protocol Data Unit", which can be found at:
             https://sentinels.copernicus.eu/documents/247904/2142675/Sentinel-1-SAR-Space-Packet-Protocol-Data-Unit.pdf
             For additional information on Level-0 product decoding, see:
             https://sentinel.esa.int/documents/247904/0/Sentinel-1-Level-0-Data-Decoding-Package.pdf/a8742c59-4914-40c4-8309-c77515649f17
*/

#pragma once

#include <vector>
#include <string>
#include <complex>
#include <unordered_map>

#include "misc_types.h"

const std::complex<float> I(0.0f, 1.0f); 

const float PI             = 3.14159235;
const float SPEED_OF_LIGHT = 299792458.0;

const double F_REF = 37.53472224;

const u_int8_t PRIMARY_HEADER_SIZE   = 6;
const u_int8_t SECONDARY_HEADER_SIZE = 62;
const u_int8_t WORD_SIZE             = 16;


// Table 2.4-1 from Page 13
const INT_VEC_1D PRIMARY_HEADER = {
    3,   // Packet Version Number
    1,   // Packet Type
    1,   // Secondary Header Flag
    7,   // Process ID
    4,   // Process Category
    2,   // Sequence Flags
    14,  // Packet Sequence Count 
    16   // Packet Data Length
};


// Table 2.4-1 from Page 13
const STRING_VEC_1D PRIMARY_HEADER_FIELDS = {
    "packet_version_number",
    "packet_type",
    "secondary_header_flag",
    "process_id",
    "process_category",
    "sequence_flags",
    "packet_sequence_count",
    "packet_data_length"
};


// Tables 3.2-1 -> 3.2-19 from Pages 15 -> 54
const INT_VEC_1D SECONDARY_HEADER = {
    32,  // Coarse Time
    16,  // Fine Time
    32,  // Sync Marker
    32,  // Data Take ID
    8,   // ECC Number
    1,   // N/A
    3,   // Test Mode
    4,   // RX Channel ID
    32,  // Instrument Configuration ID
    8,   // Sub-Commutative Ancillary Data Word Index
    16,  // Sub-Commutative Ancillary Data Word
    32,  // Counter Service
    32,  // PRI Count
    1,   // Error Flag
    2,   // N/A
    5,   // BAQ Mode
    8,   // BAQ Block Length
    8,   // N/A
    8,   // Range Decimation
    8,   // RX Gain
    1,   // TX Ramp Rate Sign
    15,  // TX Ramp Rate Mangitude
    1,   // Pulse Start Frequency Sign
    15,  // Pulse Start Frequency Magnitude
    24,  // Pulse Length
    3,   // N/A
    5,   // Rank
    24,  // PRI
    24,  // SWST
    24,  // SWL
    1,   // SAS SSB MESSAGE > SSB Flag - "Imaging Mode" is 0, "Calibration Mode" is 1
    3,   // SAS SSB MESSAGE > Polarisation
    2,   // SAS SSB MESSAGE > Temperature Compensation
    2,   // SAS SSB MESSAGE > N/A
    4,   // SAS SSB MESSAGE > Elevation Beam Address or "SAS Test"
    2,   // SAS SSB MESSAGE > N/A
    10,  // SAS SSB MESSAGE > Azimuth Beam Address or Calibration Beam Address
    2,   // SES SSB MESSAGE > Calibration Mode
    1,   // SES SSB MESSAGE > N/A
    5,   // SES SSB MESSAGE > TX Pulse Number
    4,   // SES SSB MESSAGE > Signal Type
    3,   // SES SSB MESSAGE > N/A
    1,   // SES SSB MESSAGE > Swap
    8,   // SES SSB MESSAGE > Swath Number
    16,  // Number of Quadratures
    8    // N/A
};


// Tables 3.2-1 -> 3.2-19 from Pages 15 -> 54
const STRING_VEC_1D SECONDARY_HEADER_FIELDS = {
    "coarse_time",
    "fine_time",
    "sync_marker",
    "data_take_id",
    "ecc_number",
    "na_1",
    "test_mode",
    "rx_channel_id",
    "instrument_configuration_id",
    "sc_data_word_index",
    "sc_data_word",
    "space_packet_count",
    "pri_count",
    "error_flag",
    "na_2",
    "baq_mode",
    "baq_block_length",
    "na_3",
    "range_decimation",
    "rx_gain",
    "tx_ramp_rate_sign",
    "tx_ramp_rate_mag",
    "pulse_start_frequency_sign",
    "pulse_start_frequency_mag",
    "pulse_length",
    "na_4",
    "rank",
    "pri",
    "swst",
    "swl",
    "ssb_flag",
    "polarisation",
    "temperature_compensation",
    "na_5",
    "elevation_beam_address",
    "na_6",
    "azimuth_beam_address",
    "calibration_mode",
    "na_7",
    "tx_pulse_number",
    "signal_type",
    "na_8",
    "swap",
    "swath_number",
    "num_quadratures",
    "na_9"
};


// Figures 4-7 -> 4-11 from Pages 71 -> 73
const INT_VEC_1D BRC_TO_HUFFMAN_START_BIT_LEN = {1, 1, 1, 2, 2};

// Figures 4-7 -> 4-11 from Pages 71 -> 73
const std::vector<std::unordered_map<u_int16_t, u_int8_t>> HUFFMAN_CODINGS = {
    {
        {0, 0},
        {2, 1},
        {6, 2},
        {7, 3}
    },
    {
        {0,  0},
        {2,  1},
        {6,  2},
        {14, 3},
        {15, 4}
    },
    {
        {0,  0},
        {2,  1},
        {6,  2},
        {14, 3},
        {30, 4},
        {62, 5},
        {63, 6}
    },
    {
        {0,   0},
        {1,   1},
        {2,   2},
        {6,   3},
        {14,  4},
        {30,  5}, 
        {62,  6},
        {126, 7},
        {254, 8},
        {255, 9}
    },
    {
        {0,    0},
        {2,    1},
        {3,    2},
        {4,    3},
        {5,    4},
        {12,   5},
        {13,   6},
        {14,   7},
        {30,   8},
        {62,   9},
        {252, 10},
        {253, 11},
        {508, 12},
        {509, 13},
        {510, 14},
        {511, 15}
    }
};


// Figures 4-7 -> 4-11 from Pages 71 -> 73
const std::vector<std::vector<std::unordered_map<u_int16_t, u_int8_t>>> HUFFMAN_CODINGS_WITH_LENGTH = {
    {
        {{0, 0}},
        {{2, 1}},
        {{6, 2},
         {7, 3}}
    },
    {
        {{0,  0}},
        {{2,  1}},
        {{6,  2}},
        {{14, 3},
         {15, 4}}
    },
    {
        {{0,  0}},
        {{2,  1}},
        {{6,  2}},
        {{14, 3}},
        {{30, 4}},
        {{62, 5},
         {63, 6}}
    },
    {
        {{0,   0}},
        {{0,   0},
         {1,   1},
         {2,   2}},
        {{6,   3}},
        {{14,  4}},
        {{30,  5}}, 
        {{62,  6}},
        {{126, 7}},
        {{254, 8},
         {255, 9}}
    },
    {
        {{0,    0}},
        {{0,    0}},
        {{2,    1},
         {3,    2},
         {4,    3},
         {5,    4}},
        {{12,   5},
         {13,   6},
         {14,   7}},
        {{30,   8}},
        {{62,   9}},
        {{0,    0}},
        {{252, 10},
         {253, 11}},
        {{508, 12},
         {509, 13},
         {510, 14},
         {511, 15}}
    }
};



// Table 5.2-1 from Page 78
const std::unordered_map<int, int> BAQ_MODE_TO_THIDX = {{3, 3}, {4, 7}, {5, 15}};

// Table 5.2-1 from Page 78
const UINT16_VEC_1D BRC_TO_THIDX = {3, 3, 5, 6, 8};

// Figures 4-7 -> 4-11 from Pages 71 -> 73
const std::vector<D_VEC_2D> SIMPLE_RECONSTRUCTION = {
    {   // Values for BAQ Compressed Data ***
        { 3.0000,  3.0000,  3.1200,  3.5500,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},
        { 7.0000,  7.0000,  7.0000,  7.1700,  7.4000,  7.7600,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},
        {15.0000, 15.0000, 15.0000, 15.0000, 15.0000, 15.0000, 15.4400, 15.5600, 16.1100, 16.3800, 16.6500, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},
    },  
    {   // Values for FDBAQ Compress Data ***
        { 3.0000,  3.0000,  3.1600,  3.5300,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},
        { 4.0000,  4.0000,  4.0800,  4.3700,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},
        { 6.0000,  6.0000,  6.0000,  6.1500,  6.5000,  6.8800,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},
        { 9.0000,  9.0000,  9.0000,  9.0000,  9.3600,  9.5000, 10.1000,  0.0000,  0.0000,  0.0000,  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},
        {15.0000, 15.0000, 15.0000, 15.0000, 15.0000, 15.0000, 15.2200, 15.5000, 16.0600,  0.0000,  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},
    }
};


// Table 5.2-1 from Page 78
const std::unordered_map<int, int> BAQ_MODE_TO_M_CODE = {{3, 3}, {4, 5}, {5, 10}};

// Table 5.2-2 from Page 79
const INT_VEC_1D BRC_TO_M_CODE = {3, 4, 6, 9, 15};

// Table 5.2-2 from Page 79
const std::vector<D_VEC_2D> NORMALIZED_RECONSTRUCTION = {
    {   // Values for BAQ Compressed Data ***
        {0.2490, 0.7680, 1.3655, 2.1864, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},
        {0.1290, 0.3900, 0.6601, 0.9471, 1.2623, 1.6261, 2.0793, 2.7467, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},
        {0.0660, 0.1985, 0.3320, 0.4677, 0.6061, 0.7487, 0.8964, 1.0510, 1.2143, 1.3896, 1.5800, 1.7914, 2.0329, 2.3234, 2.6971, 3.2692},
    },
    {   // Values for FDBAQ Compressed Data ***
        {0.3637, 1.0915, 1.8208, 2.6406, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},
        {0.3042, 0.9127, 1.5216, 2.1313, 2.8426, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},
        {0.2305, 0.6916, 1.1528, 1.6140, 2.0754, 2.5369, 3.1191, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},
        {0.1702, 0.5107, 0.8511, 1.1916, 1.5321, 1.8726, 2.2131, 2.5536, 2.8942, 3.3744, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},
        {0.1130, 0.3389, 0.5649, 0.7908, 1.0167, 1.2428, 1.4687, 1.6947, 1.9206, 2.1466, 2.3725, 2.5985, 2.8244, 3.0504, 3.2764, 3.6623},
    }
};


// Table 5.2-3 from Page 80
const D_VEC_1D THIDX_TO_SF = {
      0.00,   0.63,   1.25,   1.88,   2.51,   3.13,   3.76,   4.39,   5.01,   5.64,   6.27,   6.89,   7.52,   8.15,   8.77,   9.40,  10.03,
     10.65,  11.28,  11.91,  12.53,  13.16,  13.79,  14.41,  15.04,  15.67,  16.29,  16.92,  17.55,  18.17,  18.80,  19.43,  20.06,  20.68,
     21.31,  21.93,  22.56,  23.19,  23.81,  24.44,  25.06,  25.69,  26.32,  26.94,  27.57,  28.20,  28.83,  29.45,  30.08,  30.71,  31.33,
     31.96,  32.59,  33.21,  33.84,  34.47,  35.09,  35.72,  36.35,  36.97,  37.60,  38.23,  38.85,  39.48,  40.11,  40.73,  41.36,  41.99,
     42.61,  43.24,  43.87,  44.49,  45.12,  45.75,  46.37,  47.00,  47.63,  48.25,  48.88,  49.51,  50.13,  50.76,  51.39,  52.01,  52.64,
     53.27,  53.89,  54.52,  55.15,  55.77,  56.40,  57.03,  57.65,  58.28,  58.91,  59.54,  60.16,  60.79,  61.41,  62.04,  62.98,  64.24,
     65.49,  66.74,  68.00,  69.25,  70.50,  71.76,  73.01,  74.26,  75.52,  76.77,  78.02,  79.28,  80.53,  81.78,  83.04,  84.29,  85.54, 
     86.80,  88.05,  89.30,  90.56,  91.81,  93.07,  94.32,  95.57,  96.83,  98.08,  99.33, 100.58, 101.84, 103.09, 104.34, 105.60, 106.85,
    108.10, 109.36, 110.61, 111.86, 113.11, 114.37, 115.62, 116.87, 118.13, 119.38, 120.63, 121.89, 123.14, 124.39, 125.65, 126.90, 128.15,
    129.41, 130.66, 131.91, 133.17, 134.42, 135.67, 136.93, 138.18, 139.43, 140.69, 141.94, 143.19, 144.45, 145.70, 146.95, 148.21, 149.46,
    150.71, 151.97, 153.22, 154.47, 155.73, 156.98, 158.23, 159.49, 160.74, 161.99, 163.25, 164.50, 165.75, 167.01, 168.26, 169.51, 170.77,
    172.02, 173.27, 174.53, 175.78, 177.03, 178.29, 179.54, 180.79, 182.05, 183.30, 184.55, 185.81, 187.06, 188.31, 189.57, 190.82, 192.07,
    193.33, 194.58, 195.83, 197.09, 198.34, 199.59, 200.85, 202.10, 203.35, 204.61, 205.86, 207.11, 208.37, 209.62, 210.87, 212.13, 213.38,
    214.63, 215.89, 217.14, 218.39, 219.65, 220.90, 222.15, 223.41, 224.66, 225.91, 227.17, 228.42, 229.67, 230.93, 232.18, 233.43, 234.69,
    235.94, 237.19, 238.45, 239.70, 240.95, 242.21, 243.46, 244.71, 245.97, 247.22, 248.47, 249.73, 250.98, 252.23, 253.49, 254.74, 255.99, 255.99
};


// Table 3.2-4 from Page 19
const STRING_VEC_1D ECC_CODE_TO_SENSOR_MODE = {
    "contingency",
    "stripmap_1",
    "stripmap_2",
    "stripmap_3",
    "stripmap_4",
    "stripmap_5n",
    "stripmap_6",
    "contingency",
    "interferomatric_wide_swath",
    "wave_mode",
    "stripmap-5s",
    "stripmap_1_w_cal",
    "stripmap_2_w_cal",
    "stripmap_3_w_cal",
    "stripmap_4_w_cal",
    "rf_characterization_mode",
    "test_mode",
    "elevation_notch_s3",
    "azimuth_notch_s1",
    "azimuth_notch_s2",
    "azimuth_notch_s3",
    "azimuth_notch_s4",
    "azimuth_notch_s5n",
    "azimuth_notch_s5s",
    "azimuth_notch_s6",
    "stripmap_5n_wo_cal",
    "stripmap_5s_wo_cal",
    "stripmap_6_wo_cal",
    "contingency",
    "contingency",
    "contingency",
    "elevation_notch_s3_wo_cal",
    "extra_wide_swath",
    "azimuth_notch_s1_wo_cal",
    "azimuth_notch_s2_wo_cal",
    "azimuth_notch_s3_wo_cal",
    "contingency",
    "noise_characterization_s1",
    "noise_characterization_s2",
    "noise_characterization_s3",
    "noise_characterization_s4",
    "noise_characterization_s5n",
    "noise_characterization_s5s",
    "noise_characterization_s6",
    "noise_characterization_ew",
    "noise_characterization_iw",
    "noise_characterization_wave",
    "contingency"
};


// Table 3.2-4 from Page 19
const std::unordered_map<int, std::string> SWATH_NUM_TO_STRING = {
    {0, "S1"},
    {1, "S2"},
    {2, "S3"},
    {3, "S4"},
    {4, "S5_N"},
    {5, "S5_S"},
    {6, "S6"},
    {10, "IW1"},
    {11, "IW2"},
    {12, "IW3"},
    {20, "EW1"},
    {21, "EW2"},
    {22, "EW3"},
    {23, "EW4"},
    {24, "EW5"},
    {25, "EW1_100MHz_CAL"},
    {26, "EW2_100MHz_CAL"},
    {27, "EW3_100MHz_CAL"},
    {28, "EW4_100MHz_CAL"},
    {29, "EW5_100MHz_CAL"},
    {30, "WV1"},
    {31, "WV2"},
    {32, "WV1_100MHz_CAL"},
    {33, "WV2_100MHz_CAL"},
    {36, "S1_100MHz_CAL"},
    {37, "S2_100MHz_CAL"},
    {38, "S3_100MHz_CAL"},
    {39, "S4_100MHz_CAL"},
    {40, "S5_100MHz_N_CAL"},
    {41, "S5_100MHz_S_CAL"},
    {42, "S6_100MHz_CAL"},
    {43, "IW1_100MHz_CAL"},
    {44, "IW2_100MHz_CAL"},
    {45, "IW3_100MHz_CAL"},
    {50, "S1_TXCAL"},
    {51, "S2_TXCAL"},
    {52, "S3_TXCAL"},
    {53, "S4_TXCAL"},
    {54, "S5_N_TXCAL"},
    {55, "S5_S_TXCAL"},
    {56, "S6_TXCAL"},
    {60, "IW1_TXCAL"},
    {61, "IW2_TXCAL"},
    {62, "IW3_TXCAL"},
    {70, "EW1_TXCAL"},
    {71, "EW2_TXCAL"},
    {72, "EW3_TXCAL"},
    {73, "EW4_TXCAL"},
    {74, "EW5_TXCAL"},
    {75, "EW1_100MHz_TXCAL"},
    {76, "EW2_100MHz_TXCAL"},
    {77, "EW3_100MHz_TXCAL"},
    {78, "EW4_100MHz_TXCAL"},
    {79, "EW5_100MHz_TXCAL"},
    {80, "WV1_TXCAL"},
    {81, "WV2_TXCAL"},
    {82, "WV1_100MHz_TXCAL"},
    {83, "WV2_100MHz_TXCAL"},
    {86, "S1_100MHz_TXCAL"},
    {87, "S2_100MHz_TXCAL"},
    {88, "S3_100MHz_TXCAL"},
    {89, "S4_100MHz_TXCAL"},
    {90, "S5_100MHz_N_TXCAL"},
    {91, "S5_100MHz_S_TXCAL"},
    {92, "S6_100MHz_TXCAL"},
    {93, "IW1_100MHz_TXCAL"},
    {94, "IW2_100MHz_TXCAL"},
    {95, "IW3_100MHz_TXCAL"},
};


const F_VEC_1D RANGE_DECIMATION = {
    112.6041667,
    100.0925926,
    0.000000000,
    83.41049387,
    66.72839509,
    56.30208336,
    50.04629632,
    25.02314816,
    64.34523813,
    46.91840280,
    17.32371796,
    54.59595962
};