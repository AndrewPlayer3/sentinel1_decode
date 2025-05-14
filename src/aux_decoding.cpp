/*
By: Andrew Player
Name: packet_decoding.cpp
Description: Funtions for reading packets, index info, and annotation info from the Level-0 dat files. 
             See "SAR Space Packet Protocol Data Unit", for more information on the packet specification:
             https://sentinels.copernicus.eu/documents/247904/2142675/Sentinel-1-SAR-Space-Packet-Protocol-Data-Unit.pdf
             For additional information on Level-0 product decoding, see:
             https://sentinel.esa.int/documents/247904/0/Sentinel-1-Level-0-Data-Decoding-Package.pdf/a8742c59-4914-40c4-8309-c77515649f17
*/

#include "aux_decoding.h"


void print_annotation_record(
    const std::string& filename,
    const int&         index
) {
    std::unordered_map<std::string, u_int64_t> record = annotation_decoder(filename)[index];

    std::cout << "Days Uplink: "           << record["days_ul"]         << std::endl;
    std::cout << "Milliseconds Uplink: "   << record["milliseconds_ul"] << std::endl;
    std::cout << "Microseconds Uplink: "   << record["microseconds_ul"] << std::endl;
    std::cout << "Days Downlink: "         << record["days_dl"]         << std::endl;
    std::cout << "Milliseconds Downlink: " << record["milliseconds_dl"] << std::endl;
    std::cout << "Microseconds Downlink: " << record["microseconds_dl"] << std::endl;
    std::cout << "Packet Length: "         << record["packet_length"]   << std::endl;
    std::cout << "Error Flag: "            << record["error_flag"]      << std::endl;
}


void print_index_records(const std::string& filename)
{
    VEC_UNSORTEDMAP records = index_decoder(filename);

    int index = 0;
    for (std::unordered_map<std::string, u_int64_t> record : records)
    {
        std::cout << "Index Record #"  << index << ":"            << std::endl;
        std::cout << "Date/Time: "     << record["date_time"]     << std::endl;
        std::cout << "Time Delta: "    << record["time_delta"]    << std::endl;
        std::cout << "Data Size: "     << record["data_size"]     << std::endl;
        std::cout << "Unit Offset: "   << record["unit_offset"]   << std::endl;
        std::cout << "Byte Offset: "   << record["byte_offset"]   << std::endl;
        std::cout << "Variable Flag: " << record["variable_flag"] << std::endl;
        std::cout << "Spare Data: "    << record["spare_data"]    << std::endl;
        std::cout << "---"                                        << std::endl;

        index += 1;
    }
}


VEC_UNSORTEDMAP index_decoder(const std::string& filename)
{
    std::ifstream data = open_file(filename);
    return index_decoder(data);
}


/* Returns a vector of maps containing the index records */
VEC_UNSORTEDMAP index_decoder(std::ifstream& data)
{
    int record_size = 36;

    VEC_UNSORTEDMAP index_records;

    while (!data.eof())
    {
        UINT8_VEC_1D record_bytes = read_bytes(data, record_size);

        index_records.push_back(std::unordered_map<std::string, u_int64_t>({
            {"date_time",     read_n_bits(record_bytes,   0, 64)},
            {"time_delta",    read_n_bits(record_bytes,  64, 64)},
            {"data_size",     read_n_bits(record_bytes, 128, 32)},
            {"unit_offset",   read_n_bits(record_bytes, 160, 32)},
            {"byte_offset",   read_n_bits(record_bytes, 192, 64)},
            {"variable_flag", read_n_bits(record_bytes, 256,  8)},
            {"spare_data",    read_n_bits(record_bytes, 264, 24)}
        }));
    }
    return index_records;
}


VEC_UNSORTEDMAP annotation_decoder(const std::string& filename)
{
    std::ifstream data = open_file(filename);
    return annotation_decoder(data);
}


/* Returns a vector of maps containing the annotation records */
VEC_UNSORTEDMAP annotation_decoder(std::ifstream& data)
{
    int record_size = 26;

    VEC_UNSORTEDMAP annotation_records;

    while (!data.eof())
    {
        UINT8_VEC_1D record_bytes = read_bytes(data, record_size);

        annotation_records.push_back(std::unordered_map<std::string, u_int64_t>({
            {"days_ul",             read_n_bits(record_bytes,   0, 16)},
            {"milliseconds_ul",     read_n_bits(record_bytes,  16, 32)},
            {"microseconds_ul",     read_n_bits(record_bytes,  48, 16)},
            {"days_dl",             read_n_bits(record_bytes,  64, 16)},
            {"milliseconds_dl",     read_n_bits(record_bytes,  80, 32)},
            {"microseconds_dl",     read_n_bits(record_bytes, 112, 16)},
            {"packet_length",       read_n_bits(record_bytes, 128, 16)},
            {"num_transfer_frames", read_n_bits(record_bytes, 144, 16)},
            {"error_flag",          read_n_bits(record_bytes, 160,  8)},
            {"bit_1_type",          read_n_bits(record_bytes, 168,  1)},
            {"bit_2_type",          read_n_bits(record_bytes, 169,  2)},
            {"bit_6_type",          read_n_bits(record_bytes, 171,  6)},
            {"spare_field",         read_n_bits(record_bytes, 177,  8)}
        }));
    }
    return annotation_records;
}


/* Returns a vector of maps containing the annotation records */
F_VEC_1D annotation_time_decoder(const std::string& filename)
{
    std::ifstream data = open_file(filename);
    return annotation_time_decoder(data);
}

/* Returns a vector of maps containing the annotation records */
F_VEC_1D annotation_time_decoder(std::ifstream& data)
{
    int record_size = 26;

    // Seconds from January 6th 1980 and January 1st 2000. (MJD2000 Epoch and GPS Epoch)
    double epoch_difference = 7300 * 86400;

    VEC_UNSORTEDMAP annotation_records = annotation_decoder(data);

    F_VEC_1D times;

    for (std::unordered_map<std::string, u_int64_t> record : annotation_records)
    {
        double time = (record["days_ul"] * 86400) + (record["milliseconds_ul"] * 1e-3) + (record["microseconds"] * 1e-6) + epoch_difference;
        times.push_back(time);
    }

    return times;
}