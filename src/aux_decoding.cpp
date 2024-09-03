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


std::vector<std::unordered_map<std::string, std::string>> build_data_word_dicts(const std::vector<L0Packet>& packets)
{
    int num_packets = packets.size();
    std::vector<std::unordered_map<std::string, std::string>> data_word_dicts;
    std::unordered_map<std::string, std::string> sub_comm_dict = SUB_COMM_KEY_VAL;
    int initial_data_word_index = 0;
    int sc_data_word_index = 0;

    for (int i = 0; i < num_packets; i++)
    {
        L0Packet packet = packets[i];
        sc_data_word_index = packet.secondary_header("sc_data_word_index");
        if (i == 0)
        {
            initial_data_word_index = sc_data_word_index;
        }
        else
        {
            if (sc_data_word_index == initial_data_word_index)
            {
                data_word_dicts.push_back(sub_comm_dict);
                sub_comm_dict = SUB_COMM_KEY_VAL;
            }
        }
        int data_word = packet.secondary_header("sc_data_word");
        std::string key = SUB_COMM_KEY_POS[sc_data_word_index].first;
        int val = SUB_COMM_KEY_POS[sc_data_word_index].second;
        
    }
}


// def build_data_word_dict(packets):
//     data_word_dicts = []
//     sub_comm_dict = SUB_COMM_KEY_VAL.copy()
    
//     initial_data_word_index = 0
//     sc_data_word_index = 0

//     for i in range(len(packets)):
//         packet = packets[i]
        
//         secondary_header = packet.get_secondary_header()
        
//         sc_data_word_index = secondary_header['sc_data_word_index']

//         if i == 0:
//             initial_data_word_index = sc_data_word_index
//         else:
//             if sc_data_word_index == initial_data_word_index:
//                 data_word_dicts.append(sub_comm_dict)
//                 sub_comm_dict = SUB_COMM_KEY_VAL.copy()

//         data_word =  secondary_header['sc_data_word']

//         key, pos = SUB_COMM_KEY_POS[sc_data_word_index]
//         pos = pos * WORD_SIZE
//         sub_comm_dict[key] = sub_comm_dict[key][0:pos] + data_word + sub_comm_dict[key][pos+WORD_SIZE:]
    
//     if sc_data_word_index != initial_data_word_index:
//         data_word_dicts.append(sub_comm_dict)

//     for i in range(len(data_word_dicts)):
//         data_word_dicts[i]['x_axis_position'] = to_float64(data_word_dicts[i]['x_axis_position'])
//         data_word_dicts[i]['y_axis_position'] = to_float64(data_word_dicts[i]['y_axis_position'])
//         data_word_dicts[i]['z_axis_position'] = to_float64(data_word_dicts[i]['z_axis_position'])
//         data_word_dicts[i]['x_axis_velocity'] = to_float32(data_word_dicts[i]['x_axis_velocity'])
//         data_word_dicts[i]['y_axis_velocity'] = to_float32(data_word_dicts[i]['y_axis_velocity'])
//         data_word_dicts[i]['z_axis_velocity'] = to_float32(data_word_dicts[i]['z_axis_velocity'])
//         data_word_dicts[i]['q0_quaternion']   = to_float32(data_word_dicts[i]['q0_quaternion'])
//         data_word_dicts[i]['q1_quaternion']   = to_float32(data_word_dicts[i]['q1_quaternion'])
//         data_word_dicts[i]['q2_quaternion']   = to_float32(data_word_dicts[i]['q2_quaternion'])
//         data_word_dicts[i]['q3_quaternion']   = to_float32(data_word_dicts[i]['q3_quaternion'])
//         data_word_dicts[i]['omega_x']         = to_float32(data_word_dicts[i]['omega_x'])
//         data_word_dicts[i]['omega_y']         = to_float32(data_word_dicts[i]['omega_y'])
//         data_word_dicts[i]['omega_z']         = to_float32(data_word_dicts[i]['omega_z'])

//     return data_word_dicts
