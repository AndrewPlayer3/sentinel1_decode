#include "state_vectors.h"


// Time in seconds since J2000 Epoch.
double parse_time_stamp(u_int64_t raw_time_stamp)
{
    uint coarse_time = (raw_time_stamp >> 24) & 0xffffffff;
    uint raw_fine_time = raw_time_stamp & 0xffffff;
    double fine_time = 0.0;
    for (int i = 23; i >= 0; i--)
    {
        fine_time += (raw_fine_time >> i) * pow(2, i-24);
        raw_fine_time &= 1 << 8;
    }
    return double(coarse_time) + fine_time;
}


SUBCOMM_DICTS build_data_word_dicts(PACKET_VEC_1D& packets)
{
    int num_packets = packets.size();

    SUBCOMM_DICT_INT raw_subcomm_dict = SUB_COMM_KEY_VAL_INT;
    SUBCOMM_DICTS processed_subcomm_dicts;

    int start_index = 0;
    int sc_data_word_index = 0;
    while(packets[++start_index].secondary_header("sc_data_word_index") != 1);  // Skip to first to the first occurance of the 0 index

    for (int i = start_index; i < num_packets; i++)
    {
        L0Packet packet = packets[i];
        sc_data_word_index = packet.secondary_header("sc_data_word_index");
        if (sc_data_word_index == 1 and i != start_index)
        {
            SUBCOMM_DICT_DOUBLE processed_subcomm_dict;
            for(std::pair<std::string, u_int64_t> key_val : raw_subcomm_dict)
            {
                bool is_time = key_val.first == "data_time_stamp" or key_val.first == "pod_data_stamp";
                double val = is_time ? 
                    parse_time_stamp(key_val.second) 
                    : 
                    int_to_ieee754(key_val.second, SUB_COMM_DOUBLE_KEYS.contains(key_val.first));
                processed_subcomm_dict[key_val.first] = val;
            }
            processed_subcomm_dicts.push_back(processed_subcomm_dict);
            raw_subcomm_dict = SUB_COMM_KEY_VAL_INT;
        }
        u_int64_t data_word = packet.secondary_header("sc_data_word");
        std::string key = SUB_COMM_KEY_POS[sc_data_word_index].first;
        int offset = SUB_COMM_KEY_POS[sc_data_word_index].second;
        u_int64_t mask = -1 << (offset * 16);
        raw_subcomm_dict[key] = ((raw_subcomm_dict[key] << 16) & mask) | data_word;
    }

    return processed_subcomm_dicts;
}


void STATE_VECTORS::parse_subcomm_dicts(SUBCOMM_DICTS subcomm_dicts)
{
    std::set<double> unique_times = {0.0};

    for (SUBCOMM_DICT_DOUBLE subcomm_dict : subcomm_dicts)
    {
        STATE_VECTOR state_vector(subcomm_dict);

        double current_time = state_vector.time;
        double previous_time = *unique_times.rbegin();
        bool is_valid_time = current_time > previous_time;

        if (is_valid_time)
        {
            unique_times.insert(current_time);
            times.push_back(state_vector.time);
            positions.push_back(state_vector.position);
            velocities.push_back(state_vector.velocity);
            quaternions.push_back(parse_quaternions(state_vector.quaternions).to_vector(SCALAR_LAST));
            angular_rates.push_back(state_vector.angular_rate);
        }
    }
}


void STATE_VECTORS::print()
{
    std::cout << std::setw(15) << "time"
                << std::setw(15) << "position_x"
                << std::setw(15) << "position_y"
                << std::setw(15) << "position_z"
                << std::setw(15) << "velocity_x"
                << std::setw(15) << "velocity_y"
                << std::setw(15) << "velocity_z"
                << std::setw(15) << "q0_quaternion"
                << std::setw(15) << "q1_quaternion"
                << std::setw(15) << "q2_quaternion"
                << std::setw(15) << "q3_quaternion" << std::endl;

    std::cout << std::string(165, '-') << std::endl;

    for (size_t i = 0; i < positions.size(); ++i) {
        std::cout << std::setw(15) << times[i]
                    << std::setw(15) << positions[i][0]
                    << std::setw(15) << positions[i][1]
                    << std::setw(15) << positions[i][2]
                    << std::setw(15) << velocities[i][0]
                    << std::setw(15) << velocities[i][1]
                    << std::setw(15) << velocities[i][2]
                    << std::setw(15) << quaternions[i][0]
                    << std::setw(15) << quaternions[i][1]
                    << std::setw(15) << quaternions[i][2]
                    << std::setw(15) << quaternions[i][3] << std::endl;
    }
    std::cout << std::string(165, '-') << std::endl;
    std::cout << "State Vector Count: " << times.size() << std::endl;
}


Quaternion get_quaternions_from_rotation_matrix(F_VEC_2D& R)
{
    Quaternion q;

    double trace = R[0][0] + R[1][1] + R[2][2];

    if (trace > 0) 
    {
        double S = sqrt(trace + 1.0) * 2;
        q.w = 0.25 * S;
        q.x = (R[2][1] - R[1][2]) / S;
        q.y = (R[0][2] - R[2][0]) / S;
        q.z = (R[1][0] - R[0][1]) / S;
    } 
    else if ((R[0][0] > R[1][1]) && (R[0][0] > R[2][2])) 
    {
        double S = sqrt(1.0 + R[0][0] - R[1][1] - R[2][2]) * 2;
        q.w = (R[2][1] - R[1][2]) / S;
        q.x = 0.25 * S;
        q.y = (R[0][1] + R[1][0]) / S;
        q.z = (R[0][2] + R[2][0]) / S;
    } 
    else if (R[1][1] > R[2][2]) 
    {
        double S = sqrt(1.0 + R[1][1] - R[0][0] - R[2][2]) * 2;
        q.w = (R[0][2] - R[2][0]) / S;
        q.x = (R[0][1] + R[1][0]) / S;
        q.y = 0.25 * S;
        q.z = (R[1][2] + R[2][1]) / S;
    } 
    else 
    {
        double S = sqrt(1.0 + R[2][2] - R[0][0] - R[1][1]) * 2;
        q.w = (R[1][0] - R[0][1]) / S;
        q.x = (R[0][2] + R[2][0]) / S;
        q.y = (R[1][2] + R[2][1]) / S;
        q.z = 0.25 * S;
    }

    return q;
}


F_VEC_2D get_eocfi_rotation_matrix(const Quaternion& q)
{
    float e00 = q.x*q.x - q.y*q.y - q.z*q.z + q.w*q.w;
    float e01 = 2 * (q.x * q.y + q.z * q.w);
    float e02 = 2 * (q.x * q.z - q.y * q.w);

    float e10 = 2 * (q.x * q.y - q.z * q.w);
    float e11 = -q.x*q.x + q.y*q.y  - q.z*q.z + q.w*q.w;
    float e12 = 2 * (q.y * q.z + q.x * q.w);

    float e20 = 2 * (q.x * q.z + q.y * q.w);
    float e21 = 2 * (q.y * q.z - q.x * q.w);
    float e22 = -q.x*q.x - q.y*q.y + q.z*q.z + q.w*q.w;

    F_VEC_2D R = {
        {-e10, -e11, -e12},
        {-e00, -e01, -e02},
        {-e20, -e21, -e22}
    };

    return R;
}


// https://forum.step.esa.int/uploads/default/original/2X/3/315f415284743cf42be03a464121df1ca337a26d.pdf
Quaternion parse_quaternions(const F_VEC_1D& quaternions)
{
    Quaternion raw_q(quaternions);
    F_VEC_2D eocfi_rotation_matrix = get_eocfi_rotation_matrix(raw_q);
    Quaternion rotated_q = get_quaternions_from_rotation_matrix(eocfi_rotation_matrix);

    return Quaternion({rotated_q.w, -rotated_q.x, -rotated_q.y, -rotated_q.z});
}