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
            quaternions.push_back(state_vector.quaternions);
            angular_rates.push_back(state_vector.angular_rate);
            attitudes.push_back(state_vector.attitude);
        }
    }
}


void STATE_VECTORS::print()
{
    int hyphen_count = 210;
    std::cout << std::string(hyphen_count, '-') << std::endl;
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
              << std::setw(15) << "q3_quaternion"
              << std::setw(15) << "roll"
              << std::setw(15) << "pitch"
              << std::setw(15) << "yaw" << std::endl;
    std::cout << std::string(hyphen_count, '-') << std::endl;
    for (size_t i = 0; i < positions.size(); ++i) {
        std::cout   << std::fixed    << std::setprecision(4)
                    << std::setw(15) << times[i]
                    << std::setw(15) << positions[i][0]
                    << std::setw(15) << positions[i][1]
                    << std::setw(15) << positions[i][2]
                    << std::setw(15) << velocities[i][0]
                    << std::setw(15) << velocities[i][1]
                    << std::setw(15) << velocities[i][2]
                    << std::setw(15) << quaternions[i][0]
                    << std::setw(15) << quaternions[i][1]
                    << std::setw(15) << quaternions[i][2]
                    << std::setw(15) << quaternions[i][3]
                    << std::setw(15) << attitudes[i][0]
                    << std::setw(15) << attitudes[i][1]
                    << std::setw(15) << attitudes[i][2] << std::endl;
    }
    std::cout << std::string(hyphen_count, '-') << std::endl;
}


// Find the indexes a and b in times that bound time: I.e. times[a] <= time <= times[b]
std::pair<D_VEC_1D::iterator, D_VEC_1D::iterator> find_time_bounds(const double& time, D_VEC_1D::iterator times_start, D_VEC_1D::iterator times_end)
{
    int size = times_end - times_start;

    if (size == 1) return {times_start, times_end};

    int pivot = size / 2;
    double time_at_pivot = *(times_start + pivot);

    if (time_at_pivot == time)     return {times_start, times_start};
    else if (time_at_pivot > time) return find_time_bounds(time, times_start, times_start+pivot);
    else                           return find_time_bounds(time, times_start+pivot, times_end);
}


D_VEC_1D interpolate_vector(const F_VEC_2D& vec, int a_index, int b_index, double linear_extrapolant)
{
    D_VEC_1D new_vec(vec[a_index].size());

    std::transform(
        vec[a_index].begin(), vec[a_index].end(), 
            vec[b_index].begin(), new_vec.begin(),
                [linear_extrapolant] (double a, double b) { return std::lerp(a, b, linear_extrapolant); }
    );

    return new_vec;
}


STATE_VECTOR STATE_VECTORS::interpolate(const double& time)
{    
    STATE_VECTOR state_vector;

    int num_times = times.size();
    int a_index, b_index;
    double a, b, t;

    if (time > times.back())
    {
        state_vector.velocity = velocities.back();
        state_vector.position = positions.back();
        state_vector.quaternions = quaternions.back();
        state_vector.angular_rate = angular_rates.back();
        state_vector.attitude = attitudes.back();

        return state_vector;
    }
    if (time < times[0])
    {
        state_vector.velocity = velocities[0];
        state_vector.position = positions[0];
        state_vector.quaternions = quaternions[0];
        state_vector.angular_rate = angular_rates[0];
        state_vector.attitude = attitudes[0];

        return state_vector;
    }

    std::pair<D_VEC_1D::iterator, D_VEC_1D::iterator> time_bounds = find_time_bounds(time, times.begin(), times.end());

    a_index = time_bounds.first - times.begin();
    b_index = time_bounds.second - times.begin();

    a = *time_bounds.first;
    b = *time_bounds.second;
    t = (time - a) / (b - a);

    state_vector.velocity = interpolate_vector(velocities, a_index, b_index, t);
    state_vector.position = interpolate_vector(positions, a_index, b_index, t);
    state_vector.quaternions = interpolate_vector(quaternions, a_index, b_index, t);
    state_vector.angular_rate = interpolate_vector(angular_rates, a_index, b_index, t);
    state_vector.attitude = interpolate_vector(attitudes, a_index, b_index, t);

    return state_vector;
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


F_VEC_1D Quaternion::to_euler_angles()
{
    double roll, pitch, yaw;

    F_VEC_1D R = to_rotation_matrix();

    pitch = asin(-std::max(-1.0, std::min(1.0, R[7])));

    if (std::abs(R[7]) < 1)
    {
        roll = atan2(R[6], R[8]);
        yaw = atan2(R[1], R[4]);
    }
    else
    {
        roll = atan2(R[5], R[4]);
        yaw = 0.0;
    }

    return {-pitch, -roll, yaw};
}


D_VEC_1D Quaternion::to_rotation_matrix()
{
    D_VEC_1D R(9);

    R[0] = 1 - 2 * (y*y + z*z);
    R[1] = 2 * (x*y + w*z);
    R[2] = 2 * (x*z - w*y);

    R[3] = 2 * (x*y - w*z);
    R[4] = 1 - 2 * (x*x + z*z);
    R[5] = 2 * (y*z + w*x);

    R[6] = 2 * (x*z + w*y);
    R[7] = 2 * (y*z - w*x);
    R[8] = 1 - 2 * (x*x + y*y);

    return R;
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