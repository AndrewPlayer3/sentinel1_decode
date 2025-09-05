#pragma once

#include <algorithm>
#include <cmath>
#include "packet.h"

typedef struct std::unordered_map<std::string, u_int64_t> SUBCOMM_DICT_INT;
typedef struct std::unordered_map<std::string, double> SUBCOMM_DICT_DOUBLE;
typedef struct std::vector<SUBCOMM_DICT_DOUBLE> SUBCOMM_DICTS;

SUBCOMM_DICTS build_data_word_dicts(PACKET_VEC_1D& packets);

enum QUAT_ORDER {
    SCALAR_FIRST,
    SCALAR_LAST
};

// Scalar-first quaternion i.e. w, x, y, z
struct Quaternion
{
    double w, x, y, z;

    Quaternion() {}

    Quaternion(D_VEC_1D quaternions, QUAT_ORDER order = SCALAR_FIRST)
    {
        if (order == SCALAR_FIRST)
        {
            w = quaternions[0];
            x = quaternions[1];
            y = quaternions[2];
            z = quaternions[3];
        }
        else
        {
            w = quaternions[3];
            x = quaternions[0];
            y = quaternions[1];
            z = quaternions[2];
        }
    }

    D_VEC_1D to_vector(QUAT_ORDER order = SCALAR_FIRST)
    {
        if (order == SCALAR_FIRST) return {w, x, y, z};
        else                       return {x, y, z, w};
    }

    D_VEC_1D to_euler_angles();

    D_VEC_1D to_rotation_matrix();
};


Quaternion get_quaternions_from_rotation_matrix(F_VEC_2D& R);
F_VEC_2D get_eocfi_rotation_matrix(const Quaternion& q);
Quaternion parse_quaternions(const D_VEC_1D& quaternions);

struct STATE_VECTOR
{
    double time;

    D_VEC_1D position;
    D_VEC_1D velocity;
    D_VEC_1D quaternions;
    D_VEC_1D angular_rate;
    D_VEC_1D attitude;

    STATE_VECTOR() {}

    STATE_VECTOR(SUBCOMM_DICT_DOUBLE& subcomm_dict) 
    {
        time = subcomm_dict.at("data_time_stamp");
 
        position = {
            subcomm_dict.at("x_axis_position"),
            subcomm_dict.at("y_axis_position"),
            subcomm_dict.at("z_axis_position"),
        };

        velocity = {
            subcomm_dict.at("x_axis_velocity"),
            subcomm_dict.at("y_axis_velocity"),
            subcomm_dict.at("z_axis_velocity"),
        };

        angular_rate = {
            subcomm_dict.at("omega_x"),
            subcomm_dict.at("omega_y"),
            subcomm_dict.at("omega_z"),
        };

        Quaternion quat = parse_quaternions({
            subcomm_dict.at("q0_quaternion"),
            subcomm_dict.at("q1_quaternion"),
            subcomm_dict.at("q2_quaternion"),
            subcomm_dict.at("q3_quaternion"),
        });

        quaternions = quat.to_vector(SCALAR_LAST);

        attitude = quat.to_euler_angles();
    }
};

struct STATE_VECTORS
{
    D_VEC_1D times;

    D_VEC_2D positions;
    D_VEC_2D velocities;
    D_VEC_2D quaternions;
    D_VEC_2D angular_rates;
    D_VEC_2D attitudes;

    STATE_VECTORS() {}

    STATE_VECTORS(SUBCOMM_DICTS& subcomm_dicts) 
    {
        parse_subcomm_dicts(subcomm_dicts);
    }

    STATE_VECTORS(PACKET_VEC_1D& packets) 
    {
        SUBCOMM_DICTS subcomm_dicts = build_data_word_dicts(packets);
        parse_subcomm_dicts(subcomm_dicts);
    }

    void parse_subcomm_dicts(SUBCOMM_DICTS subcomm_dicts);
    void print();

    STATE_VECTOR interpolate(const double& time);
};
