/*
By: Andrew Player
Name: structs.hpp
Description: Some typedefs/structs for convinience.
*/

#pragma once

#include <complex>
#include <vector>


typedef struct std::vector<double>                    F_VEC_1D;
typedef struct std::vector<std::vector<double>>       F_VEC_2D;

typedef struct std::vector<double>                   D_VEC_1D;
typedef struct std::vector<std::vector<double>>      D_VEC_2D;

typedef struct std::vector<int>                      INT_VEC_1D;
typedef struct std::vector<std::vector<int>>         INT_VEC_2D;

typedef struct std::vector<u_int8_t>                 UINT8_VEC_1D;
typedef struct std::vector<std::vector<u_int8_t>>    UINT8_VEC_2D;

typedef struct std::vector<u_int16_t>                UINT16_VEC_1D;
typedef struct std::vector<std::vector<u_int16_t>>   UINT16_VEC_2D;

typedef struct std::vector<u_int32_t>                UINT32_VEC_1D;
typedef struct std::vector<std::vector<u_int32_t>>   UINT32_VEC_2D;

typedef struct std::vector<u_int64_t>                UINT64_VEC_1D;
typedef struct std::vector<std::vector<u_int64_t>>   UINT64_VEC_2D;

typedef struct std::vector<std::string>              STRING_VEC_1D;
typedef struct std::vector<std::vector<std::string>> STRING_VEC_2D;

typedef struct std::vector<std::complex<double>>              CF_VEC_1D;
typedef struct std::vector<std::vector<std::complex<double>>> CF_VEC_2D;


struct SIGNAL_PAIR {
    CF_VEC_2D signals;
    CF_VEC_2D replica_chirps;
};