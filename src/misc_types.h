/*
By: Andrew Player
Name: structs.hpp
Description: Some typedefs/structs for convinience.
*/

#pragma once

#include <cstdint>
#include <complex>
#include <string>
#include <vector>
#include <unordered_map>

typedef std::uint8_t u_int8_t;
typedef std::uint16_t u_int16_t;
typedef std::uint32_t u_int32_t;
typedef std::uint64_t u_int64_t;

using SUBCOMM_DICT_INT    = std::unordered_map<std::string, std::uint64_t>;
using SUBCOMM_DICT_DOUBLE = std::unordered_map<std::string, double>;
using SUBCOMM_DICTS       = std::vector<SUBCOMM_DICT_DOUBLE>;

template <typename T>
using VEC_1D = std::vector<T>;
template <typename T>
using VEC_2D = std::vector<std::vector<T>>;

using F_VEC_1D   = VEC_1D<float>;
using F_VEC_2D   = VEC_2D<float>;
using D_VEC_1D   = VEC_1D<double>;
using D_VEC_2D   = VEC_2D<double>;
using INT_VEC_1D = VEC_1D<int>;
using INT_VEC_2D = VEC_2D<int>;

using UINT8_VEC_1D  = VEC_1D<u_int8_t>;
using UINT8_VEC_2D  = VEC_2D<u_int8_t>;
using UINT16_VEC_1D = VEC_1D<u_int16_t>;
using UINT16_VEC_2D = VEC_2D<u_int16_t>;
using UINT32_VEC_1D = VEC_1D<u_int32_t>;
using UINT32_VEC_2D = VEC_2D<u_int32_t>;
using UINT64_VEC_1D = VEC_1D<u_int64_t>;
using UINT64_VEC_2D = VEC_2D<u_int64_t>;

using STRING_VEC_1D = VEC_1D<std::string>;
using STRING_VEC_2D = VEC_2D<std::string>;

using cdouble    = std::complex<double>;
using CF_VEC_1D  = VEC_1D<cdouble>;
using CF_VEC_2D  = VEC_2D<cdouble>;

struct SIGNAL_PAIR {
    CF_VEC_2D signals;
    CF_VEC_2D replica_chirps;
};
