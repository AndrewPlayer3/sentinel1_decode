#pragma once

#include <unordered_map>
#include <string>
#include <iostream>

#include "decoding_utils.h"
#include "misc_types.h"


std::unordered_map<std::string, bool> parse_options(
    const std::unordered_map<std::string, bool>& options,
          char* args[],
    const int&  arg_index
);

void validate_args(
    const std::string&   command,
    const STRING_VEC_1D& command_args,
    const STRING_VEC_1D& arg_types,
          char*          args[],
    const int&           start_index = 2
);

std::string parse_scaling_mode(
    std::unordered_map<std::string, bool> options
);

void print_help(
    const STRING_VEC_1D& help_strings
);