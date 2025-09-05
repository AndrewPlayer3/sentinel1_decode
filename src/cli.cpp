#include "cli.h"

std::unordered_map<std::string, bool> parse_options(
    const std::unordered_map<std::string, bool>& options,
          char* args[],
    const int&  arg_index
) {
    std::unordered_map<std::string, bool> selections = options;

    int index = arg_index;

    while(args[index] != __null)
    {
        std::string users_option = std::string(args[index]);

        if (selections.contains(users_option))
        {
            selections[users_option] = !selections[users_option];
        }
        index += 1;
    }
    return selections;
}


void validate_args(
    const std::string&   command,
    const STRING_VEC_1D& command_args,
    const STRING_VEC_1D& arg_types,
          char*          args[],
    const int&           start_index
) {
    int end_index = command_args.size() + start_index;
    for (int i = start_index; i < end_index; i++)
    {
        int index = i - start_index;
        if (args[i] == NULL)
        {
            std::cout << command  << " is missing the following argument: "
                 << command_args[index] << std::endl;
            exit(1);
        }
        try
        {
            if      (arg_types[index] == "string") (std::string(args[i]));
            else if (arg_types[index] == "char")   char test = (char(args[i][0]));
            else if (arg_types[index] == "int")    std::stoi(args[i]);
            else if (arg_types[index] == "double")  std::stof(args[i]);
            else if (arg_types[index] == "path")   open_file(std::string(args[i]));
        }
        catch(...)
        {
            std::cout << args[i] << " is not a valid " << command_args[index];

            if (arg_types[index] == "path")
            {
                std::cout << " or the file may not exist";
            }

            std::cout << std::endl;

            exit(1);
        }
    }
}


std::string parse_scaling_mode(
    std::unordered_map<std::string, bool> options
) {
    std::string mode;

    if (options["--norm_log"]) return "norm_log";
    if (options["--norm"])     return "norm";
    if (options["--mag"])      return "mag";
    if (options["--real"])     return "real";
    if (options["--imag"])     return "imag";

    else return "real";
}


void print_help(const STRING_VEC_1D& help_strings)
{
    for (std::string help_string : help_strings)
    {
        std::cout << help_string << std::endl;
    }
}
