#include "tiff.h"
#include "tiffio.h"
#include "signal_processing.h"
#include "image_formation.h"
#include "image_write.h"
#include "cli.h"
#include "packet.h"
#include "burst.h"
#include "swath.h"


void burst_command(char *argv[], std::unordered_map<std::string, bool>& options)
{
    STRING_VEC_1D args  = {"swath", "burst_num", "in_path", "out_path"};
    STRING_VEC_1D types = {"string", "int", "path", "string"};
    validate_args("burst", args, types, argv);

    std::string swath    = std::string(argv[2]);
    int burst_num        = std::stoi(argv[3]);
    std::string in_path  = std::string(argv[4]);
    std::string out_path = std::string(argv[5]);
    std::string scaling  = parse_scaling_mode(options);

    write_burst(in_path, out_path, swath, burst_num, scaling);
}


void swath_command(char *argv[], std::unordered_map<std::string, bool>& options)
{
    STRING_VEC_1D args  = {"swath", "in_path", "out_path"};
    STRING_VEC_1D types = {"string", "path", "string"};
    validate_args("swath", args, types, argv);

    std::string swath    = std::string(argv[2]);
    std::string in_path  = std::string(argv[3]);
    std::string out_path = std::string(argv[4]);
    std::string scaling  = parse_scaling_mode(options);

    write_swath(in_path, out_path, swath, scaling);
}


void burst_replica_chirps_command(char *argv[], std::unordered_map<std::string, bool>& options)
{
    STRING_VEC_1D args  = {"swath", "burst_num", "in_path", "out_path"};
    STRING_VEC_1D types = {"string", "int", "path", "string"};
    validate_args("burst", args, types, argv);

    std::string swath    = std::string(argv[2]);
    int burst_num        = std::stoi(argv[3]);
    std::string in_path  = std::string(argv[4]);
    std::string out_path = std::string(argv[5]);
    std::string scaling  = parse_scaling_mode(options);

    write_burst_replica_chirps(in_path, out_path, swath, burst_num, scaling);
}


void swath_replica_chirps_command(char *argv[], std::unordered_map<std::string, bool>& options)
{
    STRING_VEC_1D args  = {"swath", "in_path", "out_path"};
    STRING_VEC_1D types = {"string", "path", "string"};
    validate_args("swath", args, types, argv);

    std::string swath    = std::string(argv[2]);
    std::string in_path  = std::string(argv[3]);
    std::string out_path = std::string(argv[4]);
    std::string scaling  = parse_scaling_mode(options);

    write_swath_replica_chirps(in_path, out_path, swath, scaling);
}


void range_compressed_burst_command(char *argv[], std::unordered_map<std::string, bool>& options)
{
    STRING_VEC_1D args  = {"swath", "burst_num", "in_path", "out_path"};
    STRING_VEC_1D types = {"string", "int", "path", "string"};
    validate_args("range_compressed_burst", args, types, argv);

    std::string swath    = std::string(argv[2]);
    int burst_num        = std::stoi(argv[3]);
    std::string in_path  = std::string(argv[4]);
    std::string out_path = std::string(argv[5]);
    std::string scaling  = parse_scaling_mode(options);

    write_range_compressed_burst(in_path, out_path, swath, burst_num, scaling);
}


void range_compressed_swath_command(char *argv[], std::unordered_map<std::string, bool>& options)
{
    STRING_VEC_1D args  = {"swath", "in_path", "out_path"};
    STRING_VEC_1D types = {"string", "path", "string"};
    validate_args("range_compressed_swath", args, types, argv);

    std::string swath    = std::string(argv[2]);
    std::string in_path  = std::string(argv[3]);
    std::string out_path = std::string(argv[4]);
    std::string scaling  = parse_scaling_mode(options);

    write_range_compressed_swath(in_path, out_path, swath, scaling);
}


// void range_doppler_burst_command(char *argv[], std::unordered_map<std::string, bool>& options)
// {
//     STRING_VEC_1D args  = {"swath", "burst_num", "in_path", "out_path"};
//     STRING_VEC_1D types = {"string", "int", "path", "string"};
//     validate_args("range_doppler_burst", args, types, argv);

//     std::string swath    = std::string(argv[2]);
//     int burst_num        = std::stoi(argv[3]);
//     std::string in_path  = std::string(argv[4]);
//     std::string out_path = std::string(argv[5]);
//     std::string scaling  = parse_scaling_mode(options);

//     write_range_doppler_burst(in_path, out_path, swath, burst_num, scaling);
// }


// void range_doppler_swath_command(char *argv[], std::unordered_map<std::string, bool>& options)
// {
//     STRING_VEC_1D args  = {"swath", "in_path", "out_path"};
//     STRING_VEC_1D types = {"string", "path", "string"};
//     validate_args("range_doppler_swath", args, types, argv);

//     std::string swath    = std::string(argv[2]);
//     std::string in_path  = std::string(argv[3]);
//     std::string out_path = std::string(argv[4]);
//     std::string scaling  = parse_scaling_mode(options);

//     write_range_doppler_swath(in_path, out_path, swath, scaling);
// }


void azimuth_compressed_burst_command(char *argv[], std::unordered_map<std::string, bool>& options)
{
    STRING_VEC_1D args  = {"swath", "burst_num", "in_path", "out_path"};
    STRING_VEC_1D types = {"string", "int", "path", "string"};
    validate_args("azimuth_compressed_burst", args, types, argv);

    std::string swath    = std::string(argv[2]);
    int burst_num        = std::stoi(argv[3]);
    std::string in_path  = std::string(argv[4]);
    std::string out_path = std::string(argv[5]);
    std::string scaling  = parse_scaling_mode(options);

    write_azimuth_compressed_burst(in_path, out_path, swath, burst_num, scaling);
}


void azimuth_compressed_swath_command(char *argv[], std::unordered_map<std::string, bool>& options)
{
    STRING_VEC_1D args  = {"swath", "in_path", "out_path"};
    STRING_VEC_1D types = {"string", "path", "string"};
    validate_args("azimuth_compressed_swath", args, types, argv);

    std::string swath    = std::string(argv[2]);
    std::string in_path  = std::string(argv[3]);
    std::string out_path = std::string(argv[4]);
    std::string scaling  = parse_scaling_mode(options);

    write_azimuth_compressed_swath(in_path, out_path, swath, scaling);
}


int main(int argc, char* argv[])
{
    STRING_VEC_1D help_strings = {
        "burst [swath] [burst_num] [in_path] [out_path]",
        "swath [swath] [in_path] [out_path]",
        "burst_replica_chirps [swath] [burst_num] [in_path] [out_path]",
        "swath_replica_chirps [swath] [in_path] [out_path]",
        "range_compressed_burst [swath] [burst_num] [in_path] [out_path]",
        "range_compressed_swath [swath] [in_path] [out_path]",
        "range_doppler_burst [swath] [burst_num] [in_path] [out_path]",
        "range_doppler_swath [swath] [in_path] [out_path]",
        "azimuth_compressed_burst [swath] [burst_num] [in_path] [out_path]",
        "azimuth_compressed_swath [swath] [in_path] [out_path]",
        "Scaling Options: [--norm_log|--norm|--mag|--real|--imag]"
    };

    if(argv[1] == __null) 
    {
        std::cout << "Please enter a command:" << std::endl;
        print_help(help_strings);
        return 1;
    }

    std::string command = std::string(argv[1]);

    std::unordered_map<std::string, bool> options = {
        {"--inverse",  false},
        {"--fft",      false},
        {"--norm",     false},
        {"--norm_log", false},
        {"--mag",      false},
        {"--real",     false},
        {"--imag",     false}
    };
    options = parse_options(options, argv, 2);

    if      (command == "burst")                    burst_command(&(argv[0]), options);
    else if (command == "swath")                    swath_command(&(argv[0]), options);
    else if (command == "burst_replica_chirps")     burst_replica_chirps_command(&(argv[0]), options);
    else if (command == "swath_replica_chirps")     swath_replica_chirps_command(&(argv[0]), options);
    else if (command == "range_compressed_burst")   range_compressed_burst_command(&(argv[0]), options);
    else if (command == "range_compressed_swath")   range_compressed_swath_command(&(argv[0]), options);
    // else if (command == "range_doppler_burst")      range_doppler_burst_command(&(argv[0]), options);
    // else if (command == "range_doppler_swath")      range_doppler_swath_command(&(argv[0]), options);
    else if (command == "azimuth_compressed_burst") azimuth_compressed_burst_command(&(argv[0]), options);
    else if (command == "azimuth_compressed_swath") azimuth_compressed_swath_command(&(argv[0]), options);
    else if (command == "help" or command == "--help" or command == "-h")
    {
        print_help(help_strings);
    }
    else
    {
        std::cout << command << " is not a valid command." << std::endl;
    }

    return 0;
}
