/*
By: Andrew Player
Name: main.cpp
Description: Main function with random command-line arguments for testing purposes. 
             See "SAR Space Packet Protocol Data Unit", for more information on the packet specification:
             https://sentinels.copernicus.eu/documents/247904/2142675/Sentinel-1-SAR-Space-Packet-Protocol-Data-Unit.pdf
             For additional information on Level-0 product decoding, see:
             https://sentinel.esa.int/documents/247904/0/Sentinel-1-Level-0-Data-Decoding-Package.pdf/a8742c59-4914-40c4-8309-c77515649f17
*/
#include "plot.h"


void pulse_command(char* argv[], std::unordered_map<std::string, bool>& options)
{
    STRING_VEC_1D args  = {"packet_index", "filepath"};
    STRING_VEC_1D types = {"int", "path"};
    validate_args("pulse", args, types, argv);

    int packet_index     = std::stoi(argv[2]);
    std::string filename = std::string(argv[3]);
    std::string scaling  = parse_scaling_mode(options);

    plot_pulse(filename, packet_index, scaling);
}


void pulse_compression_command(char* argv[], std::unordered_map<std::string, bool>& options)
{
    STRING_VEC_1D args  = {"packet_index", "filepath"};
    STRING_VEC_1D types = {"int", "path"};
    validate_args("pulse_compression", args, types, argv);

    int packet_index     = std::stoi(argv[2]);
    std::string filename = std::string(argv[3]);
    std::string scaling  = parse_scaling_mode(options);

    plot_pulse_compression(filename, packet_index, scaling);
}


void pulse_image_command(char* argv[], std::unordered_map<std::string, bool>& options)
{
    STRING_VEC_1D args  = {"swath", "burst_num", "filepath"};
    STRING_VEC_1D types = {"std::string", "int", "path"};
    validate_args("pulse_img", args, types, argv);

    std::string swath    = std::string(argv[2]);
    int burst_num        = std::stoi(argv[3]);
    std::string filename = std::string(argv[4]);
    std::string scaling  = parse_scaling_mode(options);

    plot_pulse_image(filename, swath, burst_num, scaling);
}


void range_compressed_burst_command(char* argv[], std::unordered_map<std::string, bool>& options)
{
    STRING_VEC_1D args = {"swath", "burst_num", "filepath"};
    STRING_VEC_1D types = {"string", "int", "path"};
    validate_args("range_compressed_burst", args, types, argv);

    std::string swath    = std::string(argv[2]);
    int burst_num        = std::stoi(argv[3]);
    std::string filename = std::string(argv[4]);
    std::string scaling  = parse_scaling_mode(options);

    plot_range_compressed_burst(filename, swath, burst_num, scaling);
}


void range_compressed_swath_command(char* argv[], std::unordered_map<std::string, bool>& options)
{
    STRING_VEC_1D args  = {"swath", "filepath"};
    STRING_VEC_1D types = {"string", "path"};
    validate_args("range_compressed_swath", args, types, argv);

    std::string swath     = std::string(argv[2]);
    std::string filename  = std::string(argv[3]);
    std::string scaling   = parse_scaling_mode(options);

    plot_range_compressed_swath(filename, swath, scaling);
}


void range_doppler_swath_command(char* argv[], std::unordered_map<std::string, bool>& options)
{
    STRING_VEC_1D args  = {"swath", "filepath"};
    STRING_VEC_1D types = {"string", "path"};
    validate_args("range_doppler_swath", args, types, argv);

    std::string swath     = std::string(argv[2]);
    std::string filename  = std::string(argv[3]);
    std::string scaling   = parse_scaling_mode(options);

    plot_range_doppler_swath(filename, swath, scaling);
}


void azimuth_compressed_burst_command(char* argv[], std::unordered_map<std::string, bool>& options)
{
    STRING_VEC_1D args = {"swath", "burst_num", "filepath"};
    STRING_VEC_1D types = {"string", "int", "path"};
    validate_args("azimuth_compressed_burst", args, types, argv);

    std::string swath    = std::string(argv[2]);
    int burst_num        = std::stoi(argv[3]);
    std::string filename = std::string(argv[4]);
    std::string scaling  = parse_scaling_mode(options);

    plot_azimuth_compressed_burst(filename, swath, burst_num, scaling);
}


void azimuth_compressed_swath_command(char* argv[], std::unordered_map<std::string, bool>& options)
{
    STRING_VEC_1D args = {"swath", "filepath"};
    STRING_VEC_1D types = {"string", "path"};
    validate_args("azimuth_compressed_swath", args, types, argv);

    std::string swath    = std::string(argv[2]);
    std::string filename = std::string(argv[3]);
    std::string scaling  = parse_scaling_mode(options);

    plot_azimuth_compressed_swath(filename, swath, scaling);
}


void fft_axis_command(char *argv[], std::unordered_map<std::string, bool>& options)
{
    STRING_VEC_1D args  = {"swath", "burst_num", "axis", "fft_size", "filepath"};
    STRING_VEC_1D types = {"string", "int", "int", "int", "path"};
    validate_args("fft_axis", args, types, argv);

    std::string swath = std::string(argv[2]);
    int  burst_num = std::stoi(argv[3]);
    int  axis      = std::stoi(argv[4]);
    int  fft_size  = std::stoi(argv[5]);
    bool inverse   = options["--inverse"];
    std::string filepath = std::string(argv[6]);
    std::string scaling  = parse_scaling_mode(options);

    plot_fft_axis(filepath, swath, burst_num, axis, fft_size, inverse, scaling);
}


void fft2_command(char *argv[], std::unordered_map<std::string, bool>& options)
{
    STRING_VEC_1D args  = {"swath", "burst_num", "fft_rows", "fft_cols", "filepath"};
    STRING_VEC_1D types = {"string", "int", "int", "int", "path"};
    validate_args("fft2", args, types, argv);

    std::string swath = std::string(argv[2]);
    int  burst_num = std::stoi(argv[3]);
    int  fft_rows  = std::stoi(argv[4]);
    int  fft_cols  = std::stoi(argv[5]);
    bool inverse   = options["--inverse"];
    std::string filepath = std::string(argv[6]);
    std::string scaling  = parse_scaling_mode(options);

    plot_fft2d(filepath, swath, burst_num, fft_rows, fft_cols, inverse, scaling);
}


void fft_command(char *argv[], std::unordered_map<std::string, bool>& options)
{
    STRING_VEC_1D args  = {"packet_index", "fft_size", "filepath"};
    STRING_VEC_1D types = {"int", "int", "path"};
    validate_args("fft", args, types, argv);

    std::string filename = std::string(argv[4]);
    int  packet_index = std::stoi(argv[2]);
    int  fft_size     = std::stoi(argv[3]);
    bool inverse      = options["--inverse"];
    std::string scaling = parse_scaling_mode(options);

    plot_fft(filename, packet_index, fft_size, inverse, scaling);
}


void burst_command(char *argv[], std::unordered_map<std::string, bool>& options)
{
    STRING_VEC_1D args  = {"swath", "burst_num", "filepath"};
    STRING_VEC_1D types = {"string", "int", "path"};
    validate_args("burst", args, types, argv);

    std::string filepath  = std::string(argv[4]);
    std::string swath     = std::string(argv[2]);
    int burst_num = std::stoi(argv[3]); 
    std::string scaling   = parse_scaling_mode(options);

    plot_burst(filepath, swath, burst_num, scaling);
}


void swath_command(char *argv[], std::unordered_map<std::string, bool>& options)
{
    STRING_VEC_1D args = {"swath", "filepath"};
    STRING_VEC_1D types = {"std::string", "path"};
    validate_args("swath", args, types, argv);

    std::string filepath = std::string(argv[3]);
    std::string swath    = std::string(argv[2]);
    std::string scaling  = parse_scaling_mode(options);

    plot_swath(filepath, swath, scaling);
}


void signal_command(char *argv[], std::unordered_map<std::string, bool>& options)
{
    STRING_VEC_1D args = {"packet_index", "filepath"};
    STRING_VEC_1D types = {"int", "path"};
    validate_args("signal", args, types, argv);

    std::string filename = std::string(argv[3]);
    int packet_index     = std::stoi(argv[2]); 
    std::string scaling  = parse_scaling_mode(options);

    plot_signal(filename, packet_index, scaling);
}


int main(int argc, char* argv[]) 
{
    STRING_VEC_1D help_strings = {
        "signal [packet_index] [mode] [path]",
        "swath [swath] [path]",
        "burst [swath] [burst_num] [path]",
        "fft [packet_index] [fft_size] [path] [--inverse]",
        "fft2 [swath] [burst_num] [path] [fft_rows] [fft_cols] [--inverse]",
        "fft_axis [swath] [burst_num] [axis] [fft_size] [path] [--inverse]",
        "range_compressed_burst [swath] [burst_num] [path]",
        "range_compressed_swath [swath] [path]",
        "range_doppler_swath [swath] [path]"
        "azimuth_compressed_burst [swath] [burst_num] [path]",
        "azimuth_compressed_swath [swath] [path]",
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

    if      (command == "signal")                 signal_command(&(argv[0]), options);
    else if (command == "swath")                  swath_command(&(argv[0]), options);
    else if (command == "burst")                  burst_command(&(argv[0]), options);
    else if (command == "fft")                    fft_command(&(argv[0]), options); 
    else if (command == "fft2")                   fft2_command(&(argv[0]), options);
    else if (command == "fft_axis")               fft_axis_command(&(argv[0]), options);
    else if (command == "pulse")                  pulse_command(&(argv[0]), options);
    else if (command == "pulse_compression")      pulse_compression_command(&(argv[0]), options);
    else if (command == "burst_pulses")           pulse_image_command(&(argv[0]), options);
    else if (command == "range_compressed_burst") range_compressed_burst_command(&(argv[0]), options);
    else if (command == "range_compressed_swath") range_compressed_swath_command(&(argv[0]), options);
    else if (command == "range_doppler_swath")    range_doppler_swath_command(&(argv[0]), options);
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
