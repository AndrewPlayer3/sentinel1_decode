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


void range_doppler_burst_command(char *argv[], std::unordered_map<std::string, bool>& options)
{
    STRING_VEC_1D args  = {"swath", "burst_num", "in_path", "out_path"};
    STRING_VEC_1D types = {"string", "int", "path", "string"};
    validate_args("range_compressed_burst", args, types, argv);

    std::string swath    = std::string(argv[2]);
    int burst_num        = std::stoi(argv[3]);
    std::string in_path  = std::string(argv[4]);
    std::string out_path = std::string(argv[5]);
    std::string scaling  = parse_scaling_mode(options);

    write_range_doppler_burst(in_path, out_path, swath, burst_num, scaling);
}


void range_doppler_swath_command(char *argv[], std::unordered_map<std::string, bool>& options)
{
    STRING_VEC_1D args  = {"swath", "in_path", "out_path"};
    STRING_VEC_1D types = {"string", "path", "string"};
    validate_args("range_compressed_swath", args, types, argv);

    std::string swath    = std::string(argv[2]);
    std::string in_path  = std::string(argv[3]);
    std::string out_path = std::string(argv[4]);
    std::string scaling  = parse_scaling_mode(options);

    write_range_doppler_swath(in_path, out_path, swath, scaling);
}


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


void save_raw_npy_command(char *argv[], std::unordered_map<std::string, bool>& options)
{
    STRING_VEC_1D args  = {"swath", "in_path", "out_path"};
    STRING_VEC_1D types = {"string", "path", "string"};
    validate_args("azimuth_compressed_swath", args, types, argv);

    std::string swath    = std::string(argv[2]);
    std::string in_path  = std::string(argv[3]);
    std::string out_path = std::string(argv[4]);

    // Open the file in binary mode
    std::ofstream file(out_path, std::ios::binary);
    if (!file) {
        std::cerr << "Error: Unable to open file " << out_path << " for writing." << std::endl;
        return;
    }

    S1_Decoder s1(in_path);
    CF_VEC_2D data = s1.get_swath(swath);

    int rows = data.size();
    int cols = data[0].size();

    std::cout << "Num Rows: " << rows << std::endl;
    std::cout << "Num Cols: " << cols << std::endl;

    // Iterate through the 2D vector
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            std::complex<double> value = data[i][j];
            // Convert each complex<double> to complex<float>
            std::complex<float> valueFloat(static_cast<float>(value.real()), static_cast<float>(value.imag()));
            
            // Write the real and imaginary parts as interleaved floats
            file.write(reinterpret_cast<const char*>(&valueFloat), sizeof(std::complex<float>));
        }
    }

    file.close();

    if (!file) {
        std::cerr << "Error: Failed to write to file " << out_path << "." << std::endl;
    } else {
        std::cout << "Data successfully saved to " << out_path << "." << std::endl;
    }
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
    else if (command == "range_doppler_burst")      range_doppler_burst_command(&(argv[0]), options);
    else if (command == "range_doppler_swath")      range_doppler_swath_command(&(argv[0]), options);
    else if (command == "azimuth_compressed_burst") azimuth_compressed_burst_command(&(argv[0]), options);
    else if (command == "azimuth_compressed_swath") azimuth_compressed_swath_command(&(argv[0]), options);
    else if (command == "save_raw_npy")             save_raw_npy_command(&(argv[0]), options);
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
