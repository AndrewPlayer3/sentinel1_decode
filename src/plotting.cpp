/*
By: Andrew Player
Name: main.cpp
Description: Main function with random command-line arguments for testing purposes. 
             See "SAR Space Packet Protocol Data Unit", for more information on the packet specification:
             https://sentinels.copernicus.eu/documents/247904/2142675/Sentinel-1-SAR-Space-Packet-Protocol-Data-Unit.pdf
             For additional information on Level-0 product decoding, see:
             https://sentinel.esa.int/documents/247904/0/Sentinel-1-Level-0-Data-Decoding-Package.pdf/a8742c59-4914-40c4-8309-c77515649f17
*/
#include "plot.hpp"


unordered_map<string, bool> parse_options(
    const unordered_map<string, bool>& options,
          char* args[],
    const int&  arg_index
) {
    unordered_map<string, bool> selections = options;

    int index = arg_index;

    while(args[index] != __null)
    {
        string users_option = string(args[index]);

        if (selections.contains(users_option))
        {
            selections[users_option] = !selections[users_option];
        }
        index += 1;
    }
    return selections;
}


void validate_args(
    const string& command,
    const vector<string>& command_args,
    const vector<string>& arg_types,
          char*   args[],
    const int&    start_index = 2,
    const string& help = ""
) {
    int end_index = command_args.size() + start_index;
    for (int i = start_index; i < end_index; i++)
    {
        int index = i - start_index;
        if (args[i] == __null)
        {
            cout << command  << " is missing the following argument: "
                 << command_args[index] << endl;
            exit(1);
        }
        try
        {
            if      (arg_types[index] == "string") string(args[i]);
            else if (arg_types[index] == "char")   char(args[i][0]);
            else if (arg_types[index] == "int")    stoi(args[i]);
            else if (arg_types[index] == "float")  stof(args[i]);
            else if (arg_types[index] == "path")   open_file(string(args[i]));
        }
        catch(...)
        {
            cout << args[i] << " is not a valid " << command_args[index] << "." << endl;
            exit(1);
        }
    }
}


string parse_scaling_mode(unordered_map<string, bool> options)
{
    string mode;

    if (options["--norm_log"]) return "norm_log";
    if (options["--norm"])     return "norm";
    if (options["--mag"])      return "mag";
    if (options["--real"])     return "real";
    if (options["--imag"])     return "imag";

    else return "real";
}


void pulse_command(char* argv[], unordered_map<string, bool>& options)
{
    vector<string> args  = {"packet_index", "filepath"};
    vector<string> types = {"int", "path"};
    validate_args("pulse", args, types, argv);

    int packet_index = stoi(argv[2]);
    string filename  = string(argv[3]);
    string scaling   = parse_scaling_mode(options);

    plot_pulse(filename, packet_index, scaling);
}


void pulse_compression_command(char* argv[], unordered_map<string, bool>& options)
{
    vector<string> args = {"packet_index", "filepath"};
    vector<string> types = {"int", "path"};
    validate_args("pulse_compression", args, types, argv);

    int packet_index = stoi(argv[2]);
    string filename  = string(argv[3]);
    string scaling   = parse_scaling_mode(options);
    bool   do_fft    = options["--fft"];

    plot_pulse_compression(filename, packet_index, do_fft, scaling);
}


void pulse_image_command(char* argv[], unordered_map<string, bool>& options)
{
    vector<string> args = {"swath", "burst_num", "filepath"};
    vector<string> types = {"string", "int", "path"};
    validate_args("pulse_img", args, types, argv);

    string swath     = string(argv[2]);
    int    burst_num = stoi(argv[3]);
    string filename  = string(argv[4]);
    string scaling   = parse_scaling_mode(options);

    plot_pulse_image(filename, swath, burst_num, scaling);
}


void range_compressed_burst_command(char* argv[], unordered_map<string, bool>& options)
{
    vector<string> args = {"swath", "burst_num", "filepath"};
    vector<string> types = {"string", "int", "path"};
    validate_args("range_compressed_burst", args, types, argv);

    string swath     = string(argv[2]);
    int    burst_num = stoi(argv[3]);
    string filename  = string(argv[4]);
    string scaling   = parse_scaling_mode(options);

    plot_range_compressed_burst(filename, swath, burst_num, scaling);
}


void range_compressed_swath_command(char* argv[], unordered_map<string, bool>& options)
{
    vector<string> args = {"swath", "filepath"};
    vector<string> types = {"string", "path"};
    validate_args("range_compressed_swath", args, types, argv);

    string swath     = string(argv[2]);
    string filename  = string(argv[3]);
    string scaling   = parse_scaling_mode(options);

    plot_range_compressed_swath(filename, swath, scaling);
}


void fft_axis_command(char *argv[], unordered_map<string, bool>& options)
{
    vector<string> args = {"swath", "burst_num", "axis", "fft_size", "filepath"};
    vector<string> types = {"string", "int", "int", "int", "path"};
    validate_args("fft_axis", args, types, argv);

    string swath    = string(argv[2]);
    int    burst_num = stoi(argv[3]);
    int    axis     = stoi(argv[4]);
    int    fft_size = stoi(argv[5]);
    string filepath = string(argv[6]);
    bool   inverse  = options["--inverse"];
    string scaling  = parse_scaling_mode(options);

    plot_fft_axis(filepath, swath, burst_num, axis, fft_size, inverse, scaling);
}


void fft2_command(char *argv[], unordered_map<string, bool>& options)
{
    vector<string> args = {"swath", "burst_num", "fft_rows", "fft_cols", "filepath"};
    vector<string> types = {"string", "int", "int", "int", "path"};
    validate_args("fft2", args, types, argv);

    string swath    = string(argv[2]);
    int    burst_num = stoi(argv[3]);
    int    fft_rows = stoi(argv[4]);
    int    fft_cols = stoi(argv[5]);
    string filepath = string(argv[6]);
    bool   inverse  = options["--inverse"];
    string scaling  = parse_scaling_mode(options);

    plot_fft2d(filepath, swath, burst_num, fft_rows, fft_cols, inverse, scaling);
}


void fft_command(char *argv[], unordered_map<string, bool>& options)
{
    vector<string> args = {"packet_index", "fft_size", "filepath"};
    vector<string> types = {"int", "int", "path"};
    validate_args("fft", args, types, argv);

    string filename     = string(argv[4]);
    int    packet_index = stoi(argv[2]);
    int    fft_size     = stoi(argv[3]);
    bool   inverse      = options["--inverse"];
    string scaling      = parse_scaling_mode(options);

    plot_fft(filename, packet_index, fft_size, inverse, scaling);
}


void burst_command(char *argv[], unordered_map<string, bool>& options)
{
    vector<string> args = {"swath", "burst_num", "filepath"};
    vector<string> types = {"string", "int", "path"};
    validate_args("burst", args, types, argv);

    string filepath  = string(argv[4]);
    string swath     = string(argv[2]);
    int    burst_num = stoi(argv[3]); 
    string scaling   = parse_scaling_mode(options);

    plot_burst(filepath, swath, burst_num, scaling);
}


void swath_command(char *argv[], unordered_map<string, bool>& options)
{
    vector<string> args = {"swath", "filepath"};
    vector<string> types = {"string", "int", "path"};
    validate_args("swath", args, types, argv);

    string filepath = string(argv[3]);
    string swath    = string(argv[2]);
    string scaling  = parse_scaling_mode(options);

    plot_swath(filepath, swath, scaling);
}


void signal_command(char *argv[], unordered_map<string, bool>& options)
{
    vector<string> args = {"packet_index", "filepath"};
    vector<string> types = {"int", "path"};
    validate_args("signal", args, types, argv);

    string filename     = string(argv[3]);
    int    packet_index = stoi(argv[2]); 
    string scaling      = parse_scaling_mode(options);

    plot_signal(filename, packet_index, scaling);
}


void print_help(vector<string> help_strings)
{
    for (string help_string : help_strings)
    {
        cout << help_string << endl;
    }
}


int main(int argc, char* argv[]) 
{
    vector<string> help_strings = {
        "signal [packet_index] [mode] [path]",
        "swath [swath] [path]",
        "burst [swath] [burst_num] [path]",
        "fft [packet_index] [fft_size] [path] [--inverse]",
        "fft2 [swath] [burst_num] [path] [fft_rows] [fft_cols] [--inverse]",
        "fft_axis [swath] [burst_num] [axis] [fft_size] [path] [--inverse]",
        "range_compressed_burst [swath] [burst_num] [path]"
        "range_compressed_swath [swath] [path]"
        "Scaling Options: [--norm_log|--norm|--mag|--real|--imag]"
    };

    if(argv[1] == __null) 
    {
        cout << "Please enter a command:" << endl;
        print_help(help_strings);
        return 1;
    }

    fftwf_init_threads();
    fftwf_plan_with_nthreads(omp_get_max_threads());

    string command = string(argv[1]);

    unordered_map<string, bool> options = {
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

    else if (command == "help" or command == "--help" or command == "-h")
    {
        print_help(help_strings);
    }
    else
    {
        cout << command << " is not a valid command." << endl;
    }

    fftwf_cleanup_threads();

    return 0;
}
