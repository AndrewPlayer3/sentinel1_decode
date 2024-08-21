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
          char*   args[],
    const int&    start_index = 2,
    const string& help = ""
) {
    string error_string = "Error! " + command + " is missing the following argument: ";

    int end_index = command_args.size() + start_index;
    for (int i = start_index; i < end_index; i++)
    {
        if (args[i] == __null)
        {
            throw runtime_error(error_string + command_args[i]);
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


void fft_axis_command(char *argv[], unordered_map<string, bool>& options)
{
    vector<string> args = {"swath", "axis", "fft_size", "filepath"};

    validate_args("fft_axis", args, argv);

    string swath    = string(argv[2]);
    int    axis     = stoi(argv[3]);
    int    fft_size = stoi(argv[4]);
    string filepath = string(argv[5]);
    bool   inverse  = options["--inverse"];
    string scaling  = parse_scaling_mode(options);

    plot_fft_axis(filepath, swath, axis, fft_size, inverse, scaling);
}


void fft2_command(char *argv[], unordered_map<string, bool>& options)
{
    vector<string> args = {"swath", "fft_rows", "fft_cols", "filepath"};

    validate_args("fft2", args, argv);

    string swath    = string(argv[2]);
    int    fft_rows = stoi(argv[3]);
    int    fft_cols = stoi(argv[4]);
    string filepath = string(argv[5]);
    bool   inverse  = options["--inverse"];
    string scaling  = parse_scaling_mode(options);

    plot_fft2d(filepath, swath, fft_rows, fft_cols, inverse, scaling);
}


void fft_command(char *argv[], unordered_map<string, bool>& options)
{
    vector<string> args = {"packet_index", "fft_size", "filepath"};

    validate_args("fft", args, argv);

    string filename     = string(argv[4]);
    int    packet_index = stoi(argv[2]);
    int    fft_size     = stoi(argv[3]);
    bool   inverse      = options["--inverse"];
    string scaling      = parse_scaling_mode(options);

    plot_fft(filename, packet_index, fft_size, inverse, scaling);
}


void swath_command(char *argv[], unordered_map<string, bool>& options)
{
    vector<string> args = {"swath", "filepath"};

    validate_args("swath", args, argv, 2);

    string filepath = string(argv[3]);
    string swath    = string(argv[2]);
    string scaling  = parse_scaling_mode(options);

    plot_swath(filepath, swath, scaling);
}


void complex_samples_command(char *argv[], unordered_map<string, bool>& options)
{
    vector<string> args = {"packet_index", "filepath"};

    validate_args("complex_samples", args, argv);

    string filename     = string(argv[3]);
    int    packet_index = stoi(argv[2]); 
    string scaling      = parse_scaling_mode(options);

    plot_complex_samples(filename, packet_index, scaling);
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
        "complex_samples [packet_index] [mode] [path]",
        "swath [swath] [path]",
        "fft [packet_index] [fft_size] [path] [--inverse]",
        "fft2 [swath] [path] [fft_rows] [fft_cols] [--inverse]",
        "fft_axis [swath] [axis] [fft_size] [path] [--inverse]",
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
        {"--norm",     false},
        {"--norm_log", false},
        {"--mag",      false},
        {"--real",     false},
        {"--imag",     false}
    };
    options = parse_options(options, argv, 2);

    if      (command == "complex_samples") complex_samples_command(&(argv[0]), options);
    else if (command == "swath")           swath_command(&(argv[0]), options);
    else if (command == "fft")             fft_command(&(argv[0]), options); 
    else if (command == "fft2")            fft2_command(&(argv[0]), options);
    else if (command == "fft_axis")        fft_axis_command(&(argv[0]), options);

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
