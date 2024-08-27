#include "tiff.h"
#include "tiffio.h"
#include "signal_processing.h"
#include "packet.h"
#include "burst.h"
#include "swath.h"


void write_tif(std::vector<float>& img_data, int& rows, int& cols, const std::string out_filename)
{
    TIFF* tif = TIFFOpen(out_filename.c_str(), "w");
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH, rows);
    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, cols);
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 32);
    TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
    TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
    for (int i = 0; i < rows; i++)
    {
        TIFFWriteScanline(tif, &(img_data[i * cols]), i, 0);
    }
    TIFFClose(tif);
}


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
    const int&           start_index = 2,
    const std::string&   help = ""
) {
    int end_index = command_args.size() + start_index;
    for (int i = start_index; i < end_index; i++)
    {
        int index = i - start_index;
        if (args[i] == __null)
        {
            std::cout << command  << " is missing the following argument: "
                 << command_args[index] << std::endl;
            exit(1);
        }
        try
        {
            if      (arg_types[index] == "string") std::string(args[i]);
            else if (arg_types[index] == "char")   char(args[i][0]);
            else if (arg_types[index] == "int")    std::stoi(args[i]);
            else if (arg_types[index] == "float")  std::stof(args[i]);
            else if (arg_types[index] == "path")   open_file(std::string(args[i]));
        }
        catch(...)
        {
            std::cout << args[i] << " is not a valid " << command_args[index] << "." << std::endl;
            exit(1);
        }
    }
}


std::string parse_scaling_mode(std::unordered_map<std::string, bool> options)
{
    std::string mode;

    if (options["--norm_log"]) return "norm_log";
    if (options["--norm"])     return "norm";
    if (options["--mag"])      return "mag";
    if (options["--real"])     return "real";
    if (options["--imag"])     return "imag";

    else return "real";
}


SIGNAL_PAIR get_signal_pairs_from_swath(
    const std::string& filename,
    const std::string& swath_name
) {
    Swath swath(filename, swath_name);

    SIGNAL_PAIR signal_pair;
    signal_pair.signals = std::move(swath.get_all_signals());
    signal_pair.replica_chirps = std::move(swath.get_all_replica_chirps());

    return signal_pair;
}


F_VEC_1D scale(const CF_VEC_2D& signal, const std::string& scaling_mode)
{
    int rows = signal.size();
    int cols = signal[0].size();
    
    F_VEC_1D samples(rows*cols);

    if      (scaling_mode == "norm_log") samples = flatten(norm_2d(signal, true));
    else if (scaling_mode == "norm"    ) samples = flatten(norm_2d(signal, false));
    else if (scaling_mode == "mag"     ) samples = flatten(magnitude_2d(signal));    
    else if (scaling_mode == "real" or scaling_mode == "imag")
    {
        bool real = scaling_mode == "real";
        int  size = signal.size() * signal[0].size();
        
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                samples[i*cols+j] = real ? signal[i][j].real() : signal[i][j].imag();
            }
        }
    }
    else
    {
        throw std::invalid_argument(scaling_mode + " is not a valid scaling mode.");
    }
    return samples;
}


void range_compressed_swath(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath_name,
    const std::string& scaling_mode
) {
    std::cout << "Decoding Complex Samples and Replica Chirps" << std::endl;

    SIGNAL_PAIR signal_pair = get_signal_pairs_from_swath(in_filename, swath_name);

    int num_signals = signal_pair.signals.size();

    std::chrono::time_point compression_start = std::chrono::high_resolution_clock::now();

    CF_VEC_2D pulse_compressed(num_signals);

    for (int i = 0; i < num_signals; i++)
    {
        pulse_compressed[i] = pulse_compression(signal_pair.signals[i], signal_pair.replica_chirps[i]);
    }

    std::chrono::time_point compression_end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<float> compression_time = compression_end - compression_start;

    std::cout << "Pulse compression completed in " << compression_time.count() << "s." << std::endl;

    int rows = pulse_compressed.size();
    int cols = pulse_compressed[0].size();

    F_VEC_1D scaled = scale(pulse_compressed, scaling_mode);

    write_tif(scaled, rows, cols, out_filename);
}


void range_compressed_burst(
    const std::string& in_filename,
    const std::string& out_filename,
    const std::string& swath,
    const int&         burst_num,
    const std::string& scaling_mode
) {
    std::cout << "Decoding Burst" << std::endl;

    Burst burst(in_filename, swath, burst_num);

    int num_packets = burst.get_num_packets();

    std::cout << "Decoded " << num_packets << " signals and replica chirps." << std::endl;

    CF_VEC_2D pulse_compressed(num_packets);

    for (int i = 0; i < num_packets; i++)
    {
        pulse_compressed[i] = pulse_compression(burst.get_signal(i), burst.get_replica_chirp(i));
    }
    std::cout << "Pulse compression completed." << std::endl;

    int rows = pulse_compressed.size();
    int cols = pulse_compressed[0].size();

    F_VEC_1D scaled = scale(pulse_compressed, scaling_mode);

    write_tif(scaled, rows, cols, out_filename);
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

    range_compressed_burst(in_path, out_path, swath, burst_num, scaling);
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

    range_compressed_swath(in_path, out_path, swath, scaling);
}


void print_help(STRING_VEC_1D help_strings)
{
    for (std::string help_string : help_strings)
    {
        std::cout << help_string << std::endl;
    }
}


int main(int argc, char* argv[])
{
    STRING_VEC_1D help_strings = {
        "range_compressed_burst [swath] [burst_num] [in_path] [out_path]",
        "range_compressed_swath [swath] [in_path] [out_path]",
        "Scaling Options: [--norm_log|--norm|--mag|--real|--imag]"
    };

    if(argv[1] == __null) 
    {
        std::cout << "Please enter a command:" << std::endl;
        print_help(help_strings);
        return 1;
    }

    fftwf_init_threads();
    fftwf_plan_with_nthreads(omp_get_max_threads());

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

    if      (command == "range_compressed_burst") range_compressed_burst_command(&(argv[0]), options);
    else if (command == "range_compressed_swath") range_compressed_swath_command(&(argv[0]), options);
    else if (command == "help" or command == "--help" or command == "-h")
    {
        print_help(help_strings);
    }
    else
    {
        std::cout << command << " is not a valid command." << std::endl;
    }

    fftwf_cleanup_threads();

    return 0;
}