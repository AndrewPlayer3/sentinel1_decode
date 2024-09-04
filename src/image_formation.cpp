#include "image_formation.h"


CF_VEC_1D get_reference_function(const CF_VEC_1D& replica_chirp, const int& range_samples)
{
    int num_samples = replica_chirp.size();

    CF_VEC_1D reference = replica_chirp;

    reference.resize(range_samples);

    std::for_each(
        reference.begin() + num_samples, reference.end(),
            [] (std::complex<float>& n) { n = 0.0; }
    );

    CF_VEC_1D rep_fft = conjugate(reference); // compute_1d_dft(conjugate(reference), 0, false);

    F_VEC_1D norm = magnitude_1d(replica_chirp);

    float norm_size = norm.size();
    std::for_each(
        norm.begin(), norm.end(),
            [norm_size](std::complex<float> n) { n *= n / norm_size; }
    );
    float energy = std::accumulate(norm.begin(), norm.end(), 0.0);

    apply_hanning_window_in_place(rep_fft);

    std::for_each(
        rep_fft.begin(), rep_fft.end(),
            [energy](std::complex<float> &n) { n /= energy; }
    );

    return compute_1d_dft(rep_fft, 0, false);
}


CF_VEC_1D pulse_compression(
    const CF_VEC_1D& signal,
    const CF_VEC_1D& replica_chirp
) {
    int num_samples     = signal.size();
    int replica_samples = replica_chirp.size();

    CF_VEC_1D signal_fft = compute_1d_dft(signal,  0, false);
    CF_VEC_1D reference  = get_reference_function(replica_chirp, num_samples);

    std::transform(
        reference.begin(), reference.end(),
            signal_fft.begin(), signal_fft.begin(),
                [] (std::complex<float>& n, std::complex<float>& r) { return n * r;}
    );
    return compute_1d_dft(signal_fft, num_samples, true);
}


SIGNAL_PAIR get_signal_pairs_from_swath(
    const std::string& filename,
    const std::string& swath_name
) {
    std::ifstream data = open_file(filename);
    return get_signal_pairs_from_swath(data, swath_name);
}


SIGNAL_PAIR get_signal_pairs_from_swath(
    std::ifstream&     data,
    const std::string& swath_name
) {
    Swath swath(data, swath_name);

    SIGNAL_PAIR signal_pair;
    signal_pair.signals = std::move(swath.get_all_signals());
    signal_pair.replica_chirps = std::move(swath.get_all_replica_chirps());

    return signal_pair;
}


// TODO: These functions have quite a bit of duplication now.
CF_VEC_2D range_compress_swath(
    const std::string& filename,
    const std::string& swath_name
) {
    std::ifstream data = open_file(filename);
    return range_compress_swath(data, swath_name);
}


CF_VEC_2D range_compress_swath(
    std::ifstream&     data,
    const std::string& swath_name
) {
    SIGNAL_PAIR signal_pair = get_signal_pairs_from_swath(data, swath_name);

    int num_signals = signal_pair.signals.size();

    CF_VEC_2D pulse_compressed(num_signals);

    for (int i = 0; i < num_signals; i++)
    {
        pulse_compressed[i] = pulse_compression(signal_pair.signals[i], signal_pair.replica_chirps[i]);
    }
    return pulse_compressed;
}


CF_VEC_2D range_compress_burst(
    const std::string& filename,
    const std::string& swath,
    const int&         burst_num
) {
    std::ifstream data = open_file(filename);
    return range_compress_burst(data, swath, burst_num);
}


CF_VEC_2D range_compress_burst(
    std::ifstream&     data,
    const std::string& swath,
    const int&         burst_num
) {
    Burst burst(data, swath, burst_num);

    int num_packets = burst.get_num_packets();

    CF_VEC_2D pulse_compressed(num_packets);

    for (int i = 0; i < num_packets; i++)
    {
        pulse_compressed[i] = pulse_compression(burst.get_signal(i), burst.get_replica_chirp(i));
    }
    return pulse_compressed;
}


CF_VEC_2D range_compress_burst(
    Burst& burst
) {
    int num_packets = burst.get_num_packets();

    CF_VEC_2D pulse_compressed(num_packets);

    for (int i = 0; i < num_packets; i++)
    {
        pulse_compressed[i] = pulse_compression(burst.get_signal(i), burst.get_replica_chirp(i));
    }
    return pulse_compressed;
}


CF_VEC_2D range_doppler_swath(
    const std::string& filename,
    const std::string& swath_name
) {
    std::ifstream data = open_file(filename);
    return range_doppler_swath(data, swath_name);
}


CF_VEC_2D range_doppler_swath(
    std::ifstream&     data,
    const std::string& swath_name
) {
    CF_VEC_2D range_compressed = range_compress_swath(data, swath_name);
    return compute_axis_dft(range_compressed, 0, 0, false);
}


CF_VEC_2D range_doppler_burst(
    const std::string& filename,
    const std::string& swath_name,
    const int&         burst_num
) {
    std::ifstream data = open_file(filename);
    return range_doppler_burst(data, swath_name, burst_num);
}


CF_VEC_2D range_doppler_burst(
    std::ifstream&     data,
    const std::string& swath_name,
    const int&         burst_num
) {
    CF_VEC_2D range_compressed = range_compress_burst(data, swath_name, burst_num);
    return compute_axis_dft(range_compressed, 0, 0, false);
}


std::vector<std::complex<double>> get_average_cross_correlation(
    const CF_VEC_2D& signals
) {
    int rows = signals.size();
    int cols = signals[0].size();

    std::vector<std::complex<double>> out_row(cols);

    for (int i = 0; i < rows - 1; i++)
    {
        CF_VEC_1D conj = conjugate(signals[i+1]);
        CF_VEC_1D signal = signals[i];
        for (int j = 0; j < cols; j++)
        {
            out_row[j] += signal[j] * conj[j];
        }
        
    }
    return out_row;
}


float get_fine_dc_burst(const CF_VEC_2D& signals, const float& prf)
{
    std::vector<std::complex<double>> accc = get_average_cross_correlation(signals);
    std::complex<double> average_accc = std::complex<double>(0.0, 0.0);
    average_accc = std::accumulate(accc.begin(), accc.end(), average_accc);
    average_accc /= accc.size();
    std::cout << "Average ACCC: " << average_accc << std::endl;
    double angle = std::tan(average_accc.imag() / average_accc.real());
    float fine_dc = -1.0f * ( prf / (2.0f * PI) ) * angle;
    return fine_dc;
}


double get_velocity(
    PACKET_VEC_1D& packets,
    const int& dict_index = 0
) {
    std::vector<std::unordered_map<std::string, double>> dicts = build_data_word_dicts(packets);

    double v_x = dicts[dict_index].at("x_axis_velocity");
    double v_y = dicts[dict_index].at("y_axis_velocity");
    double v_z = dicts[dict_index].at("z_axis_velocity");

    return std::sqrt(v_x*v_x + v_y*v_y + v_z*v_z);
}


std::vector<double> get_velocities(
    PACKET_VEC_1D& packets
) { 
    std::vector<std::unordered_map<std::string, double>> dicts = build_data_word_dicts(packets);

    std::vector<double> velocities(dicts.size() - 1);

    for (int i = 0; i < dicts.size() - 1; i++)
    {
        double v_x = dicts[i].at("x_axis_velocity");
        double v_y = dicts[i].at("y_axis_velocity");
        double v_z = dicts[i].at("z_axis_velocity");

        velocities[i] = std::sqrt(v_x*v_x + v_y*v_y + v_z*v_z);
    }

    return velocities;
}


CF_VEC_2D azimuth_compress_swath(
    Swath& swath
) {
    CF_VEC_2D azimuth_compressed_swath;
    int num_bursts = swath.get_num_bursts();
    for (int i = 0; i < num_bursts; i++)
    {
        Burst burst = swath.get_burst(i);
        PACKET_VEC_1D packets = burst.get_packets();

        std::cout << "Range compressing Burst " << i << " of " << num_bursts << std::endl;
        CF_VEC_2D range_compressed_burst = range_compress_burst(burst);

        std::cout << "Azimuth Compressing Burst " << i << " of " << num_bursts << std::endl;
        CF_VEC_2D azimuth_compressed_burst = azimuth_compress(packets, range_compressed_burst);

        std::cout << "Adding Azimuth Compressed Burst " << i << " of " << num_bursts 
                  << " to the Output Vector."<< std::endl;
        for (int j = 0; j < azimuth_compressed_burst.size(); j++)
        {
            azimuth_compressed_swath.push_back(azimuth_compressed_burst[j]);
        }
    }
    return azimuth_compressed_swath;
}


CF_VEC_2D azimuth_compress(
    PACKET_VEC_1D& packets,
    CF_VEC_2D& signals
) {

    std::cout << "Getting Azimuth FM Descriptors from Header Information" << std::endl;
    int num_packets = packets.size();

    L0Packet initial_packet = packets[0];

    int   rank = initial_packet.secondary_header("rank");
    float pl   = initial_packet.get_pulse_length() * 0.000001;
    float pri  = initial_packet.get_pri() * 0.000001;
    float prf  = 1 / pri;
    float swst = initial_packet.get_swst() * 0.000001;

    std::cout << "Converting to Range Doppler Domain" << std::endl;
    CF_VEC_2D signals_rd = compute_axis_dft(signals, 0, 0, false);

    std::cout << "Computing Fine DC" << std::endl;
    float fine_dc = get_fine_dc_burst(signals_rd, prf);
    float abs_dc  = std::abs(fine_dc);

    std::cout << "Converting back to Range Azimuth Domain" << std::endl;
    CF_VEC_2D signals_ra = compute_axis_dft(signals_rd, 0, 0, true);

    float slant_range_time_near = rank * pri + swst + DELTA_T_SUPPRESSED;
    float slant_range_time_far  = slant_range_time_near + pl;
    float slant_range_near      = (slant_range_time_near * SPEED_OF_LIGHT) / 2;
    float slant_range_far       = (slant_range_time_far  * SPEED_OF_LIGHT) / 2;
    float slant_range           = (slant_range_near + slant_range_far) / 2;
    float sample_range_start    = fine_dc - (prf / 2);
    float sample_range_end      = fine_dc + (prf / 2);

    float delta = (sample_range_start + sample_range_end) / num_packets;

    std::cout << "Rank: "             << rank             << std::endl;
    std::cout << "Pulse Length: "     << pl               << std::endl;
    std::cout << "PRI: "              << pri              << std::endl;
    std::cout << "PRF: "              << prf              << std::endl;
    std::cout << "SWST: "             << swst             << std::endl;
    std::cout << "Fine DC: "          << fine_dc          << std::endl;
    std::cout << "ABS Fine DC: "      << abs_dc           << std::endl;
    std::cout << "Slant Range Near: " << slant_range_near << std::endl;
    std::cout << "Slant Range Far: "  << slant_range_far  << std::endl;
    std::cout << "Slant Range: "      << slant_range      << std::endl;

    F_VEC_1D azimuth_sample_range(num_packets);
    for (int i = 0; i < num_packets; i++)
    {
        if (i == 0) azimuth_sample_range[i] = sample_range_start;
        else azimuth_sample_range[i] = azimuth_sample_range[i-1] + delta;
    }
    double V = get_velocity(packets);

    std::cout << "Computing D(f, V)" << std::endl;
    CF_VEC_1D D(num_packets);
    for (int i = 0; i < num_packets; i++)
    {
        float range_migration_factor = 1.0f - (
            std::pow(SPEED_OF_LIGHT, 2) * std::pow(azimuth_sample_range[i], 2) / (4.0f * V * V * abs_dc * abs_dc)
        );
        D[i] = std::sqrt(std::complex<float>(range_migration_factor, 0.0f));
    }

    std::cout << "Computing the Azimuth Match Filters" << std::endl;
    CF_VEC_1D azimuth_match_filters(num_packets);
    for (int i = 0; i < num_packets; i++)
    {
        azimuth_match_filters[i] = std::exp(
            (I * 4.0f * PI * slant_range * std::complex<float>(std::pow(D[i], 2)) * abs_dc) / SPEED_OF_LIGHT
        );
    }
    azimuth_match_filters = conjugate(azimuth_match_filters);

    std::cout << "Transposing and Computing Range (Azimuth) DFT" << std::endl;
    CF_VEC_2D range_azimuth = transpose(signals);
    CF_VEC_1D match_filter  = compute_1d_dft(azimuth_match_filters, 0, false);

    std::cout << "Applying Match Filter" << std::endl;
    int rows = range_azimuth.size();
    int cols = range_azimuth[0].size();

    CF_VEC_2D azimuth_compressed = CF_VEC_2D(rows, CF_VEC_1D(cols));
    
    for (int row = 0; row < rows; row++)
    {
        CF_VEC_1D match_filtered(cols);
        CF_VEC_1D range_az = compute_1d_dft(range_azimuth[row], 0, false);

        std::transform(
            range_az.begin(), range_az.end(),
                match_filter.begin(), match_filtered.begin(),
                    [] (std::complex<float>& n, std::complex<float>& m) { return n * m;}
        );    
        match_filtered =  compute_1d_dft(match_filtered, 0, true);

        std::transform(
            match_filtered.begin(), match_filtered.end(),
                azimuth_compressed[row].begin(),
                    [] (std::complex<float>& n) {return n;}
        );
        azimuth_compressed[row] = compute_1d_dft(azimuth_compressed[row], 0, true);
    }
    std::cout << "Transposing after Azimuth Match Filtering" << std::endl;
    CF_VEC_2D azimuth_compressed_trans = transpose(azimuth_compressed);

    return azimuth_compressed_trans;
}
