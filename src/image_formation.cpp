#include "image_formation.h"


CF_VEC_1D get_reference_function(const CF_VEC_1D& replica_chirp)
{
    int num_samples = replica_chirp.size();

    CF_VEC_1D weighted_chirp = apply_hanning_window(conjugate(replica_chirp));
    F_VEC_1D norm = magnitude_1d(weighted_chirp);

    float energy = 0.0;
    for (int i = 0; i < num_samples; i++)
    {
        energy += (norm[i] * norm[i]);
    }
    energy /= norm.size();

    CF_VEC_1D reference(num_samples);
    for (int i = 0; i < num_samples; i++)
    {
        weighted_chirp[i] /= energy;
    }
    return weighted_chirp;
}


CF_VEC_1D pulse_compression(
    const CF_VEC_1D& signal,
    const CF_VEC_1D& replica_chirp
) {
    int num_samples     = signal.size();
    int replica_samples = replica_chirp.size();

    CF_VEC_1D signal_fft = compute_1d_dft(signal,  0, false);
    CF_VEC_1D reference  = get_reference_function(replica_chirp);

    reference.resize(num_samples);

    for (int i = replica_samples; i < num_samples; i++)
    {
        reference[i] = 0.0;
    }

    reference = compute_1d_dft(reference,  0, false);

    for (int i = 0; i < num_samples; i++)
    {
        signal_fft[i] *= reference[i];
    }    
    return compute_1d_dft(signal_fft, num_samples, true);
}


std::vector<std::complex<double>> get_average_cross_correlation(
    const CF_VEC_2D& signals
) {
    int rows = signals.size();
    int cols = signals[0].size();
    std::vector<std::complex<double>> out_row(cols);
    for (int i = 0; i < cols; i++)
    {
        out_row[i] = std::complex<double>(0.0, 0.0);
    }
    for (int i = 0; i < rows - 1; i++)
    {
        std::vector<std::complex<float>> conj = conjugate(signals[i+1]);
        for (int j = 0; j < cols; j++)
        {
            out_row[j] += std::complex<double>(signals[i][j]) * std::complex<double>(conj[j]);
        }
    }
    return out_row;
}


float get_fine_dc_burst(const CF_VEC_2D& signals, const float& prf)
{
    std::vector<std::complex<double>> accc = get_average_cross_correlation(signals);
    std::complex<double> average_accc = std::complex<double>(0.0, 0.0);
    for (int i = 0; i < accc.size(); i++)
    {
        average_accc += accc[i];
    }
    std::cout << "Average ACCC: " << average_accc << std::endl;
    average_accc /= accc.size();
    std::cout << "Average ACCC: " << average_accc << std::endl;
    double angle = std::tan(average_accc.imag() / average_accc.real());
    float fine_dc = -1.0f * ( prf / (2.0f * PI) ) * angle;
    return fine_dc;
}


CF_VEC_2D azimuth_compress(
    PACKET_VEC_1D& packets,
    CF_VEC_2D& signals
) {

    std::cout << "Getting Azimuth FM Descriptors from Header Information" << std::endl;
    int num_packets = packets.size();

    L0Packet initial_packet = packets[0];

    int  rank = initial_packet.secondary_header("rank");
    float pl  = initial_packet.get_pulse_length() * 0.000001;
    float pri = initial_packet.get_pri() * 0.000001;
    float prf = 1 / pri;
    float swst = initial_packet.get_swst() * 0.000001;

    std::cout << "Converting to Range Doppler Domain" << std::endl;
    CF_VEC_2D signals_rd = compute_axis_dft(signals, 0, 0, false);

    std::cout << "Computing Fine DC" << std::endl;
    float fine_dc = get_fine_dc_burst(signals_rd, prf);

    std::cout << "Converting back to Range Azimuth Domain" << std::endl;
    CF_VEC_2D signals_ra = compute_axis_dft(signals_rd, 0, 0, true);

    float slant_range_time_near = rank * pri + swst + DELTA_T_SUPPRESSED;
    float slant_range_time_far  = slant_range_time_near + pl;
    float slant_range_near = (slant_range_time_near * SPEED_OF_LIGHT) / 2;
    float slant_range_far  = (slant_range_time_far  * SPEED_OF_LIGHT) / 2;
    float slant_range = (slant_range_near + slant_range_far) / 2;

    float abs_dc = std::abs(fine_dc);
    float sample_range_start = fine_dc - (prf / 2);
    float sample_range_end   = fine_dc + (prf / 2);

    float delta = (sample_range_start + sample_range_end) / num_packets;

    std::cout << "Rank: " << rank << std::endl;
    std::cout << "Pulse Length: " << pl << std::endl;
    std::cout << "PRI: " << pri << std::endl;
    std::cout << "PRF: " << prf << std::endl;
    std::cout << "SWST: " << swst << std::endl;

    std::cout << "Fine DC: " << fine_dc << std::endl;
    std::cout << "ABS Fine DC: " << abs_dc << std::endl;

    std::cout << "Slant Range Near: " << slant_range_near << std::endl;
    std::cout << "Slant Range Far: " << slant_range_far << std::endl;
    std::cout << "Slant Range: " << slant_range << std::endl;

    F_VEC_1D azimuth_sample_range(num_packets);
    for (int i = 0; i < num_packets; i++)
    {
        if (i == 0) azimuth_sample_range[i] = sample_range_start;
        else azimuth_sample_range[i] = azimuth_sample_range[i-1] + delta;
    }
    float V = 7594.502898410107;

    std::cout << "Computing D(f, V)" << std::endl;
    CF_VEC_1D D(num_packets);
    for (int i = 0; i < num_packets; i++)
    {
        float value = 1.0f - (std::powf(SPEED_OF_LIGHT, 2) * std::powf(azimuth_sample_range[i], 2)) / (4.0f * std::powf(V, 2) * std::powf(abs_dc, 2));
        D[i] = std::sqrt(std::complex<float>(value, 0.0f));
    }

    std::cout << "Computing the Azimuth Match Filters" << std::endl;
    CF_VEC_1D azimuth_match_filters(num_packets);
    for (int i = 0; i < num_packets; i++)
    {
        azimuth_match_filters[i] = std::exp((I * 4.0f * PI * slant_range * std::complex<float>(std::pow(D[i], 2)) * abs_dc) / SPEED_OF_LIGHT);
    }

    std::cout << "Transposing and Computing Range (Azimuth) DFT" << std::endl;
    CF_VEC_2D range_azimuth_t = transpose(signals);
    CF_VEC_2D range_azimuth = compute_axis_dft(range_azimuth_t, 0, 1, false);
    CF_VEC_1D match_filter = compute_1d_dft(conjugate(azimuth_match_filters), 0, false);


    std::cout << "Applying Match Filter" << std::endl;
    int rows = range_azimuth.size();
    int cols = range_azimuth[0].size();
    CF_VEC_2D azimuth_compressed = CF_VEC_2D(rows, CF_VEC_1D(cols));
    for (int row = 0; row < rows; row++)
    {
        CF_VEC_1D match_filtered(cols);

        // std::cout << match_filtered.size() << std::endl;
        // std::cout << match_filter.size() << std::endl;
        // std::cout << azimuth_compressed.size() << std::endl;
        // std::cout << range_azimuth.size() << std::endl;

        for (int col = 0; col < cols; col++)
        {
            match_filtered[col] = range_azimuth[row][col] * match_filter[col];
        }
        CF_VEC_1D match_filtered_ifft =  compute_1d_dft(match_filtered, 0, true);

        // std::cout << std::endl;
        // std::cout << match_filtered.size() << std::endl;
        // std::cout << std::endl;

        for (int col = 0; col < cols; col++)
        {
           azimuth_compressed[row][col] = match_filtered_ifft[col];
           // std::cout << azimuth_compressed[row][col] << " ";
        }
        // std::cout << std::endl;
    }
    std::cout << "Transposing after Azimuth Match Filtering" << std::endl;
    CF_VEC_2D azimuth_compressed_trans = transpose(azimuth_compressed);

    std::cout << "Final Inverse DFT for Azimuth Compression" << std::endl;
    CF_VEC_2D output = compute_axis_dft(azimuth_compressed_trans, 0, 0, true);

    return output;
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
