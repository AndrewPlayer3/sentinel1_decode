#include "image_formation.h"


CF_VEC_1D get_reference_function(const CF_VEC_1D& replica_chirp, const int& range_samples)
{
    int num_samples = replica_chirp.size();

    CF_VEC_1D reference = replica_chirp;
    apply_hanning_window_in_place(reference);
    reference.resize(range_samples);

    F_VEC_1D norm = magnitude_1d(replica_chirp);

    double norm_size = norm.size();
    std::for_each(
        norm.begin(), norm.end(),
            [norm_size](std::complex<double> n) { n *= n / norm_size; }
    );
    double energy = std::accumulate(norm.begin(), norm.end(), 0.0);

    compute_1d_dft_in_place(reference, 0, false);
    conjugate_in_place(reference);

    std::for_each(
        reference.begin(), reference.end(),
            [energy](std::complex<double> &n) { n /= energy; }
    );

    return reference;
}


CF_VEC_1D pulse_compression(
    const CF_VEC_1D& signal,
    const CF_VEC_1D& replica_chirp
) {
    int num_samples     = signal.size();
    int replica_samples = replica_chirp.size();

    int fft_pad_amount = int(std::pow(2, std::ceil(std::log2(num_samples))) - num_samples);
    int num_fft_samples = num_samples + fft_pad_amount;

    CF_VEC_1D signal_ = signal;
    signal_.resize(num_fft_samples);
    compute_1d_dft_in_place(signal_,  0, false);

    CF_VEC_1D reference = get_reference_function(replica_chirp, num_fft_samples);

    std::transform(
        reference.begin(), reference.end(),
            signal_.begin(), signal_.begin(),
                [] (std::complex<double>& n, std::complex<double>& r) { return n * r;}
    );
    compute_1d_dft_in_place(signal_, 0, true);

    std::rotate(signal_.begin(), signal_.end()-(replica_samples / 2), signal_.end());

    return CF_VEC_1D(signal_.begin(), signal_.begin() + num_samples);
}


CF_VEC_1D pulse_compression_in_place(
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
                [] (std::complex<double>& n, std::complex<double>& r) { return n * r;}
    );
    compute_1d_dft_in_place(signal_fft, num_samples, true);
    return signal_fft;
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
    Swath swath(filename, swath_name);
    return range_compress_swath(swath);
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


CF_VEC_2D range_compress_swath(
    Swath& swath
) {
    CF_VEC_2D range_compressed_swath;

    int num_bursts = swath.get_num_bursts();

    for (int i = 0; i < num_bursts; i++)
    {
        Burst burst = swath.get_burst(i);
        PACKET_VEC_1D packets = burst.get_packets();

        std::cout << "Range compressing Burst " << i << " of " << num_bursts << std::endl;
        CF_VEC_2D range_compressed_burst = range_compress_burst(burst);

        std::cout << "Adding Range Compressed Burst " << i << " of " << num_bursts 
                  << " to the Output Vector."<< std::endl;
        range_compressed_swath.insert(range_compressed_swath.end(), range_compressed_burst.begin(), range_compressed_burst.end());
    }
    return range_compressed_swath;
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
    Swath swath(data, swath_name);

    CF_VEC_2D range_doppler_swath;

    int num_bursts = swath.get_num_bursts();

    for (int i = 0; i < num_bursts; i++)
    {
        Burst burst = swath.get_burst(i);

        std::cout << "Range compressing Burst " << i << " of " << num_bursts << std::endl;
        CF_VEC_2D range_compressed_burst = range_compress_burst(burst);

        CF_VEC_2D range_doppler_burst = range_compressed_burst;
        int num_az_samples = range_compressed_burst.size();
        int num_range_samples = range_compressed_burst[0].size();
        int num_az_fft_samples = int(std::pow(2, std::ceil(std::log2(num_az_samples))));
        int fft_pad_amount = num_az_fft_samples - num_az_samples;

        //std::cout << "Num Samples Calculated: " << num_az_fft_samples << " " << num_range_samples << std::endl;

        //std::cout << "Azimuth Padding and FFT - Burst " << i << " of " << num_bursts << std::endl;
        range_doppler_burst.resize(num_az_fft_samples, CF_VEC_1D(num_range_samples));
        //std::cout << "After Resize: " << range_doppler_burst.size() << " " << range_doppler_burst[0].size() << std::endl;
        range_doppler_burst = transpose(range_doppler_burst);
        //std::cout << "After Transpose: " << range_doppler_burst.size() << " " << range_doppler_burst[0].size() << std::endl;
        range_doppler_burst = compute_axis_dft(range_doppler_burst, 0, 1, false);
        //std::cout << "After FFT: " << range_doppler_burst.size() << " " << range_doppler_burst[0].size() << std::endl;

        std::for_each(
            range_doppler_burst.begin(), range_doppler_burst.end(),
                [] (CF_VEC_1D& row) { fftshift(row); }
        );

        range_doppler_burst = transpose(range_doppler_burst);

        //std::cout << "After FFTSHIFT: " << range_doppler_burst.size() << " " << range_doppler_burst[0].size() << std::endl;

        std::cout << "Adding Range Doppler Burst " << i << " of " << num_bursts 
                  << " to the Output Vector."<< std::endl;
        range_doppler_swath.insert(
            range_doppler_swath.end(), range_doppler_burst.begin(), range_doppler_burst.end()
        );
    }

    return range_doppler_swath;
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


CF_VEC_1D get_accc_range_blocks(
    const CF_VEC_1D& avg_cross_corr,
    const int&       num_azimuth
) {
    CF_VEC_2D accc_blocks;
    CF_VEC_1D block;

    int block_size = static_cast<int>(std::ceil(
        static_cast<double>(avg_cross_corr.size()) / static_cast<double>(num_azimuth)
    ));

    for (int i = 0; i < avg_cross_corr.size(); i++)
    {
        if (i % block_size == 0 and i != 0)
        {
            accc_blocks.push_back(block);
            block = {};
        }
        block.push_back(avg_cross_corr[i]);
    }
    accc_blocks.push_back(block);

    std::cout << "Block Size: " << block_size << std::endl;
    std::cout << "ACCC Blocks Size: " << accc_blocks.size() << std::endl;

    CF_VEC_1D averaged_blocks(accc_blocks.size());

    std::transform(
        accc_blocks.begin(), accc_blocks.end(), averaged_blocks.begin(),
            [] (const CF_VEC_1D& v) {
                return std::accumulate(v.begin(), v.end(), std::complex<double>(0.0, 0.0)) / static_cast<double>(v.size());
            }
    );

    return averaged_blocks;
}


CF_VEC_1D get_average_cross_correlation(
    const CF_VEC_2D& signals
) {
    int num_azimuth = signals.size();
    int num_range = signals[0].size();

    std::vector<std::complex<double>> accc(num_range, std::complex<double>(0, 0));

    for (int i = 0; i < num_azimuth - 1; i++)
    {
        const CF_VEC_1D& signal = signals[i];
        const CF_VEC_1D& next_signal = signals[i + 1];

        for (int j = 0; j < num_range; j++)
        {
            accc[j] += (signal[j] * std::conj(next_signal[j])) / (std::abs(signal[j]) * std::abs(next_signal[j]));
        }
    }

    return get_accc_range_blocks(accc, num_azimuth);
};


// CF_VEC_1D unwrap_fine_dc_estimates(
//     const CF_VEC_1D& fine_dcs,
//     const double&    burst_time,
//     const double&     pulse_length,
//     const double&     prf,
//     const int&       num_azimuth,
//     const int&       num_range
// ) {
//     CF_VEC_1D u(num_azimuth);
//     for (int i = 0; i < num_azimuth; i++)
//     {
//         u[i] = std::exp(I * 2.0 * PI * fine_dcs[i] / prf);
//     }
//     compute_1d_dft_in_place(u, 0, false);

//     double u_max = 0.0;
//     for (const std::complex<double>& u_val : u)
//     {
//         if (std::norm(u_val) > u_max) u_max = std::norm(u_val);
//     }
//     double a = u_max / burst_time;
//     double b = std::tan(u_max / (2.0 * PI));

//     F_VEC_1D time = linspace(0.0, pulse_length, num_azimuth);

//     CF_VEC_1D residual(num_azimuth);
//     std::transform(
//         time.begin(), time.end(), u.begin(),
//             residual.begin(),
//                 [a, b] (const double& t, const std::complex<double>& n) 
//                 { 
//                     return std::atan(n * std::exp(-1.0 * I * (a * t + b)) / (2 * PI));
//                 }
//     );

//     CF_VEC_1D unwrapped_fine_dcs(num_azimuth);
//     std::transform(
//         residual.begin(), residual.end(), time.begin(),
//             unwrapped_fine_dcs.begin(),
//                 [a, b, prf] (const std::complex<double>& r, const double& t) 
//                 { 
//                     return (a * t + b + r) / prf;
//                 }
//     );
//     return unwrapped_fine_dcs;
// }


CF_VEC_1D unwrap_fine_dc_estimates(
    const CF_VEC_1D& fine_dcs,
    const double&    burst_time,
    const float&     pulse_length,
    const float&     prf,
    const int&       num_azimuth,
    const int&       num_range
) {
    std::complex<double> imag(0.0, 1.0);
    
    // Compute complex exponentials
    CF_VEC_1D u(num_azimuth);
    for (int i = 0; i < num_azimuth; i++) {
        u[i] = std::exp(imag * 2.0 * PI * fine_dcs[i] / static_cast<double>(prf));
    }
    compute_1d_dft_in_place(u, 0, false);

    // Find the maximum phase (unwrap these phases first)
    CF_VEC_1D unwrapped_phases(num_azimuth);
    for (int i = 0; i < num_azimuth; i++) {
        unwrapped_phases[i] = std::arg(u[i]); // Phase of the DFT result
    }

    // Calculate a (rate of phase change) and b (initial phase)
    F_VEC_1D time = linspace(0.0, pulse_length, num_azimuth);
    std::complex<double> a = (unwrapped_phases.back() - unwrapped_phases.front()) / (time.back() - time.front());
    std::complex<double> b = unwrapped_phases[0];

    // Compute residual (optional adjustment)
    CF_VEC_1D residual(num_azimuth);
    std::transform(
        time.begin(), time.end(), unwrapped_phases.begin(),
        residual.begin(),
        [a, b] (const double& t, const std::complex<double>& phase) {
            return phase - (a * t + b); // Compute residual
        }
    );

    // Final unwrapped Doppler centroids
    CF_VEC_1D unwrapped_fine_dcs(num_azimuth);
    std::transform(
        residual.begin(), residual.end(), time.begin(),
        unwrapped_fine_dcs.begin(),
        [a, b, prf] (const std::complex<double>& r, const double& t) {
            return (a * t + b + r) / static_cast<double>(prf); // Adjust phase
        }
    );

    return unwrapped_fine_dcs;
}


CF_VEC_1D get_fine_dcs_for_burst(
    const CF_VEC_2D& signals,
    const double& burst_time,
    const double&  pulse_length,
    const double&  prf,
    const int&    num_azimuth,
    const int&    num_range
) {
    CF_VEC_1D accc = get_average_cross_correlation(signals);

    CF_VEC_1D fine_dcs(accc.size());
    std::transform(
        accc.begin(), accc.end(),
            fine_dcs.begin(),
                [prf] (std::complex<double>& c) { 
                    return -1.0 * ( prf / (2.0 * PI)) * std::tan(c.imag() / c.real());
                }
    );

    fine_dcs.resize(num_azimuth);

    CF_VEC_1D unwrapped_fine_dcs = unwrap_fine_dc_estimates(
        fine_dcs, burst_time, pulse_length,
        prf, num_azimuth, num_range
    );
    return unwrapped_fine_dcs;
}


double get_velocity(
    PACKET_VEC_1D& packets,
    const int&     dict_index = 0
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
    const std::string& filename,
    const std::string& swath_name
) {
    Swath swath(filename, swath_name);
    
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
        azimuth_compressed_swath.insert(
            azimuth_compressed_swath.end(), azimuth_compressed_burst.begin(), azimuth_compressed_burst.end()
        );
    }
    return azimuth_compressed_swath;
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
        azimuth_compressed_swath.insert(
            azimuth_compressed_swath.end(), azimuth_compressed_burst.begin(), azimuth_compressed_burst.end()
        );
    }
    return azimuth_compressed_swath;
}


CF_VEC_1D get_azimuth_matched_filter(
    const int&   num_azimuth,
    const double& doppler_centroid_est,
    const double& prf,
    const double& velocity, 
    const double& slant_range
) {
    double sample_range_start = doppler_centroid_est - (prf);
    double sample_range_end   = doppler_centroid_est + (prf);
    F_VEC_1D azimuth_sample_range = linspace(sample_range_start, sample_range_end, num_azimuth);

    CF_VEC_1D D(num_azimuth);

    double denominator = (4.0 * velocity * velocity * doppler_centroid_est * doppler_centroid_est);

    for (int i = 0; i < num_azimuth; i++)
    {
        double range_migration_factor = 1.0 - (
            std::pow(SPEED_OF_LIGHT, 2) * std::pow(azimuth_sample_range[i], 2) / denominator
        );
        D[i] = std::sqrt(std::complex<double>(range_migration_factor, 0.0));
    }
    CF_VEC_1D azimuth_match_filter(num_azimuth);

    std::complex<double> phase = I * 4.0 * PI * slant_range * doppler_centroid_est;
    for (int i = 0; i < num_azimuth; i++)
    {
        azimuth_match_filter[i] = (1.0 / num_azimuth) * std::exp(
            (phase * D[i] * D[i] * doppler_centroid_est) / SPEED_OF_LIGHT
        );
    }
    // conjugate_in_place(azimuth_match_filter);
    // apply_hanning_window_in_place(azimuth_match_filter);

    return azimuth_match_filter; // compute_1d_dft(azimuth_match_filter, 0, false);
}


CF_VEC_2D get_azimuth_matched_filters(
    const CF_VEC_1D& doppler_centroid_ests,
    const double& velocity,
    const double& pulse_length,
    const double& pri,
    const double& window_start_time,
    const int&   rank,
    const int&   num_azimuth,
    const int&   num_range
) {
    double prf = 1 / pri;

    double slant_range_time_near = rank * pri + window_start_time + DELTA_T_SUPPRESSED;
    double slant_range_time_far  = slant_range_time_near + pulse_length;
    double slant_range_near      = (slant_range_time_near * SPEED_OF_LIGHT) / 2;
    double slant_range_far       = (slant_range_time_far  * SPEED_OF_LIGHT) / 2;
    double slant_range           = (slant_range_near + slant_range_far) / 2;

    int block_size = int(std::ceil(double(num_range) / double(num_azimuth)));

    CF_VEC_1D current_match_filter(num_azimuth);
    CF_VEC_2D azimuth_match_filters(num_range, CF_VEC_1D(num_azimuth));

    int fine_dc_index = 0;
    for (int i = 0; i < num_range; i++)
    {
        if (i % block_size == 0)
        {
            double fine_dc = std::abs(doppler_centroid_ests[fine_dc_index]);
            current_match_filter = get_azimuth_matched_filter(
                num_azimuth,
                fine_dc,
                prf,
                velocity, 
                slant_range
            );
            fine_dc_index++;
        }
        azimuth_match_filters[i] = current_match_filter;
    }
    return azimuth_match_filters;
}


double gps_time_to_seconds_since_midnight(double gps_time_seconds) {
    std::chrono::system_clock::time_point gps_epoch = 
        std::chrono::system_clock::from_time_t(315964800); // Unix timestamp for Jan 6, 1980

    long long whole_seconds = static_cast<long long>(gps_time_seconds);
    double fractional_seconds = gps_time_seconds - whole_seconds;

    std::chrono::system_clock::time_point current_time = gps_epoch + std::chrono::seconds(whole_seconds);

    std::time_t current_time_t = std::chrono::system_clock::to_time_t(current_time);
    std::tm* current_time_tm = std::gmtime(&current_time_t);

    double seconds_since_midnight = 
        current_time_tm->tm_hour * 3600 +
        current_time_tm->tm_min * 60 +
        current_time_tm->tm_sec;

    seconds_since_midnight += fractional_seconds;

    return seconds_since_midnight;
}


// FFT shift function
void fftshift(std::vector<std::complex<double>>& data) {
    size_t N = data.size();

    // If N is even
    if (N % 2 == 0) {
        // Rotate the vector by N/2 to shift the FFT output
        std::rotate(data.begin(), data.begin() + N / 2, data.end());
    }
    // If N is odd
    else {
        // Rotate the vector by (N+1)/2 to shift the FFT output
        std::rotate(data.begin(), data.begin() + (N + 1) / 2, data.end());
    }
}


CF_VEC_2D azimuth_compress(
    PACKET_VEC_1D& packets,
    CF_VEC_2D& signals
) {
    const int& num_azimuth = signals.size();
    const int& num_range = signals[0].size();

    std::cout << "Getting Azimuth FM Descriptors from Header Information" << std::endl;
    int num_packets = packets.size();

    L0Packet initial_packet = packets[0];

    int   rank = initial_packet.secondary_header("rank");
    double pl   = initial_packet.get_pulse_length() * 0.000001;
    double pri  = initial_packet.get_pri() * 0.000001;
    double prf  = 1 / pri;
    double swst = initial_packet.get_swst() * 0.000001;

    double burst_time = initial_packet.get_time();

    burst_time = gps_time_to_seconds_since_midnight(burst_time);

    std::cout << "Converting to Range Doppler Domain" << std::endl;
    CF_VEC_2D signals_rd = compute_axis_dft(signals, 0, 0, false);

    std::cout << "Computing Fine DC" << std::endl;
    CF_VEC_1D fine_dcs = get_fine_dcs_for_burst(
        signals_rd,
        burst_time,
        pl,
        prf,
        num_azimuth,
        num_range
    );

    double V = get_velocity(packets);
    std::cout << "Velocity: " << V << std::endl;

    std::cout << "Computing the Azimuth Match Filters" << std::endl;
    CF_VEC_2D reference_functions = get_azimuth_matched_filters(
        fine_dcs, V, pl, pri, swst, rank, num_azimuth, num_range
    );

    std::cout << "First 10 AZMFs: " << std::endl;
    for (int i = 0; i < 10; i++) std::cout << reference_functions[0][i] << " ";
    std::cout << std::endl;

    std::cout << "Reference Functions Rows: " << reference_functions.size() << std::endl;
    std::cout << "Reference Functions Cols: " << reference_functions[0].size() << std::endl;

    std::cout << "Transposing and Computing Range (Azimuth) DFT" << std::endl;
    CF_VEC_2D range_azimuth = transpose(signals);

    std::cout << "Applying Match Filter" << std::endl;    
    for (int row = 0; row < num_range; row++)
    {
        CF_VEC_1D& azimuth_row = range_azimuth[row];
        CF_VEC_1D& reference_function = reference_functions[row];

        compute_1d_dft_in_place(azimuth_row, 0, false);

        std::transform(
            azimuth_row.begin(), azimuth_row.end(),
                reference_function.begin(), azimuth_row.begin(),
                    [] (std::complex<double>& signal, std::complex<double>& match) { return signal * match; }
        );

        compute_1d_dft_in_place(azimuth_row, 0, true);
        fftshift(azimuth_row);
        compute_1d_dft_in_place(azimuth_row, 0, true);
        fftshift(azimuth_row);
    }

    std::cout << "Transposing after Azimuth Match Filtering" << std::endl;
    return transpose(range_azimuth);
}