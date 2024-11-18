#include "image_formation.h"


CF_VEC_1D get_reference_function(const CF_VEC_1D& replica_chirp)
{
    int num_samples = replica_chirp.size();

    CF_VEC_1D reference = replica_chirp;
    apply_hanning_window_in_place(reference);

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
    const CF_VEC_1D& reference
) {
    int num_samples = signal.size();

    CF_VEC_1D signal_    = signal;

    std::transform(
        reference.begin(), reference.end(),
            signal_.begin(), signal_.begin(),
                [] (const std::complex<double>& n, std::complex<double>& r) { return n * r;}
    );

    return CF_VEC_1D(signal_.begin(), signal_.begin() + num_samples);
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
    
    if (N % 2 == 0) {
        std::rotate(data.begin(), data.begin() + N / 2, data.end());
    }
    else {
        std::rotate(data.begin(), data.begin() + (N + 1) / 2, data.end());
    }
}


// CF_VEC_1D get_average_cross_correlation(
//     const CF_VEC_2D& signals
// ) {
//     int num_azimuth = signals.size();
//     int num_range = signals[0].size();

//     std::vector<std::complex<double>> accc(num_range, std::complex<double>(0, 0));

//     for (int i = 0; i < num_azimuth - 1; i++)
//     {
//         const CF_VEC_1D& signal = signals[i];
//         const CF_VEC_1D& next_signal = signals[i + 1];

//         for (int j = 0; j < num_range; j++)
//         {
//             accc[j] += (signal[j] * std::conj(next_signal[j])) / (std::abs(signal[j]) * std::abs(next_signal[j]));
//         }
//     }

//     return get_accc_range_blocks(accc, num_azimuth);
// };


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


// CF_VEC_1D unwrap_fine_dc_estimates(
//     const CF_VEC_1D& fine_dcs,
//     const double&    burst_time,
//     const float&     pulse_length,
//     const float&     prf,
//     const int&       num_azimuth,
//     const int&       num_range
// ) {
//     std::complex<double> imag(0.0, 1.0);
    
//     // Compute complex exponentials
//     CF_VEC_1D u(num_azimuth);
//     for (int i = 0; i < num_azimuth; i++) {
//         u[i] = std::exp(imag * 2.0 * PI * fine_dcs[i] / static_cast<double>(prf));
//     }
//     compute_1d_dft_in_place(u, 0, false);

//     // Find the maximum phase (unwrap these phases first)
//     CF_VEC_1D unwrapped_phases(num_azimuth);
//     for (int i = 0; i < num_azimuth; i++) {
//         unwrapped_phases[i] = std::arg(u[i]); // Phase of the DFT result
//     }

//     // Calculate a (rate of phase change) and b (initial phase)
//     F_VEC_1D time = linspace(0.0, pulse_length, num_azimuth);
//     std::complex<double> a = (unwrapped_phases.back() - unwrapped_phases.front()) / (time.back() - time.front());
//     std::complex<double> b = unwrapped_phases[0];

//     // Compute residual (optional adjustment)
//     CF_VEC_1D residual(num_azimuth);
//     std::transform(
//         time.begin(), time.end(), unwrapped_phases.begin(),
//         residual.begin(),
//         [a, b] (const double& t, const std::complex<double>& phase) {
//             return phase - (a * t + b); // Compute residual
//         }
//     );

//     // Final unwrapped Doppler centroids
//     CF_VEC_1D unwrapped_fine_dcs(num_azimuth);
//     std::transform(
//         residual.begin(), residual.end(), time.begin(),
//         unwrapped_fine_dcs.begin(),
//         [a, b, prf] (const std::complex<double>& r, const double& t) {
//             return (a * t + b + r) / static_cast<double>(prf); // Adjust phase
//         }
//     );

//     return unwrapped_fine_dcs;
// }


// CF_VEC_1D get_fine_dcs_for_burst(
//     const CF_VEC_2D& signals,
//     const double& burst_time,
//     const double&  pulse_length,
//     const double&  prf,
//     const int&    num_azimuth,
//     const int&    num_range
// ) {
//     CF_VEC_1D accc = get_average_cross_correlation(signals);

//     CF_VEC_1D fine_dcs(accc.size());
//     std::transform(
//         accc.begin(), accc.end(),
//             fine_dcs.begin(),
//                 [prf] (std::complex<double>& c) { 
//                     return -1.0 * ( prf / (2.0 * PI)) * std::tan(c.imag() / c.real());
//                 }
//     );

//     fine_dcs.resize(num_azimuth);

//     CF_VEC_1D unwrapped_fine_dcs = unwrap_fine_dc_estimates(
//         fine_dcs, burst_time, pulse_length,
//         prf, num_azimuth, num_range
//     );
//     return unwrapped_fine_dcs;
// }


// double get_velocity(
//     PACKET_VEC_1D& packets,
//     const int&     dict_index = 0
// ) {
//     std::vector<std::unordered_map<std::string, double>> dicts = build_data_word_dicts(packets);

//     double v_x = dicts[dict_index].at("x_axis_velocity");
//     double v_y = dicts[dict_index].at("y_axis_velocity");
//     double v_z = dicts[dict_index].at("z_axis_velocity");

//     return std::sqrt(v_x*v_x + v_y*v_y + v_z*v_z);
// }


// std::vector<double> get_velocities(
//     PACKET_VEC_1D& packets
// ) { 
//     std::vector<std::unordered_map<std::string, double>> dicts = build_data_word_dicts(packets);

//     std::vector<double> velocities(dicts.size() - 1);

//     for (int i = 0; i < dicts.size() - 1; i++)
//     {
//         double v_x = dicts[i].at("x_axis_velocity");
//         double v_y = dicts[i].at("y_axis_velocity");
//         double v_z = dicts[i].at("z_axis_velocity");

//         velocities[i] = std::sqrt(v_x*v_x + v_y*v_y + v_z*v_z);
//     }

//     return velocities;
// }


// CF_VEC_1D get_azimuth_matched_filter(
//     const int&   num_azimuth,
//     const double& doppler_centroid_est,
//     const double& prf,
//     const double& velocity, 
//     const double& slant_range
// ) {
//     double sample_range_start = doppler_centroid_est - (prf);
//     double sample_range_end   = doppler_centroid_est + (prf);
//     F_VEC_1D azimuth_sample_range = linspace(sample_range_start, sample_range_end, num_azimuth);

//     CF_VEC_1D D(num_azimuth);

//     double denominator = (4.0 * velocity * velocity * doppler_centroid_est * doppler_centroid_est);

//     for (int i = 0; i < num_azimuth; i++)
//     {
//         double range_migration_factor = 1.0 - (
//             std::pow(SPEED_OF_LIGHT, 2) * std::pow(azimuth_sample_range[i], 2) / denominator
//         );
//         D[i] = std::sqrt(std::complex<double>(range_migration_factor, 0.0));
//     }
//     CF_VEC_1D azimuth_match_filter(num_azimuth);

//     std::complex<double> phase = I * 4.0 * PI * slant_range * doppler_centroid_est;
//     for (int i = 0; i < num_azimuth; i++)
//     {
//         azimuth_match_filter[i] = (1.0 / num_azimuth) * std::exp(
//             (phase * D[i] * D[i] * doppler_centroid_est) / SPEED_OF_LIGHT
//         );
//     }
//     // conjugate_in_place(azimuth_match_filter);
//     // apply_hanning_window_in_place(azimuth_match_filter);

//     return azimuth_match_filter; // compute_1d_dft(azimuth_match_filter, 0, false);
// }


// CF_VEC_2D get_azimuth_matched_filters(
//     const CF_VEC_1D& doppler_centroid_ests,
//     const double& velocity,
//     const double& pulse_length,
//     const double& pri,
//     const double& window_start_time,
//     const int&   rank,
//     const int&   num_azimuth,
//     const int&   num_range
// ) {
//     double prf = 1 / pri;

//     double slant_range_time_near = rank * pri + window_start_time + DELTA_T_SUPPRESSED;
//     double slant_range_time_far  = slant_range_time_near + pulse_length;
//     double slant_range_near      = (slant_range_time_near * SPEED_OF_LIGHT) / 2;
//     double slant_range_far       = (slant_range_time_far  * SPEED_OF_LIGHT) / 2;
//     double slant_range           = (slant_range_near + slant_range_far) / 2;

//     int block_size = int(std::ceil(double(num_range) / double(num_azimuth)));

//     CF_VEC_1D current_match_filter(num_azimuth);
//     CF_VEC_2D azimuth_match_filters(num_range, CF_VEC_1D(num_azimuth));

//     int fine_dc_index = 0;
//     for (int i = 0; i < num_range; i++)
//     {
//         if (i % block_size == 0)
//         {
//             double fine_dc = std::abs(doppler_centroid_ests[fine_dc_index]);
//             current_match_filter = get_azimuth_matched_filter(
//                 num_azimuth,
//                 fine_dc,
//                 prf,
//                 velocity, 
//                 slant_range
//             );
//             fine_dc_index++;
//         }
//         azimuth_match_filters[i] = current_match_filter;
//     }
//     return azimuth_match_filters;
// }
