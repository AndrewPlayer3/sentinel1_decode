#include "doppler.h"


CF_VEC_1D get_wrapped_estimates(CF_VEC_2D& range_compressed, const double& prf)
{
    int num_az = range_compressed.size();
    int num_rng = range_compressed[0].size();

    CF_VEC_1D phase_difference_sum(num_rng);

    for (int rng_line = 0; rng_line < num_rng; rng_line++)
    {
        for (int az_line = 0; az_line < num_az - 1; az_line++)
        {
            phase_difference_sum[rng_line] += range_compressed[az_line][rng_line] * std::conj(range_compressed[az_line+1][rng_line]);
        }
    }

    std::cout << "Wrapped F_DC: " << phase_difference_sum.front() << ", " << phase_difference_sum.back() << std::endl;

    CF_VEC_1D wrapped_estimates(num_rng);

    std::transform(
        phase_difference_sum.begin(), phase_difference_sum.end(),
            wrapped_estimates.begin(),
                [prf] (std::complex<double>& v) {
                    return -prf * std::arg(v) / (2 * PI);
                }
    );

    return wrapped_estimates;
}


F_VEC_1D get_unwrapped_estimates(CF_VEC_1D& fine_dc_estimates, F_VEC_1D& offset_rng_times, const double& prf, const int& dense_grid_scalar = 8)
{
    int num_est = fine_dc_estimates.size();
    float dt = offset_rng_times[1] - offset_rng_times[0];

    CF_VEC_1D F(num_est * dense_grid_scalar);
    std::transform(
        fine_dc_estimates.begin(), fine_dc_estimates.end(),
            F.begin(),
                [prf] (std::complex<double>& dc_est) {
                    return std::floor(std::abs(std::exp(2.0 * I * PI * dc_est / prf)));
                }
    );

    compute_1d_dft_in_place(F, 0, false);
    F_VEC_1D F_norm_squared(num_est * dense_grid_scalar);
    std::transform(
        F.begin(), F.end(),
            F_norm_squared.begin(),
                [] (std::complex<double>& n) {
                    return std::abs(n) * std::abs(n);
                }
    );

    F_VEC_1D::iterator max = std::max_element(F_norm_squared.begin(), F_norm_squared.end());
    int max_index = std::distance(F_norm_squared.begin(), max);

    F_VEC_1D freqs = fftfreq(num_est * dense_grid_scalar, dt);
    
    double v = freqs[max_index];
    double a = v / dt;
    double b = std::arg(F[max_index]) / (2.0 * PI);

    std::complex<double> r = std::arg(std::exp(2.0 * I * PI * fine_dc_estimates[0] / prf) * std::exp(-1.0 * I * (a * 0.01 + b))) / (2.0 * PI);
    std::complex<double> u = (a * 0.0 + b + r) * prf;

    F_VEC_1D unwrapped_estimates(num_est);
    std::transform(
        fine_dc_estimates.begin(), fine_dc_estimates.end(),
            offset_rng_times.begin(), unwrapped_estimates.begin(),
                [prf, a, b] (std::complex<double>& dc_est, double& time) {
                    std::complex<double> res = std::exp(2.0 * I * PI * dc_est / prf) * std::exp(-1.0 * I * (a * time + b));
                    double r = std::arg(res) / (2.0 * PI);
                    return (a * time + b + r) * prf;
                }
    );

    return unwrapped_estimates;
}


std::vector<double> polyfit(const std::vector<double>& x, const std::vector<double>& y) {
    const int N = x.size();

    // Precompute sums
    double Sx[7] = {}; // Sx[i] = sum(x^i), i from 0 to 6
    double Sy[4] = {}; // Sy[i] = sum(x^i * y)

    for (int i = 0; i < N; ++i) {
        double xi = 1.0;
        for (int j = 0; j <= 6; ++j) {
            if (j <= 3) Sy[j] += xi * y[i];
            Sx[j] += xi;
            xi *= x[i];
        }
    }

    // Construct the normal equations matrix and RHS
    double A[4][4] = {
        {Sx[0], Sx[1], Sx[2], Sx[3]},
        {Sx[1], Sx[2], Sx[3], Sx[4]},
        {Sx[2], Sx[3], Sx[4], Sx[5]},
        {Sx[3], Sx[4], Sx[5], Sx[6]}
    };
    double b[4] = {Sy[0], Sy[1], Sy[2], Sy[3]};

    // Solve Ax = b using Gaussian elimination
    for (int i = 0; i < 4; ++i) {
        // Pivot
        for (int k = i + 1; k < 4; ++k) {
            if (std::abs(A[k][i]) > std::abs(A[i][i])) {
                for (int j = 0; j < 4; ++j) std::swap(A[i][j], A[k][j]);
                std::swap(b[i], b[k]);
            }
        }

        // Eliminate below
        for (int k = i + 1; k < 4; ++k) {
            double factor = A[k][i] / A[i][i];
            for (int j = i; j < 4; ++j)
                A[k][j] -= factor * A[i][j];
            b[k] -= factor * b[i];
        }
    }

    // Back substitution
    std::vector<double> coeffs(4);
    for (int i = 3; i >= 0; --i) {
        coeffs[i] = b[i];
        for (int j = i + 1; j < 4; ++j)
            coeffs[i] -= A[i][j] * coeffs[j];
        coeffs[i] /= A[i][i];
    }

    return coeffs; // [d, c, b, a]
}


double polyval(const std::vector<double>& coeffs, double x) {
    // Evaluate polynomial: a*x^3 + b*x^2 + c*x + d
    double result = 0.0;
    for (int i = coeffs.size() - 1; i >= 0; --i)
        result = result * x + coeffs[i];
    return result;
}


CF_VEC_1D get_deramping_signal(const int& num_az_samples, const double& doppler_centroid_rate, const double& burst_duration)
{
    F_VEC_1D az_time = linspace(0.0, burst_duration, num_az_samples);

    CF_VEC_1D deramping_signal(num_az_samples);

    std::transform(
        az_time.begin(), az_time.end(),
            deramping_signal.begin(),
                [doppler_centroid_rate, burst_duration] (double& t) {
                    return std::exp(I * PI * doppler_centroid_rate * std::pow(t - (0.5 * burst_duration), 2.0));
                }
    );

    return deramping_signal;
}


CF_VEC_2D dce_preconditioning(CF_VEC_2D& range_compressed, const double& doppler_centroid_rate, const double& burst_duration)
{
    int num_az = range_compressed.size();
    int num_rng = range_compressed[0].size();

    CF_VEC_1D deramping_signal = get_deramping_signal(num_az, doppler_centroid_rate, burst_duration);

    CF_VEC_2D preconditioned = transpose(range_compressed);

    std::transform(
        preconditioned.begin(), preconditioned.end(),
            preconditioned.begin(),
                [deramping_signal] (CF_VEC_1D& precon) {
                    std::transform(
                        deramping_signal.begin(), deramping_signal.end(),
                            precon.begin(), precon.begin(),
                                [] (std::complex<double> deramp, std::complex<double> rc) {
                                    return rc * deramp;
                                }
                    );
                    return precon;
                }
    );

    return transpose(preconditioned);
}


F_VEC_1D get_doppler_centroid(CF_VEC_2D& range_compressed, const double& doppler_centroid_rate, const double& burst_duration, L0Packet& first_packet)
{
    int num_az = range_compressed.size();
    int num_rng = range_compressed[0].size();

    double prf = 1.0 / (first_packet.get_pri() * 1e-6);

    CF_VEC_2D precoditioned_burst = dce_preconditioning(range_compressed, doppler_centroid_rate, 1.06);

    CF_VEC_1D wrapped_estimates = get_wrapped_estimates(precoditioned_burst, prf);

    // TODO: Make this an argument once things are working
    int num_rng_blocks = 15;
    int rng_block_size = num_rng / num_rng_blocks;

    CF_VEC_1D fine_dc_estimates(num_rng_blocks);

    for (int rng_block_index; rng_block_index < num_rng_blocks; rng_block_index++)
    {
        int start_index = rng_block_index * rng_block_size;

        int end_index = start_index + rng_block_size < num_rng ?
                            start_index + rng_block_size 
                            :
                            num_rng - 1;

        CF_VEC_1D rng_block(wrapped_estimates.begin() + start_index, wrapped_estimates.begin() + end_index);

        fine_dc_estimates[rng_block_index] = std::accumulate(rng_block.begin(), rng_block.end(), 0.0 + 0.0 * I) / (double(end_index) - double(start_index));
    }

    std::cout << "Wrapped F_DC: " << fine_dc_estimates.front() << ", " << fine_dc_estimates.back() << std::endl;

    F_VEC_1D rng_times = first_packet.get_slant_range_times(num_rng_blocks);
    double initial_time = rng_times[0];
    std::transform(
        rng_times.begin(), rng_times.end(),
            rng_times.begin(),
                [initial_time] (double& v) {
                    return v - initial_time;
                }
    );

    F_VEC_1D unwrapped_estimates = get_unwrapped_estimates(fine_dc_estimates, rng_times, prf, 8);

    std::cout << "Unwrapped F_DC: " << unwrapped_estimates.front() << ", " << unwrapped_estimates.back() << std::endl;

    F_VEC_1D poly = polyfit(rng_times, unwrapped_estimates);

    rng_times = first_packet.get_slant_range_times(num_rng);
    initial_time = rng_times[0];
    std::transform(
        rng_times.begin(), rng_times.end(),
            rng_times.begin(),
                [initial_time] (double& v) {
                    return v - initial_time;
                }
    );

    unwrapped_estimates = F_VEC_1D(num_rng);

    std::transform(
        rng_times.begin(), rng_times.end(),
            unwrapped_estimates.begin(),
                [poly] (const double& time) {
                    return polyval(poly, time);
                } 
    );

    std::cout << "Unwrapped F_DC: " << unwrapped_estimates.front() << ", " << unwrapped_estimates.back() << std::endl;

    return unwrapped_estimates;
}


double get_doppler_centroid_rate(PACKET_VEC_1D& burst_packets, const double& velocity)
{
    double pri = burst_packets[0].get_pri() * 1e-6;

    double antenna_steering_rate = (
        (burst_packets.back().get_azimuth_beam_angle() - burst_packets.front().get_azimuth_beam_angle()) / (pri * burst_packets.size())
    );

    return (-2 * velocity * antenna_steering_rate) / WAVELENGTH;
}


// Mosaic a Signal for UFR
CF_VEC_1D get_tiled_signal(CF_VEC_1D& signal, const double& num_replicas)
{
    int num_az = signal.size();
    int shape = std::floor(num_replicas * num_az);
    int offset = std::floor(num_replicas);

    CF_VEC_1D tiled_signal(shape);

    // Mosaic `floor(num_replicas)` times
    for (int i = 0; i < offset; i++)
    {
        std::transform(
            signal.begin(), signal.end(),
                tiled_signal.begin() + i * num_az,
                    [] (std::complex<double>& sample) {
                        return sample;
                    }
        );
    }

    // If `num_replicas` isn't a whole number, tile the remaining percentage after `floor(num_replicas)`
    if (std::floor(num_replicas) != num_replicas)
    {
        int num_samples = std::floor((num_replicas - std::floor(num_replicas)) * num_az);

        std::transform(
            signal.begin(), signal.begin() + num_samples,
                tiled_signal.begin() + offset * num_az,
                    [] (std::complex<double>& sample) {
                        return sample;
                    }
        );
    }

    return tiled_signal;
}


CF_VEC_2D azimuth_frequency_ufr(
    CF_VEC_2D& range_compressed,
    F_VEC_1D&  dc_estimates,
    L0Packet&  initial_packet,
    const double& dc_rate,
    const double& burst_duration,
    const double& prf,
    const double& processing_bandwidth  // B_d
) {
    int num_az = range_compressed.size();
    int num_rng = range_compressed[0].size();

    range_compressed = transpose(range_compressed);

    compute_axis_dft_in_place(range_compressed, 0, 1, false);
    fftshift_in_place(range_compressed);

    double num_replicas = std::abs(dc_rate * burst_duration / prf);
    double bandwidth = num_replicas * prf / 2.0;  // B_s

    int shape = std::floor(num_replicas * num_az);
    int mid = shape / 2;
    int start = mid - (std::floor(processing_bandwidth) / 2);
    int end = mid + (std::floor(processing_bandwidth) / 2);

    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << "Computing Azimuth Frequency UFR with Parameters: " << std::endl;
    std::cout << "PRF: " << prf << std::endl;
    std::cout << "DC Rate: " << dc_rate << std::endl;
    std::cout << "Burst Duration: " << burst_duration << std::endl;
    std::cout << "Number of Spectral Replicas: " << num_replicas << std::endl;
    std::cout << "Ramp Signal Bandwidth: " << bandwidth << std::endl;
    std::cout << "Output Shape: " << shape << std::endl;
    std::cout << "Low Pass Filter Center: " << mid << std::endl;
    std::cout << "Low Pass Filter Start Index: " << start << std::endl;
    std::cout << "Low Pass Filter End Index: " << end << std::endl;
    std::cout << "-------------------------------------------------" << std::endl;

    F_VEC_1D az_freqs = linspace(-bandwidth, bandwidth, shape);

    CF_VEC_2D ufr_output(num_rng, CF_VEC_1D(shape));

    #pragma omp parallel for
    for (int rng_line = 0; rng_line < num_rng; rng_line++)
    {
        // Mosaic Signal
        ufr_output[rng_line] = get_tiled_signal(range_compressed[rng_line], num_replicas);

        //Apply De-ramping
        for (int az_line = 0; az_line < shape; az_line++)
        {
            std::complex<double> sample = ufr_output[rng_line][az_line];
            std::complex<double> deramp = std::exp(-I * PI * (1 / dc_rate) * std::pow(az_freqs[az_line], 2.0));

            ufr_output[rng_line][az_line] = deramp * sample;
        }
    }

    // Low-pass Filter
    compute_axis_dft_in_place(ufr_output, 0, 1, false);

    #pragma omp parallel for collapse(2)
    for (int rng_line = 0; rng_line < num_rng; rng_line++)
    {
        for (int az_line = 0; az_line < shape; az_line++)
        {
            if (az_line < end-start)
            {
                ufr_output[rng_line][az_line] = ufr_output[rng_line][start + az_line];
            }
            else
            {
                ufr_output[rng_line][az_line] = 0.0;
            }
        }
    }

    compute_axis_dft_in_place(ufr_output, 0, 1, true);

    #pragma omp parallel for collapse(2)
    for (int rng_line = 0; rng_line < num_rng; rng_line++)
    {
        // Apply Re-ramping
        for (int az_line = 0; az_line < shape; az_line++)
        {
            std::complex<double> sample = ufr_output[rng_line][az_line];
            std::complex<double> reramp = std::exp(I * PI * (1 / dc_rate) * std::pow(az_freqs[az_line], 2.0));

            ufr_output[rng_line][az_line] = reramp * sample;
        }
    }

    return transpose(ufr_output);
}


CF_VEC_2D azimuth_time_ufr(
    CF_VEC_2D& azimuth_compressed,
    F_VEC_1D&  dc_estimates,
    F_VEC_2D& ka,
    L0Packet&  initial_packet,
    const double& dc_rate,
    const double& burst_duration,
    const double& prf,
    const double& processing_bandwidth  // B_d
) {
    ka = transpose(ka);
    azimuth_compressed = transpose(azimuth_compressed);

    int num_az = azimuth_compressed[0].size();
    int num_rng = azimuth_compressed.size();

    double num_replicas = std::abs(dc_rate * burst_duration / prf);
    double bandwidth = num_replicas * prf;

    double ka_mean = std::abs(std::accumulate(ka[0].begin(), ka[0].end(), 0.0)) / num_az;
    double focused_time = burst_duration + ((bandwidth - (2 * processing_bandwidth)) / ka_mean);

    int num_tiles = std::ceil((focused_time / 2) / burst_duration) - std::floor((-focused_time / 2) / burst_duration);

    int output_shape = std::floor(prf * burst_duration);
    int downsample_shape = output_shape + std::floor(processing_bandwidth);
    int shape = num_tiles * num_az;
    int mid = shape / 2;
    int start = mid - bandwidth / 2;
    int end = mid + bandwidth / 2;

    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << "Computing Azimuth Time UFR with Parameters: " << std::endl;
    std::cout << "PRF: " << prf << std::endl;
    std::cout << "DC Rate: " << dc_rate << std::endl;
    std::cout << "Burst Duration: " << burst_duration << std::endl;
    std::cout << "Doppler Bandwidth: " << processing_bandwidth << std::endl;
    std::cout << "Ramp Signal Bandwidth: " << bandwidth << std::endl;
    std::cout << "Focused Time: " << focused_time << std::endl;
    std::cout << "Mean Azimuth FM Rate: " << ka_mean << std::endl;
    std::cout << "Number of Freq Spectral Replicas: " << num_replicas << std::endl;
    std::cout << "Number of Time Spectral Replicas: " << num_tiles << std::endl;
    std::cout << "Tiled Shape: " << shape << std::endl;
    std::cout << "Low Pass Filter Center: " << mid << std::endl;
    std::cout << "Low Pass Filter Start Index: " << start << std::endl;
    std::cout << "Low Pass Filter End Index: " << end << std::endl;
    std::cout << "Downsample Shape: " << downsample_shape << std::endl;
    std::cout << "Output Shape: " << output_shape << std::endl;
    std::cout << "-------------------------------------------------" << std::endl;

    F_VEC_1D az_times = linspace(-focused_time / 2, focused_time / 2, shape);

    CF_VEC_2D ufr_intermediate(num_rng, CF_VEC_1D(shape));

    #pragma omp parallel for
    for (int rng_line = 0; rng_line < num_rng; rng_line++)
    {
        // Mosaic Signal
        ufr_intermediate[rng_line] = get_tiled_signal(azimuth_compressed[rng_line], num_tiles);

        //Apply De-ramping
        double f_dc = dc_estimates[rng_line];
        F_VEC_1D ka_row = linear_resample(ka[rng_line], shape);

        for (int az_line = 0; az_line < shape; az_line++)
        {
            double kt = -ka_row[az_line] * dc_rate / (dc_rate - ka_row[az_line]);
            double dc_time = -f_dc / ka_row[az_line];
            double az_time = dc_time + az_times[az_line];

            std::complex<double> deramp = std::exp( (-1.0 * I * PI * kt * std::pow(az_time, 2.0)) + (2.0 * I * PI * (kt + ka_row[az_line]) * dc_time * az_time) );
            std::complex<double> sample = ufr_intermediate[rng_line][az_line];

            ufr_intermediate[rng_line][az_line] = deramp * sample;
        }
    }

    // Low-pass Filter
    compute_axis_dft_in_place(ufr_intermediate, 0, 1, false);

    #pragma omp parallel for collapse(2)
    for (int rng_line = 0; rng_line < num_rng; rng_line++)
    {
        for (int az_line = 0; az_line < shape; az_line++)
        {
            if (az_line < end-start)
            {
                ufr_intermediate[rng_line][az_line] = ufr_intermediate[rng_line][start + az_line];
            }
            else
            {
                ufr_intermediate[rng_line][az_line] = 0.0;
            }
        }
    }

    compute_axis_dft_in_place(ufr_intermediate, 0, 1, true);

    // Resample
    CF_VEC_2D ufr_output(num_rng, CF_VEC_1D(downsample_shape));

    std::transform(
        ufr_intermediate.begin(), ufr_intermediate.end(),
            ufr_output.begin(),
                [downsample_shape] (CF_VEC_1D& ufr_row) {
                    return linear_resample(ufr_row, downsample_shape);
                }
    );

    ufr_intermediate.clear();
    ufr_intermediate.shrink_to_fit();

    // Apply Re-ramping
    #pragma omp parallel for
    for (int rng_line = 0; rng_line < num_rng; rng_line++)
    {
        double f_dc = dc_estimates[rng_line];
        F_VEC_1D ka_row = linear_resample(ka[rng_line], downsample_shape);

        for (int az_line = 0; az_line < downsample_shape; az_line++)
        {
            double kt = -ka_row[az_line] * dc_rate / (dc_rate - ka_row[az_line]);
            double dc_time = -f_dc / ka_row[az_line];
            double az_time = dc_time + az_times[az_line];

            std::complex<double> reramp = std::exp( (1.0 * I * PI * kt * std::pow(az_time, 2.0)) - (2.0 * I * PI * (kt + ka_row[az_line]) * dc_time * az_time) );
            std::complex<double> sample = ufr_output[rng_line][az_line];

            ufr_output[rng_line][az_line] = reramp * sample;
        }
    }

    ufr_output = transpose(ufr_output);

    // Crop to Valid Portion
    mid = downsample_shape / 2;
    start = mid - (output_shape / 2);
    end = mid + (output_shape / 2);

    return CF_VEC_2D(ufr_output.begin() + start, ufr_output.begin() + end);
}