#include "doppler.h"


F_VEC_1D get_wrapped_estimates(CF_VEC_2D& range_compressed, const double& prf)
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

    F_VEC_1D wrapped_estimates(num_rng);

    std::transform(
        phase_difference_sum.begin(), phase_difference_sum.end(),
            wrapped_estimates.begin(),
                [prf] (const std::complex<double>& v) {
                    return -prf * std::arg(v) / (2 * PI);
                }
    );

    return wrapped_estimates;
}


F_VEC_1D get_unwrapped_estimates(F_VEC_1D& fine_dc_estimates, F_VEC_1D& offset_rng_times, const double& prf, const int& dense_grid_scalar = 8)
{
    int num_est = fine_dc_estimates.size();
    float dt = offset_rng_times[1] - offset_rng_times[0];

    CF_VEC_1D F(num_est * dense_grid_scalar);
    std::transform(
        fine_dc_estimates.begin(), fine_dc_estimates.end(),
            F.begin(),
                [prf] (const double& dc_est) {
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

    F_VEC_1D unwrapped_estimates(num_est);
    std::transform(
        fine_dc_estimates.begin(), fine_dc_estimates.end(),
            offset_rng_times.begin(), unwrapped_estimates.begin(),
                [prf, a, b] (double& dc_est, double& time) {
                    std::complex<double> res = std::exp(2.0 * I * PI * dc_est / prf) * std::exp(-1.0 * I * (a * time + b));
                    double r = std::arg(res) / (2.0 * PI);
                    return (a * time + b + r) * prf;
                }
    );

    return unwrapped_estimates;
}


F_VEC_1D polyfit(const F_VEC_1D& x, const F_VEC_1D& y) {
    const int N = x.size();

    F_VEC_1D Sx(7);
    F_VEC_1D Sy(4);

    for (int i = 0; i < N; ++i) {
        double xi = 1.0;
        for (int j = 0; j <= 6; ++j) {
            if (j <= 3) Sy[j] += xi * y[i];
            Sx[j] += xi;
            xi *= x[i];
        }
    }

    // Construct the normal equations matrix and RHS
    F_VEC_2D A = {
        {Sx[0], Sx[1], Sx[2], Sx[3]},
        {Sx[1], Sx[2], Sx[3], Sx[4]},
        {Sx[2], Sx[3], Sx[4], Sx[5]},
        {Sx[3], Sx[4], Sx[5], Sx[6]}
    };
    F_VEC_1D b = {Sy[0], Sy[1], Sy[2], Sy[3]};

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
    F_VEC_1D coeffs(4);
    for (int i = 3; i >= 0; --i) {
        coeffs[i] = b[i];
        for (int j = i + 1; j < 4; ++j)
            coeffs[i] -= A[i][j] * coeffs[j];
        coeffs[i] /= A[i][i];
    }

    return coeffs; // [d, c, b, a]
}


double polyval(const F_VEC_1D& coeffs, double x) {
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

    for (int rng_line = 0; rng_line < num_rng; rng_line++)
    {
        std::transform(
            deramping_signal.begin(), deramping_signal.end(),
                preconditioned[rng_line].begin(), preconditioned[rng_line].begin(),
                    [] (std::complex<double> deramp, std::complex<double> rc) {
                        return rc * deramp;
                    }
        );
    }

    return transpose(preconditioned);
}


F_VEC_1D get_doppler_centroid(CF_VEC_2D& range_compressed, const double& doppler_centroid_rate, const double& burst_duration, L0Packet& first_packet)
{
    int num_az = range_compressed.size();
    int num_rng = range_compressed[0].size();

    double prf = 1.0 / (first_packet.get_pri() * 1e-6);

    CF_VEC_2D precoditioned_burst = dce_preconditioning(range_compressed, doppler_centroid_rate, 1.06);

    F_VEC_1D wrapped_estimates = get_wrapped_estimates(precoditioned_burst, prf);

    // TODO: Make this an argument once things are working
    int num_rng_blocks = 15;
    int rng_block_size = num_rng / num_rng_blocks;

    F_VEC_1D fine_dc_estimates(num_rng_blocks);

    for (int rng_block_index = 0; rng_block_index < num_rng_blocks; rng_block_index++)
    {
        int start_index = rng_block_index * rng_block_size;

        int end_index = start_index + rng_block_size < num_rng ?
                            start_index + rng_block_size 
                            :
                            num_rng - 1;

        F_VEC_1D rng_block(wrapped_estimates.begin() + start_index, wrapped_estimates.begin() + end_index);

        fine_dc_estimates[rng_block_index] = std::accumulate(rng_block.begin(), rng_block.end(), 0.0) / (double(end_index) - double(start_index));
    }

    std::cout << "Wrapped F_DC: " << fine_dc_estimates.front() << ", " << fine_dc_estimates.back() << std::endl;

    F_VEC_1D rng_times = first_packet.get_slant_range_times(num_rng_blocks);
    double initial_time = rng_times[0];
    std::transform(
        rng_times.begin(), rng_times.end(),
            rng_times.begin(),
                [initial_time] (const double& t) {
                    return t - initial_time;
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
                [initial_time] (const double& t) {
                    return t - initial_time;
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
