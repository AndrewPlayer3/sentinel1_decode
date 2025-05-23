#include "image_formation.h"


CF_VEC_1D get_reference_function(const CF_VEC_1D& replica_chirp)
{
    int num_samples = replica_chirp.size();

    CF_VEC_1D reference = replica_chirp;

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

    CF_VEC_1D signal_ = signal;

    std::transform(
        reference.begin(), reference.end(),
            signal_.begin(), signal_.begin(),
                [] (const std::complex<double>& n, std::complex<double>& r) { return n * r;}
    );

    return CF_VEC_1D(signal_.begin(), signal_.begin() + num_samples);
}


std::pair<PACKET_VEC_2D, int> get_azimuth_blocks(PACKET_VEC_1D& packets)
{
    PACKET_VEC_2D azimuth_blocks;
    PACKET_VEC_1D azimuth_block;
    
    int previous_size = 2 * packets[0].get_num_quads();
    int max_size = previous_size;

    for (int i = 0; i < packets.size(); i++)
    {
        L0Packet packet = packets[i];
        int size = 2 * packet.get_num_quads();
        if (size != previous_size or i == packets.size() - 1)
        {
            if (size > max_size) max_size = size;
            previous_size = size;
            azimuth_blocks.push_back(azimuth_block);
            azimuth_block = {};
        }
        azimuth_block.push_back(packet);
    }

    return std::pair<PACKET_VEC_2D, int>(azimuth_blocks, max_size);
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
    const double& doppler_bandwidth  // B_d
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
    int start = mid - (std::floor(doppler_bandwidth) / 2);
    int end = mid + (std::floor(doppler_bandwidth) / 2);

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

    F_VEC_1D  az_freqs = linspace(-bandwidth, bandwidth, shape);
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
            std::complex<double> deramp = 
                std::exp(-I * PI * (1 / dc_rate) * std::pow(az_freqs[az_line], 2.0));

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
            std::complex<double> reramp = 
                std::exp(I * PI * (1 / dc_rate) * std::pow(az_freqs[az_line], 2.0));

            ufr_output[rng_line][az_line] = reramp * sample;
        }
    }

    return transpose(ufr_output);
}


CF_VEC_2D azimuth_time_ufr(
    CF_VEC_2D& azimuth_compressed,
    F_VEC_1D&  dc_estimates,
    F_VEC_2D& az_fm_rate,
    L0Packet&  initial_packet,
    const double& dc_rate,
    const double& burst_duration,
    const double& prf,
    const double& doppler_bandwidth,  // B_d
    const int&    swath_number
) {
    az_fm_rate = transpose(az_fm_rate);
    azimuth_compressed = transpose(azimuth_compressed);

    int num_az = azimuth_compressed[0].size();
    int num_rng = azimuth_compressed.size();

    double num_replicas = std::abs(dc_rate * burst_duration / prf);
    double bandwidth = (num_replicas * prf) - doppler_bandwidth;

    double ka_mean = std::abs(std::accumulate(az_fm_rate[0].begin(), az_fm_rate[0].end(), 0.0)) / num_az;
    double focused_time = burst_duration + ((bandwidth - (2 * doppler_bandwidth)) / ka_mean);

    int num_pos_tiles = std::ceil((focused_time / 2) / burst_duration);
    int num_neg_tiles = std::floor((-focused_time / 2) / burst_duration);
    int num_tiles = num_pos_tiles - num_neg_tiles;

    int output_extension = swath_number == 11 ? 
                             std::floor(doppler_bandwidth * 1.15) 
                             : 
                             std::floor(doppler_bandwidth * 0.26);
    int output_shape = std::floor(prf * burst_duration);
    int downsample_shape = output_shape + output_extension;
    int shape = num_tiles * num_az;
    int mid = shape / 2;
    int start = mid - bandwidth / 2;
    int end = mid + bandwidth / 2;

    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << "Computing Azimuth Time UFR with Parameters: " << std::endl;
    std::cout << "PRF: " << prf << std::endl;
    std::cout << "DC Rate: " << dc_rate << std::endl;
    std::cout << "Burst Duration: " << burst_duration << std::endl;
    std::cout << "Doppler Bandwidth: " << doppler_bandwidth << std::endl;
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
        F_VEC_1D ka = linear_resample(az_fm_rate[rng_line], shape);

        // Mosaic Signal
        ufr_intermediate[rng_line] = get_tiled_signal(azimuth_compressed[rng_line], num_tiles);

        //Apply De-ramping
        double f_dc = dc_estimates[rng_line];

        for (int az_line = 0; az_line < shape; az_line++)
        {
            double kt = -ka[az_line] * dc_rate / (dc_rate - ka[az_line]);
            double dc_time = -f_dc / ka[az_line];
            double az_time = dc_time + az_times[az_line];

            std::complex<double> deramp_phase_1 = -1.0 * I * PI * kt * std::pow(az_time, 2.0);
            std::complex<double> deramp_phase_2 = 2.0 * I * PI * (kt + ka[az_line]) * dc_time * az_time;
            std::complex<double> deramp = std::exp(deramp_phase_1 + deramp_phase_2);
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
                    return quadratic_resample(ufr_row, downsample_shape);
                }
    );

    ufr_intermediate.clear();
    ufr_intermediate.shrink_to_fit();

    // Apply Re-ramping
    #pragma omp parallel for
    for (int rng_line = 0; rng_line < num_rng; rng_line++)
    {
        F_VEC_1D ka = linear_resample(az_fm_rate[rng_line], downsample_shape);

        double f_dc = dc_estimates[rng_line];

        for (int az_line = 0; az_line < downsample_shape; az_line++)
        {
            double kt = -ka[az_line] * dc_rate / (dc_rate - ka[az_line]);
            double dc_time = -f_dc / ka[az_line];
            double az_time = dc_time + az_times[az_line];

            std::complex<double> reramp_phase_1 = 1.0 * I * PI * kt * std::pow(az_time, 2.0);
            std::complex<double> reramp_phase_2 = 2.0 * I * PI * (kt + ka[az_line]) * dc_time * az_time;
            std::complex<double> reramp = std::exp(reramp_phase_1 - reramp_phase_2);
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
