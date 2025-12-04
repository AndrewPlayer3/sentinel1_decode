#include "s1_decoder.h"


void S1_Decoder::_set_packets()
{
    _flat_packets = L0Packet::get_packets(_filename, 0);

    for (L0Packet packet : _flat_packets) 
    {
        _swath_counts[packet.get_swath()]++;
    }

    for (std::pair<std::string, int> swath_count : _swath_counts)
    {
        std::string name = swath_count.first;
        if (ECHO_SWATHS.contains(name))  
        {
            _echo_packets[name] = L0Packet::get_packets_in_bursts(_flat_packets, name);
        }
        else if (CAL_SWATHS.contains(name))
        {
            _cal_packets[name] = L0Packet::get_packets_in_bursts(_flat_packets, name, true);
        }
    }
}


bool is_sm(const std::string& swath)
{
    return STRIPMAP_SWATHS.contains(swath);
}


bool is_iw(const std::string& swath)
{
    return IW_SWATHS.contains(swath);
}


bool is_ew(const std::string& swath)
{
    return EW_SWATHS.contains(swath);
}


bool is_wv(const std::string& swath)
{
    return WV_SWATHS.contains(swath);
}

bool is_cal(const std::string& swath)
{
    return CAL_SWATHS.contains(swath);
}

STATE_VECTORS S1_Decoder::get_state_vectors()
{
    return _state_vectors;
}

void S1_Decoder::apply_eccm(const int& detection_threshold, const int& mitigation_threshold)
{
    _use_eccm = true;
    _eccm_detection_threshold = detection_threshold;
    _eccm_mitigation_threshold = mitigation_threshold;
}

void S1_Decoder::_validate_request(const std::string& swath, const int& burst)
{
    try
    {
        _swath_counts.at(swath);
    }
    catch(const std::exception& e)
    {
        std::cout << "The requested swath, " << swath << ", does not exist in the data." 
                  << " The available swaths are: " << std::endl;

        for (std::pair<std::string, int> swath_count : _swath_counts)
        {
            std::cout << swath_count.first << std::endl;
        }
    
        std::exit(1);
    }

    if (burst)
    {
        int num_bursts = 0;
        try
        {
            if (is_cal(swath))
            {
                num_bursts = _cal_packets.at(swath).size();
                _cal_packets.at(swath).at(burst);
            }
            else
            {
                num_bursts = _echo_packets.at(swath).size();
                _echo_packets.at(swath).at(burst);
            }
        }
        catch(const std::exception& e)
        {
            std::cout << "The requested burst, " << burst << ", does not exist in the data." 
                    << " There are "<< num_bursts << " bursts in swath " << swath << "."
                    << std::endl;

            std::exit(1);
        }
    }
}


CF_VEC_2D S1_Decoder::get_burst(
    const std::string& swath,
    const int& burst
) {
    _validate_request(swath, burst);
    
    PACKET_VEC_1D burst_packets = _echo_packets[swath][burst];

    int num_packets = burst_packets.size();
    int num_samples = 2 * burst_packets[0].get_num_quads();

    CF_VEC_2D signals(num_packets, CF_VEC_1D(num_samples));

    #pragma omp parallel for 
    for (int i = 0; i < num_packets; i++)
    {
        L0Packet packet = burst_packets[i];
        signals[i] = packet.get_signal();
    }

    return signals;
}


CF_VEC_2D S1_Decoder::_get_cal_swath(const std::string& swath)
{
    PACKET_VEC_2D swath_packets = _cal_packets[swath];

    int num_packets = swath_packets[0].size();
    int num_samples = 2 * swath_packets[0][0].get_num_quads();

    CF_VEC_2D signals(num_packets, CF_VEC_1D(num_samples));

    #pragma omp parallel for
    for (int i = 0; i < num_packets; i++)
    {
        L0Packet packet = swath_packets[0][i];
        CF_VEC_1D signal = packet.get_signal();
        signals[i] = signal;
    }

    return signals;
}


CF_VEC_2D S1_Decoder::get_swath(const std::string& swath)
{
    _validate_request(swath);

    PACKET_VEC_2D swath_packets;

    if (is_iw(swath)) swath_packets = _echo_packets[swath];
    else if (is_sm(swath)) swath_packets = get_azimuth_blocks(_echo_packets[swath][0]).first;
    else if (is_cal(swath)) return _get_cal_swath(swath);
    else throw(std::domain_error("Only IW and SM modes are currently supported."));

    int num_bursts = swath_packets.size();
    int num_packets_total = 0;
    int max_samples = 0;

    INT_VEC_1D first_packet_index(num_bursts);

    for (int i = 0; i < num_bursts; i++)
    {
        int num_packets = swath_packets[i].size();
        num_packets_total += num_packets;
        L0Packet packet = swath_packets[i][0];
        int num_samples = 2 * packet.get_num_quads();
        first_packet_index[i] = i == 0 ? 0 : first_packet_index[i-1] + swath_packets[i-1].size();
        if (num_samples > max_samples) max_samples = num_samples;
    }

    std::cout << "Number of Bursts: "  << num_bursts  << "\n"
              << "Number of Samples: " << max_samples << std::endl;

    CF_VEC_2D signals(num_packets_total, CF_VEC_1D(max_samples));

    for (int i = 0; i < num_bursts; i++)
    {
        std::cout << "Decoding Burst #" << i << " of " << num_bursts << "." << std::endl;

        PACKET_VEC_1D& burst = swath_packets[i];

        int num_packets = burst.size();

        #pragma omp parallel for
        for (int j = 0; j < num_packets; j++)
        {
            L0Packet packet = burst[j];
            CF_VEC_1D signal = packet.get_signal();
            int index = first_packet_index[i] + j;
            signals[index] = signal;
        }
    }
    return signals;
}


CF_VEC_2D S1_Decoder::get_range_compressed_burst(
    const std::string& swath,
    const int& burst,
    const bool& range_doppler
) {
    _validate_request(swath, burst);

    PACKET_VEC_1D burst_packets = _echo_packets[swath][burst];

    return _range_compress(burst_packets, true, range_doppler);
}


CF_VEC_2D S1_Decoder::_get_range_compressed_swath_sm(
    const std::string& swath,
    const bool& range_doppler
) {    
    PACKET_VEC_1D burst_packets = _echo_packets[swath][0];

    std::pair<PACKET_VEC_2D, int> azimuth_block_pair = get_azimuth_blocks(burst_packets);
    int max_size = azimuth_block_pair.second;
    PACKET_VEC_2D azimuth_blocks = std::move(azimuth_block_pair.first);

    std::cout << "Number of Azimuth Blocks: " << azimuth_blocks.size() << std::endl;

    if (azimuth_blocks.size() == 1)
    {
        return _range_compress(azimuth_blocks[0], true, range_doppler);
    }

    CF_VEC_2D range_compressed;

    for (unsigned int i = 0; i < azimuth_blocks.size(); i++)
    {        
        CF_VEC_2D range_compressed_az_block = _range_compress(azimuth_blocks[i], true, range_doppler);

        int rows = range_compressed_az_block.size();
        int size = range_compressed_az_block[0].size();
        int diff = max_size - size;

        for (int i = 0; i < rows; i++)
        {
            CF_VEC_1D row = range_compressed_az_block.front();
            range_compressed_az_block.erase(range_compressed_az_block.begin());
            for (int j = 0; j < diff; j++)
            {
                row.emplace(row.begin(), std::complex<double>(0.0));
            }
            range_compressed.push_back(row);
        }
    }

    return range_compressed;
}


CF_VEC_2D S1_Decoder::_get_range_compressed_swath_iw(
    const std::string& swath,
    const bool& range_doppler
) {
    PACKET_VEC_2D packets = _echo_packets[swath];

    int num_bursts = packets.size();

    CF_VEC_2D range_compressed;

    for (int i = 0; i < num_bursts; i++)
    {
        std::cout << "Range Compressing Burst #" << i << std::endl;
        CF_VEC_2D range_compressed_burst = _range_compress(packets[i], true, range_doppler);
        int num_azimuth_lines = range_compressed_burst.size();
        for (int j = 0; j < num_azimuth_lines; j++)
        {
            CF_VEC_1D row = range_compressed_burst.front();
            range_compressed_burst.erase(range_compressed_burst.begin());
            range_compressed.push_back(row);
        }
    }

    return range_compressed;
}


CF_VEC_2D S1_Decoder::get_range_compressed_swath(
    const std::string& swath,
    const bool& range_doppler
) {
    _validate_request(swath);    

    if (is_sm(swath))
    {
        return _get_range_compressed_swath_sm(swath, range_doppler);
    }
    else if (is_iw(swath))
    {
        return _get_range_compressed_swath_iw(swath, range_doppler);
    }
    else
    {
        std::string msg = 
            swath + " is not a supported swath type. Only IW and SM imaging modes are supported.";
        throw(std::invalid_argument(msg));
    }
}


CF_VEC_2D S1_Decoder::get_azimuth_compressed_burst(
    const std::string& swath,
    const int& burst
) {
    _validate_request(swath, burst);    

    PACKET_VEC_1D packets = _echo_packets[swath][burst];
    return _azimuth_compress(packets, true);
}


CF_VEC_2D S1_Decoder::_get_azimuth_compressed_swath_iw(const std::string& swath)
{
    PACKET_VEC_2D packets = _echo_packets[swath];

    int num_bursts = packets.size();

    CF_VEC_2D azimuth_compressed;

    for (int i = 0; i < num_bursts; i++)
    {
        std::cout << "Azimuth Compressing Burst #" << i 
                  << std::endl
                  << "Burst #" << i << " contains " << packets[i].size() << " packets." 
                  << std::endl; 

        if (packets[i].size() < 1000)
        {
            std::cout << "Burst #" << i 
                      << " does not contain enough packets for image formation, skipping..." 
                      << std::endl;
            continue;
        }

        CF_VEC_2D azimuth_compressed_burst = _azimuth_compress(packets[i], true);

        int num_lines = azimuth_compressed_burst.size();

        for (int j = 0; j < num_lines; j++)
        {
            CF_VEC_1D row = azimuth_compressed_burst.front();
            azimuth_compressed_burst.erase(azimuth_compressed_burst.begin());
            azimuth_compressed.push_back(row);
        }

        CF_VEC_1D zero_padding(azimuth_compressed[0].size());
        for (int j = 0; j < 20; j++)
        {
            azimuth_compressed.push_back(zero_padding);
        }
    }

    return azimuth_compressed;
}


CF_VEC_2D S1_Decoder::_get_azimuth_compressed_swath_sm(const std::string& swath)
{
    PACKET_VEC_1D burst_packets = _echo_packets[swath][0];

    std::pair<PACKET_VEC_2D, int> azimuth_block_pair = get_azimuth_blocks(burst_packets);
    int max_size = azimuth_block_pair.second;
    PACKET_VEC_2D azimuth_blocks = std::move(azimuth_block_pair.first);

    std::cout << "Number of Azimuth Blocks: " << azimuth_blocks.size() << std::endl;

    if (azimuth_blocks.size() == 1)
    {
        return _azimuth_compress(azimuth_blocks[0]);
    }

    CF_VEC_2D azimuth_compressed;

    for (unsigned int i = 0; i < azimuth_blocks.size(); i++)
    {
        CF_VEC_2D azimuth_compressed_az_block = _azimuth_compress(azimuth_blocks[i]);

        int rows = azimuth_compressed_az_block.size();
        int size = azimuth_compressed_az_block[0].size();
        int diff = max_size - size;

        for (int j = 0; j < rows; j++)
        {
            CF_VEC_1D row = azimuth_compressed_az_block.front();
            azimuth_compressed_az_block.erase(azimuth_compressed_az_block.begin());
            for (int j = 0; j < diff; j++)
            {
                row.emplace(row.begin(), std::complex<double>(0.0));
            }
            azimuth_compressed.push_back(row);
        }

        CF_VEC_1D zero_padding(azimuth_compressed[0].size());
        for (int j = 0; j < 20; j++)
        {
            azimuth_compressed.push_back(zero_padding);
        }
    }

    return azimuth_compressed;
}


CF_VEC_2D S1_Decoder::get_azimuth_compressed_swath(const std::string& swath)
{
    _validate_request(swath);    

    if (is_sm(swath))
    {
        return _get_azimuth_compressed_swath_sm(swath);
    }
    else if (is_iw(swath))
    {
        return _get_azimuth_compressed_swath_iw(swath);
    }
    else
    {
        std::string msg = 
            swath + " is not a supported swath type. Only IW and SM imaging modes are supported.";
        throw(std::invalid_argument(msg));
    }
}


CF_VEC_2D S1_Decoder::_range_compress(
    PACKET_VEC_1D& packets,
    const bool& do_ifft,
    const bool& do_azimuth_fft
) {
    L0Packet first_packet = packets[0];

    int num_packets = packets.size();
    int num_samples = 2 * packets[0].get_num_quads();

    CF_VEC_2D range_compressed(num_packets, CF_VEC_1D(num_samples));

    std::cout << "Reading Complex Data" << std::endl;

    #pragma omp parallel for
    for (int i = 0; i < num_packets; i++)
    {
        L0Packet packet = packets[i];
        CF_VEC_1D signal = packet.get_signal();
        range_compressed[i] = signal;
    }

    if (_use_eccm)
    {
        char rx_pol = first_packet.get_rx_polarization();
        eccm(range_compressed, 128, 64, _eccm_detection_threshold, _eccm_mitigation_threshold, rx_pol);
    }

    compute_axis_dft_in_place(range_compressed, 0, 1, false);

    std::cout << "Range Compressing" << std::endl;

    CF_VEC_1D reference_function = get_reference_function(packets[0].get_replica_chirp());

    for (int i = 0; i < num_packets; i++)
    {
        range_compressed[i] = pulse_compression(range_compressed[i], reference_function);
    }

    if (do_ifft)
    {
        compute_axis_dft_in_place(range_compressed, 0, 1, true);
        fftshift_in_place(range_compressed);
    }

    std::for_each(
        range_compressed.begin(), range_compressed.end(),
            [num_samples] (CF_VEC_1D& row) { 
                std::rotate(row.begin(), row.end()-(num_samples / 2), row.end());
            }
    );

    if (do_azimuth_fft)
    {
        range_compressed = transpose(range_compressed);
        compute_axis_dft_in_place(range_compressed, 0, 1, false);
        fftshift_in_place(range_compressed);
        range_compressed = transpose(range_compressed);
    }

    return range_compressed;
}


CF_VEC_2D S1_Decoder::_azimuth_compress(PACKET_VEC_1D& packets, const bool& tops_mode)
{
    L0Packet first_packet = packets[0];

    int num_packets = packets.size();
    int num_samples = 2 * packets[0].get_num_quads();

    CF_VEC_2D radar_data;

    if (tops_mode) radar_data = _range_compress(packets, true, false);
    else radar_data = _range_compress(packets, false, true);

    D_VEC_1D slant_ranges = packets[0].get_slant_ranges();

    std::cout << "Slant Ranges: "
              << slant_ranges[0] << ", " << slant_ranges.back() << std::endl;

    D_VEC_1D v_0 = _state_vectors.velocities[0];

    int swath_number = packets[0].secondary_header("swath_number");
    double v_norm = std::sqrt(std::pow(v_0[0], 2.0) + std::pow(v_0[1], 2.0) + std::pow(v_0[2], 2.0));
    double pri = packets[0].get_pri() * 1e-6;
    double prf = 1 / pri;
    double burst_length_seconds = double(num_packets) / prf;
    double doppler_bandwidth = swath_number == 11 ? 675 : 600;
    double t0 = first_packet.get_time();
    double time_delta = burst_length_seconds / prf;
    double range_dec_sample_rate = first_packet.get_range_sample_rate() * 1e+6;
    double dc_rate;

    D_VEC_1D slant_range_times = first_packet.get_slant_range_times();

    std::cout << range_dec_sample_rate << std::endl;
    std::cout << slant_range_times.front() << " " << slant_range_times.back() << std::endl;

    D_VEC_1D range_freqs = linspace(-range_dec_sample_rate / 2, range_dec_sample_rate / 2, num_samples);
    D_VEC_1D az_freqs(num_packets);
    D_VEC_1D doppler_centroid(num_samples);
    D_VEC_2D az_fm_rate;

    if (tops_mode)
    {
        dc_rate = get_doppler_centroid_rate(packets, v_norm);

        doppler_centroid = get_doppler_centroid(
            radar_data,
            dc_rate,
            burst_length_seconds,
            first_packet,
            30
        );

        radar_data = azimuth_frequency_ufr(
            radar_data,
            doppler_centroid,
            first_packet,
            dc_rate,
            burst_length_seconds,
            prf,
            doppler_bandwidth
        );

        num_packets = radar_data.size();
        time_delta /= double(num_packets) / double(packets.size());
    }

    D_VEC_2D positions(num_packets, D_VEC_1D(3));
    D_VEC_1D velocities_norm(num_packets);

    for (int i = 0; i < num_packets; i++)
    {
        double t = t0 + time_delta * i;

        STATE_VECTOR state_vector = _state_vectors.interpolate(t);

        D_VEC_1D v = state_vector.velocity;
        D_VEC_1D p = state_vector.position;

        positions[i]  = p;
        velocities_norm[i] = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    }

    std::vector<fftw_plan> forward_plans = get_fftw_plans(radar_data);
    std::vector<fftw_plan> inverse_plans = get_fftw_plans(radar_data, FFTW_BACKWARD);

    if (tops_mode)
    {
        double num_replicas = std::abs(dc_rate * burst_length_seconds / prf);
        az_fm_rate = D_VEC_2D(num_packets, D_VEC_1D(num_samples));
        az_freqs = linspace(-num_replicas*prf/2, num_replicas*prf/2, num_packets);
    }
    else
    {
        az_freqs = linspace(-prf/2, prf/2, num_packets);
    }

    std::cout << "Azimuth Compressing" << std::endl;

    D_VEC_1D time_corrections = first_packet.get_timing_corrections();

    CF_VEC_1D az_filter(num_packets);
    CF_VEC_1D filter_energies(num_samples);

    #pragma omp parallel for
    for (int i = 0; i < num_packets; i++)
    {
        if (tops_mode)
        {
            fftw_execute(forward_plans[i]);
        }

        CF_VEC_1D& range_line = radar_data[i];

        D_VEC_1D effective_velocities = get_effective_velocities(
            positions[i],
            velocities_norm[i],
            slant_ranges
        );

        D_VEC_1D rcmc_factors = apply_secondary_range_compression(
            range_line,
            effective_velocities,
            slant_ranges,
            range_freqs,
            doppler_centroid,
            az_freqs[i]
        );

        fftw_execute(inverse_plans[i]);
        fftshift_in_place(range_line);

        if (tops_mode)
        {
            apply_range_cell_migration_correction(
                range_line,
                slant_ranges,
                range_freqs,
                rcmc_factors
            );
        }

        for (int j = 0; j < num_samples; j++)
        {
            double v_rel = effective_velocities[j];
            double rcmc_factor = rcmc_factors[j];
            double slant_range = slant_ranges[j];

            std::complex<double> time_correction = std::exp( 2.0 * I * PI * (az_freqs[i] + doppler_centroid[j]) * time_corrections[j] );

            std::complex<double> az_filter = 
                std::exp(4.0 * I * PI * slant_range * rcmc_factor / WAVELENGTH) * time_correction;

            // Unnecessary until windowing is added
            filter_energies[j] += std::pow( std::abs(az_filter), 2.0 );

            if (tops_mode) 
            {
                az_fm_rate[i][j] = 
                    -2 * std::pow(v_rel, 2.0) * std::pow(rcmc_factors[j], 3.0) / (WAVELENGTH * slant_range);
            }

            range_line[j] *= az_filter;
        }
    }

    destroy_fftw_plans(forward_plans);
    destroy_fftw_plans(inverse_plans);

    compute_axis_dft_in_place(radar_data, 0, 0, true);

    if (tops_mode)
    {
        radar_data = azimuth_time_ufr(
            radar_data,
            doppler_centroid,
            az_fm_rate,
            first_packet,
            dc_rate,
            burst_length_seconds,
            prf,
            doppler_bandwidth,
            swath_number
        );
        num_packets = radar_data.size();
    }

    return radar_data;
}
