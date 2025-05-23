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


CF_VEC_2D S1_Decoder::get_swath(const std::string& swath)
{
    _validate_request(swath);

    PACKET_VEC_2D swath_packets;

    if (is_iw(swath)) swath_packets = _echo_packets[swath];
    else if (is_sm(swath)) swath_packets = get_azimuth_blocks(_echo_packets[swath][0]).first;
    else if (is_cal(swath)) swath_packets = _cal_packets[swath];
    else throw(std::domain_error("Only IW and SM modes are currently supported."));

    int num_bursts = swath_packets.size();
    int num_packets_total = 0;
    INT_VEC_1D first_packet_index(num_bursts);
    int max_samples = 0;

    if (is_cal(swath))
    {
        int num_packets = swath_packets.size();
        int num_samples = 2 * swath_packets[0][0].get_num_quads();
        int signals_index = 0;
        CF_VEC_2D signals(num_packets, CF_VEC_1D(num_samples));

        #pragma omp parallel for
        for (int i = 0; i < num_packets; i++)
        {
            L0Packet packet = swath_packets[i][0];
            CF_VEC_1D signal = packet.get_signal();
            signals[i] = signal;
        }

        return signals;
    }

    std::cout << "Initializing Values" << std::endl;

    for (int i = 0; i < num_bursts; i++)
    {
        int num_packets = swath_packets[i].size();
        num_packets_total += num_packets;
        L0Packet packet = swath_packets[i][0];
        int num_samples = 2 * packet.get_num_quads();
        first_packet_index[i] = i == 0 ? 0 : first_packet_index[i-1] + swath_packets[i-1].size();
        if (num_samples > max_samples) max_samples = num_samples;
    }

    std::cout << "Number of Bursts: " << num_bursts << std::endl;
    std::cout << "Number of Samples: " << max_samples << std::endl;

    int signals_index = 0;
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

    for (int i = 0; i < azimuth_blocks.size(); i++)
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
    bool range_doppler
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

    for (int i = 0; i < azimuth_blocks.size(); i++)
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

    std::cout << "Converting the Signals into the Range/Azimuth Domain" << std::endl;

    compute_axis_dft_in_place(range_compressed, 0, 1, false);

    std::cout << "Range Compressing" << std::endl;

    CF_VEC_1D reference_function = get_reference_function(packets[0].get_replica_chirp());

    for (int i = 0; i < num_packets; i++)
    {
        range_compressed[i] = pulse_compression(range_compressed[i], reference_function);
    }

    std::cout << "Inverse FFT" << std::endl;

    if (do_ifft)
    {
        compute_axis_dft_in_place(range_compressed, 0, 1, true);
    }

    std::cout << "FFT Shift" << std::endl;

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

    F_VEC_1D v_0 = _state_vectors.velocities[0];

    double v_norm = std::sqrt(std::pow(v_0[0], 2.0) + std::pow(v_0[1], 2.0) + std::pow(v_0[2], 2.0));
    double range_sample_rate = packets[0].get_range_sample_rate();
    double pulse_length = packets[0].get_pulse_length() * 1e-6;
    double pri = packets[0].get_pri() * 1e-6;
    double prf = 1 / pri;
    double burst_length_seconds = double(num_packets) / prf;
    double dc_rate = get_doppler_centroid_rate(packets, v_norm);
    double doppler_bandwidth = prf * 0.4;

    double oversample_factor;
    double t0;
    double t1;
    double time_delta;

    F_VEC_1D doppler_centroid;

    if (tops_mode)
    {
        doppler_centroid = get_doppler_centroid(
            radar_data,
            dc_rate,
            burst_length_seconds,
            first_packet
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

        oversample_factor = double(num_packets) / double(packets.size());
        t0 = packets[0].get_time();
        t1 = packets[1].get_time();
        time_delta = (t1 - t0) / oversample_factor;
    }
    else
    {
        t0 = packets[0].get_time();
        t1 = packets[1].get_time();
        time_delta = (t1 - t0);
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

    D_VEC_1D range_freqs = linspace(-range_sample_rate/2, range_sample_rate/2, num_samples);

    std::vector<fftw_plan> forward_plans = get_fftw_plans(radar_data);
    std::vector<fftw_plan> inverse_plans = get_fftw_plans(radar_data, FFTW_BACKWARD);
    F_VEC_2D az_fm_rate;
    D_VEC_1D az_freqs;

    if (tops_mode)
    {
        double num_replicas = std::abs(dc_rate * burst_length_seconds / prf);
        az_fm_rate = F_VEC_2D(num_packets, F_VEC_1D(num_samples));
        az_freqs = linspace(-num_replicas*prf/2, num_replicas*prf/2, num_packets);
    }
    else
    {
        az_freqs = linspace(-prf/2, prf/2, num_packets);
    }

    std::cout << "Azimuth Compressing" << std::endl;

    double a = 6378137.0;     // WGS84 Semi-Major
    double b = 6356752.3142;  // WGS84 Semi-Minor
    double e = 0.0067395;     // WGS84 Eccentricity

    #pragma omp parallel for
    for (int i = 0; i < num_packets; i++)
    {
        fftw_execute(forward_plans[i]);

        CF_VEC_1D& radar_data_row = radar_data[i];

        fftshift_in_place(radar_data_row);

        D_VEC_1D& pos = positions[i];

        double az_freq = az_freqs[i];
        double lat = pos[2] / pos[0];
        double earth_rad_num = std::pow(a * a * std::cos(lat), 2.0) + std::pow(b * b * std::sin(lat), 2.0);
        double earth_rad_denom = std::pow(a * std::cos(lat), 2.0) + std::pow(b * std::sin(lat), 2.0);
        double earth_radius = sqrt(earth_rad_num / earth_rad_denom);
        double sat_alt = std::sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
        double v_sat = velocities_norm[i];

        F_VEC_1D rcmc_factors(num_samples);
        F_VEC_1D effective_velocities(num_samples);

        for (int j = 0; j < num_samples; j++)
        {
            double slant_range = slant_ranges[j];
            double numerator = earth_radius*earth_radius + sat_alt*sat_alt - slant_range*slant_range;
            double denominator = 2 * earth_radius * sat_alt;
            double beta = numerator / denominator;
            double v_ground = earth_radius * v_sat * beta / sat_alt;
            double v_rel = sqrt(v_sat * v_ground);
            double rcmc = sqrt(
                1 - (std::pow(WAVELENGTH, 2.0) * std::pow(az_freq, 2.0)) / (4 * std::pow(v_rel, 2.0))
            );

            effective_velocities[j] = v_rel;
            rcmc_factors[j] = rcmc;

            double src_fm_rate  = 2.0 * std::pow(v_rel, 2.0) * std::pow(CENTER_FREQ, 3.0) * std::pow(rcmc, 2.0);
                   src_fm_rate /= SPEED_OF_LIGHT * slant_range * std::pow(az_freq, 2.0);
  
            std::complex<double> src_filter = 
                std::exp(-1.0 * I * PI * std::pow(range_freqs[j], 2.0) / src_fm_rate);

            radar_data_row[j] *= src_filter;
        }

        if (tops_mode) fftw_execute(inverse_plans[i]);

        for (int j = 0; j < num_samples; j++)
        {
            // TODO: Range migration correction using sinc-based interpolation

            if (tops_mode) az_freq += doppler_centroid[j];

            double v_rel = effective_velocities[j];
            double rcmc  = rcmc_factors[j];
            double slant_range = slant_ranges[j];

            std::complex<double> az_filter = 
                std::exp(4.0 * I * PI * slant_range * rcmc * CENTER_FREQ / SPEED_OF_LIGHT);

            if (tops_mode) 
            {
                az_fm_rate[i][j] = -(2 * std::pow(v_rel, 2.0) * std::pow(rcmc, 3.0)) / (WAVELENGTH * slant_range);
            }

            // TODO: Time correction and antenna pattern correction

            radar_data_row[j] *= (1 / double(num_samples)) * az_filter;
        }
    }

    destroy_fftw_plans(forward_plans);
    destroy_fftw_plans(inverse_plans);

    compute_axis_dft_in_place(radar_data, 0, 0, true);

    if (tops_mode)
    {
        int swath_number = packets[0].secondary_header("swath_number");
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
    }

    return radar_data;
}
