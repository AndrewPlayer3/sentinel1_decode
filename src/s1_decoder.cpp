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


CF_VEC_2D S1_Decoder::get_burst(const std::string& swath, const int& burst)
{
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


std::pair<PACKET_VEC_2D, int> S1_Decoder::get_azimuth_blocks(PACKET_VEC_1D& packets)
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


CF_VEC_2D S1_Decoder::get_range_compressed_burst(const std::string& swath, const int& burst, bool range_doppler)
{
    PACKET_VEC_1D burst_packets = _echo_packets[swath][burst];

    return _range_compress(burst_packets, true, range_doppler);
}


CF_VEC_2D S1_Decoder::_get_range_compressed_swath_sm(const std::string& swath, bool range_doppler)
{
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


CF_VEC_2D S1_Decoder::_get_range_compressed_swath_iw(const std::string& swath, bool range_doppler)
{
    PACKET_VEC_2D packets = _echo_packets[swath];

    int num_bursts = packets.size();

    CF_VEC_2D range_compressed;

    for (int i = 0; i < num_bursts; i++)
    {
        std::cout << "Range Compressing Burst #" << i << std::endl;
        CF_VEC_2D range_compressed_burst = _range_compress(packets[i], true, range_doppler);
        int num_azimuth_lines = range_compressed_burst.size();
        for (int i = 0; i < num_azimuth_lines; i++)
        {
            CF_VEC_1D row = range_compressed_burst.front();
            range_compressed_burst.erase(range_compressed_burst.begin());
            range_compressed.push_back(row);
        }
    }

    return range_compressed;
}


CF_VEC_2D S1_Decoder::get_range_compressed_swath(const std::string& swath, bool range_doppler)
{
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
        std::string msg = swath + " is not a supported swath type. Only IW and SM imaging modes are supported.";
        throw(std::invalid_argument(msg));
    }
}


CF_VEC_2D S1_Decoder::get_azimuth_compressed_burst(const std::string& swath, const int& burst)
{
    PACKET_VEC_1D packets = _echo_packets[swath][burst];
    return _azimuth_compress(packets);
}


CF_VEC_2D S1_Decoder::_get_azimuth_compressed_swath_iw(const std::string& swath)
{
    PACKET_VEC_2D packets = _echo_packets[swath];

    int num_bursts = packets.size();

    CF_VEC_2D azimuth_compressed;

    for (int i = 0; i < num_bursts; i++)
    {
        std::cout << "Azimuth Compressing Burst #" << i << std::endl;
        CF_VEC_2D azimuth_compressed_burst = _azimuth_compress(packets[i]);
        int num_lines = azimuth_compressed_burst.size();
        for (int i = 0; i < num_lines; i++)
        {
            CF_VEC_1D row = azimuth_compressed_burst.front();
            azimuth_compressed_burst.erase(azimuth_compressed_burst.begin());
            azimuth_compressed.push_back(row);
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

        for (int i = 0; i < rows; i++)
        {
            CF_VEC_1D row = azimuth_compressed_az_block.front();
            azimuth_compressed_az_block.erase(azimuth_compressed_az_block.begin());
            for (int j = 0; j < diff; j++)
            {
                row.emplace(row.begin(), std::complex<double>(0.0));
            }
            azimuth_compressed.push_back(row);
        }
    }

    return azimuth_compressed;
}


CF_VEC_2D S1_Decoder::get_azimuth_compressed_swath(const std::string& swath)
{
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
        std::string msg = swath + " is not a supported swath type. Only IW and SM imaging modes are supported.";
        throw(std::invalid_argument(msg));
    }
}


CF_VEC_2D S1_Decoder::_range_compress(PACKET_VEC_1D& packets, bool do_ifft, bool do_azimuth_fft)
{
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


CF_VEC_2D S1_Decoder::_azimuth_compress(PACKET_VEC_1D& packets)
{
    L0Packet first_packet = packets[0];
    int num_packets = packets.size();
    int num_samples = 2 * packets[0].get_num_quads();

    CF_VEC_2D range_compressed = _range_compress(packets, true, false);

    D_VEC_1D slant_ranges = packets[0].get_slant_ranges();

    F_VEC_1D v_0 = _state_vectors.velocities[0];

    double v_norm = std::sqrt(std::pow(v_0[0], 2.0) + std::pow(v_0[1], 2.0) + std::pow(v_0[2], 2.0));
    double range_sample_rate = packets[0].get_range_sample_rate();
    double pulse_length = packets[0].get_pulse_length() * 1e-6;
    double pri = packets[0].get_pri() * 1e-6;
    double prf = 1 / pri;
    double burst_length_seconds = double(num_packets) / prf;
    double dc_rate = get_doppler_centroid_rate(packets, v_norm);

    F_VEC_1D doppler_centroid = get_doppler_centroid(
        range_compressed,
        dc_rate,
        burst_length_seconds,
        first_packet
    );

    range_compressed = azimuth_frequency_ufr(
        range_compressed,
        doppler_centroid,
        first_packet,
        dc_rate,
        burst_length_seconds,
        prf
    );
    num_packets = range_compressed.size();

    double oversample_factor = double(num_packets) / double(packets.size());
    double t0 = packets[0].get_time();
    double t1 = packets[1].get_time();
    double time_delta = (t1 - t0) / oversample_factor;

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

    double num_replicas = std::abs(dc_rate * burst_length_seconds / prf);
    D_VEC_1D range_freqs = linspace(-range_sample_rate/2, range_sample_rate/2, num_samples);
    D_VEC_1D az_freqs = linspace(-num_replicas*prf/2, num_replicas*prf/2, num_packets);

    F_VEC_2D ka(num_packets, F_VEC_1D(num_samples));

    std::cout << "Azimuth Compressing" << std::endl;

    double a = 6378137.0;     // WGS84 Semi-Major
    double b = 6356752.3142;  // WGS84 Semi-Minor
    double e = 0.0067395;     // WGS84 Eccentricity

    #pragma omp parallel for
    for (int i = 0; i < num_packets; i++)
    {
        D_VEC_1D& position = positions[i];
        F_VEC_1D& ka_row = ka[i];
        CF_VEC_1D& range_compressed_row = range_compressed[i];

        double latitude = position[2] / position[0];

        double earth_radius = sqrt(
            (pow(a*a*cos(latitude), 2) + pow(b*b*sin(latitude), 2)) 
            /
            (pow(a*cos(latitude), 2) + pow(b*sin(latitude), 2))
        );

        double pos = sqrt(position[0]*position[0] + position[1]*position[1] + position[2]*position[2]);
        double v_sat = velocities_norm[i];

        for (int j = 0; j < num_samples; j++)
        {
            double numerator = earth_radius*earth_radius + pos*pos - slant_ranges[j]*slant_ranges[j];
            double denominator = 2 * earth_radius * pos;
            double beta = numerator / denominator;
            double v_ground = earth_radius * v_sat * beta / pos;
            double v_rel = sqrt(v_sat * v_ground);
            double rcmc_factor = sqrt(1 - (WAVELENGTH*WAVELENGTH * az_freqs[i]*az_freqs[i]) / (4 * v_rel*v_rel));

            range_compressed_row[j] *= (1 / double(num_samples)) * exp(4.0 * I * PI * slant_ranges[j] * rcmc_factor * CENTER_FREQ / SPEED_OF_LIGHT);
            ka_row[j] = -(2 * std::pow(v_rel, 2.0) * std::pow(rcmc_factor, 3.0)) / (WAVELENGTH * slant_ranges[j]); 
        }
    }

    compute_axis_dft_in_place(range_compressed, 0, 0, true);

    range_compressed = azimuth_time_ufr(range_compressed, doppler_centroid, ka, first_packet, dc_rate, burst_length_seconds, prf);

    return range_compressed;
}
