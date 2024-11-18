#include "s1_decoder.h"


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

    int num_bursts = swath_packets.size();
    int num_packets_total = 0;
    INT_VEC_1D first_packet_index(num_bursts);
    int max_samples = 0;

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


CF_VEC_2D S1_Decoder::get_range_compressed_burst(const std::string& swath, const int& burst)
{
    PACKET_VEC_1D burst_packets = _echo_packets[swath][burst];

    return range_compress(burst_packets);
}


CF_VEC_2D S1_Decoder::get_range_compressed_swath_sm(const std::string& swath)
{
    PACKET_VEC_1D burst_packets = _echo_packets[swath][0];

    std::pair<PACKET_VEC_2D, int> azimuth_block_pair = get_azimuth_blocks(burst_packets);
    int max_size = azimuth_block_pair.second;
    PACKET_VEC_2D azimuth_blocks = std::move(azimuth_block_pair.first);

    std::cout << "Number of Azimuth Blocks: " << azimuth_blocks.size() << std::endl;

    if (azimuth_blocks.size() == 1)
    {
        return range_compress(azimuth_blocks[0]);
    }

    CF_VEC_2D range_compressed;

    for (int i = 0; i < azimuth_blocks.size(); i++)
    {
        CF_VEC_2D range_compressed_az_block = range_compress(azimuth_blocks[i]);

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


CF_VEC_2D S1_Decoder::get_range_compressed_swath_iw(const std::string& swath)
{
    PACKET_VEC_2D packets = _echo_packets[swath];

    int num_bursts = packets.size();

    CF_VEC_2D range_compressed;

    for (int i = 0; i < num_bursts; i++)
    {
        std::cout << "Range Compressing Burst #" << i << std::endl;
        CF_VEC_2D range_compressed_burst = range_compress(packets[i]);
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


CF_VEC_2D S1_Decoder::get_range_compressed_swath(const std::string& swath)
{
    if (is_sm(swath))
    {
        return get_range_compressed_swath_sm(swath);
    }
    else if (is_iw(swath))
    {
        return get_range_compressed_swath_iw(swath);
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
    return azimuth_compress(packets);
}


CF_VEC_2D S1_Decoder::get_azimuth_compressed_swath_iw(const std::string& swath)
{
    PACKET_VEC_2D packets = _echo_packets[swath];

    int num_bursts = packets.size();

    CF_VEC_2D azimuth_compressed;

    for (int i = 0; i < num_bursts; i++)
    {
        std::cout << "Azimuth Compressing Burst #" << i << std::endl;
        CF_VEC_2D azimuth_compressed_burst = azimuth_compress(packets[i]);
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


CF_VEC_2D S1_Decoder::get_azimuth_compressed_swath_sm(const std::string& swath)
{
    PACKET_VEC_1D burst_packets = _echo_packets[swath][0];

    std::pair<PACKET_VEC_2D, int> azimuth_block_pair = get_azimuth_blocks(burst_packets);
    int max_size = azimuth_block_pair.second;
    PACKET_VEC_2D azimuth_blocks = std::move(azimuth_block_pair.first);

    std::cout << "Number of Azimuth Blocks: " << azimuth_blocks.size() << std::endl;

    if (azimuth_blocks.size() == 1)
    {
        return azimuth_compress(azimuth_blocks[0]);
    }

    CF_VEC_2D azimuth_compressed;

    for (int i = 0; i < azimuth_blocks.size(); i++)
    {
        CF_VEC_2D azimuth_compressed_az_block = azimuth_compress(azimuth_blocks[i]);

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
        return get_azimuth_compressed_swath_sm(swath);
    }
    else if (is_iw(swath))
    {
        return get_azimuth_compressed_swath_iw(swath);
    }
    else
    {
        std::string msg = swath + " is not a supported swath type. Only IW and SM imaging modes are supported.";
        throw(std::invalid_argument(msg));
    }
}


CF_VEC_2D S1_Decoder::range_compress(PACKET_VEC_1D& packets)
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

    compute_axis_dft_in_place(range_compressed, 0, 1, true);

    std::cout << "FFT Shift" << std::endl;

    std::for_each(
        range_compressed.begin(), range_compressed.end(),
            [num_samples] (CF_VEC_1D& row) { 
                std::rotate(row.begin(), row.end()-(num_samples / 2), row.end());
            }
    );

    return range_compressed;
}


CF_VEC_2D S1_Decoder::azimuth_compress(PACKET_VEC_1D& packets)
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

    std::cout << "Converting Domains" << std::endl;

    range_compressed = transpose(range_compressed);

    compute_axis_dft_in_place(range_compressed, 0, 1, false);

    fftshift_in_place(range_compressed);

    range_compressed = transpose(range_compressed);

    std::cout << "Range Compressing" << std::endl;

    CF_VEC_1D replica = get_reference_function(packets[0].get_replica_chirp());

    std::vector<fftw_plan> plans = get_fftw_plans(range_compressed);

    #pragma omp parallel for
    for (int i = 0; i < num_packets; i++)
    {
        CF_VEC_1D& signal = range_compressed[i];

        fftw_execute(plans[i]);

        std::transform(
            replica.begin(), replica.end(),
                signal.begin(), signal.begin(),
                    [] (std::complex<double>& n, std::complex<double>& r) { return n * r;}
        );
    }

    destroy_fftw_plans(plans);

    std::cout << "Getting Auxillary Data" << std::endl;

    double burst_length_seconds = packets.back().get_time() - packets.front().get_time();

    D_VEC_1D times(num_packets);
    D_VEC_2D positions(num_packets, D_VEC_1D(3));
    D_VEC_1D velocities_norm(num_packets);

    for (int i = 0; i < num_packets; i++)
    {
        times[i] = packets[i].get_time();

        STATE_VECTOR state_vector = _state_vectors.interpolate(times[i]);

        D_VEC_1D v = state_vector.velocity;
        D_VEC_1D p = state_vector.position;

        positions[i]  = p;
        velocities_norm[i] = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    }

    double a = 6378137.0;     // WGS84 Semi-Major
    double b = 6356752.3142;  // WGS84 Semi-Minor
    double e = 0.0067395;     // WGS84 Eccentricity

    D_VEC_1D slant_ranges = packets[0].get_slant_ranges();

    double range_sample_rate = packets[0].get_range_sample_rate();
    double pulse_length = packets[0].get_pulse_length() * 1e-6;
    double pri = packets[0].get_pri() * 1e-6;
    double prf = 1 / pri;

    D_VEC_1D range_freqs = linspace(-range_sample_rate/2, range_sample_rate/2, num_samples);
    D_VEC_1D az_freqs = linspace(-prf/2, prf/2, num_packets);

    D_VEC_2D rcmc_factor_D(num_packets, D_VEC_1D(num_samples));

    std::cout << "Azimuth Compressing" << std::endl;

    #pragma omp parallel for
    for (int i = 0; i < num_packets; i++)
    {
        double latitude = positions[i][2] / positions[i][0];

        double earth_radius = sqrt(
            (pow(a*a*cos(latitude), 2) + pow(b*b*sin(latitude), 2))
            /
            (pow(a*cos(latitude), 2) + pow(b*sin(latitude), 2))
        );

        D_VEC_1D& p = positions[i];
        double pos = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
        double v_sat = velocities_norm[i];

        D_VEC_1D& rcmc_factor_row = rcmc_factor_D[i];
        CF_VEC_1D& range_compressed_row = range_compressed[i];

        for (int j = 0; j < num_samples; j++)
        {
            double numerator = earth_radius*earth_radius + pos*pos - slant_ranges[j]*slant_ranges[j];
            double denominator = 2 * earth_radius * pos;
            double beta = numerator / denominator;
            double v_ground = earth_radius * v_sat * beta / pos;

            double v_rel = sqrt(v_sat * v_ground);

            double rcmc_factor = sqrt(1 - 
                (WAVELENGTH*WAVELENGTH * az_freqs[i]*az_freqs[i]) / (4 * v_rel*v_rel)
            );
            double rcmc_shift  = slant_ranges[0] * ((1 / rcmc_factor) - 1);

            std::complex<double> rcmc_filter = exp(4.0 * I * PI * range_freqs[j] * rcmc_shift / SPEED_OF_LIGHT);

            rcmc_factor_row[j] = rcmc_factor;
            range_compressed_row[j] *= rcmc_filter;
        }
    }

    plans = get_fftw_plans(range_compressed);

    #pragma omp parallel for
    for (int i = 0; i < num_packets; i++)
    {
        fftw_execute(plans[i]);

        D_VEC_1D& rcmc_factor_row = rcmc_factor_D[i];
        CF_VEC_1D& range_compressed_row = range_compressed[i];

        fftshift_in_place(range_compressed_row);

        for (int j = 0; j < num_samples; j++)
        {
            range_compressed_row[j] *= (1 / double(num_samples)) * exp(4.0 * I * PI * slant_ranges[j] * rcmc_factor_row[j] / WAVELENGTH);
        }
    }

    rcmc_factor_D.clear();
    rcmc_factor_D.shrink_to_fit();

    destroy_fftw_plans(plans);

    compute_axis_dft_in_place(range_compressed, 0, 0, true);

    return range_compressed;
}
