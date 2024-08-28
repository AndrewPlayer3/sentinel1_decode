#include "image_formation.h"


CF_VEC_1D get_reference_function(const CF_VEC_1D& replica_chirp)
{
    int num_samples = replica_chirp.size();

    CF_VEC_1D replica_chirp_conj = conjugate(replica_chirp);
    CF_VEC_1D weighted_chirp     = apply_hanning_window(replica_chirp_conj);
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
        reference[i] = weighted_chirp[i] / energy;
    }
    return reference;
}


CF_VEC_1D pulse_compression(
    const CF_VEC_1D& signal,
    const CF_VEC_1D& replica_chirp
) {
    int num_samples     = signal.size();
    int replica_samples = replica_chirp.size();

    CF_VEC_1D pulse_compressed(num_samples);

    CF_VEC_1D signal_fft = compute_1d_dft(signal,  0, false);
    CF_VEC_1D reference  = get_reference_function(replica_chirp);

    CF_VEC_1D chirp(num_samples);

    for (int i = 0; i < num_samples; i++)
    {
        if (i < replica_samples) chirp[i] = reference[i];
        else                     chirp[i] = 0.0;
    }

    chirp = compute_1d_dft(chirp,  0, false);

    for (int i = 0; i < num_samples; i++)
    {
        pulse_compressed[i] = signal_fft[i] * chirp[i];
    }    

    return compute_1d_dft(pulse_compressed, num_samples, true);
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