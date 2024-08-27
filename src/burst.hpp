#pragma once

#include "packet.hpp"

/* Class representing a burst as a vector of packets, complex_samples, and chirps. */
class Burst {
private:
    int _burst_num;
    int _num_packets;

    string _swath;

    vector<L0Packet> _packets;
    vector<vector<complex<float>>> _signals;
    vector<vector<complex<float>>> _replica_chirps;

    void _set_pulses()
    {
        _signals.resize(_num_packets);
        _replica_chirps.resize(_num_packets);

        #pragma omp parallel for
        for (int i = 0; i < _num_packets; i++)
        {
            L0Packet packet    = _packets[i];
            _signals[i]        = packet.get_signal();
            _replica_chirps[i] = packet.get_replica_chirp();
        }
    }

public:
    Burst(
        const string& filename,
        const string& swath,
        const int&    burst_num
    ) {
        _packets     = L0Packet::get_packets_in_bursts(filename, swath)[burst_num];
        _num_packets = _packets.size();
        _swath       = swath;
        _burst_num   = burst_num;

        _set_pulses();
    }

    Burst(
        ifstream&     data,
        const string& swath,
        const int&    burst_num
    ) {
        _packets     = L0Packet::get_packets_in_bursts(data, swath)[burst_num];
        _num_packets = _packets.size();
        _swath       = swath;
        _burst_num   = burst_num;

        _set_pulses();
    }

    Burst(
        const vector<L0Packet>& packets,
        const string& swath,
        const int&    burst_num)
    {
        _packets     = packets;
        _num_packets = packets.size();
        _burst_num   = burst_num;
        _swath       = swath;

        _set_pulses();
    }

    int get_num_packets() { return _num_packets; }
    int get_burst_num()   { return _burst_num;   }

    string get_swath() { return _swath; }

    vector<L0Packet>               get_packets()        { return _packets;         }
    vector<vector<complex<float>>> get_signals()        { return _signals;         }
    vector<vector<complex<float>>> get_replica_chirps() { return _replica_chirps;  }

    vector<complex<float>> get_signal(const int& index) 
    {
        if (index < 0 or index >= _num_packets)
        {
            throw out_of_range("Signal index is out-of-range.");
        }
        return _signals[index];
    }

    vector<complex<float>> get_replica_chirp(const int& index)
    {
        if (index < 0 or index >= _num_packets)
        {
            throw out_of_range("Signal index is out-of-range.");
        }
        return _replica_chirps[index];
    }
};