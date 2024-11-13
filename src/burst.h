#pragma once

#include "packet.h"
#include "misc_types.h"


/* Class representing a burst as a vector of packets, complex_samples, and chirps. */
class Burst {
private:
    int _burst_num;
    int _num_packets;

    std::string _swath;

    PACKET_VEC_1D _packets;

    CF_VEC_2D _signals;
    CF_VEC_2D _replica_chirps;

    void _set_pulses()
    {
        _signals.resize(_num_packets);
        _replica_chirps.resize(_num_packets);

        #pragma omp parallel for
        for (int i = 0; i < _num_packets; i++)
        {
            L0Packet packet    = _packets[i];
            _signals[i]        = packet.get_signal();
        }
    }

public:
    Burst(
        const std::string& filename,
        const std::string& swath,
        const int&         burst_num
    ) {
        _packets     = L0Packet::get_packets_in_bursts(filename, swath)[burst_num];
        _num_packets = _packets.size();
        _swath       = swath;
        _burst_num   = burst_num;

        _set_pulses();
    }

    Burst(
        std::ifstream&     data,
        const std::string& swath,
        const int&         burst_num
    ) {
        _packets     = L0Packet::get_packets_in_bursts(data, swath)[burst_num];
        _num_packets = _packets.size();
        _swath       = swath;
        _burst_num   = burst_num;

        _set_pulses();
    }

    Burst(
        const PACKET_VEC_1D& packets,
        const std::string&   swath,
        const int&           burst_num)
    {
        _packets     = packets;
        _num_packets = packets.size();
        _burst_num   = burst_num;
        _swath       = swath;

        _set_pulses();
    }

    int get_num_packets() { return _num_packets; }
    int get_burst_num()   { return _burst_num;   }

    std::string get_swath() { return _swath; }

    PACKET_VEC_1D get_packets()    { return _packets;         }
    CF_VEC_2D& get_signals()        { return _signals;         }
    CF_VEC_2D get_replica_chirps() { return _replica_chirps;  }

    L0Packet operator[](int index) { 
        if (index >= 0) return _packets[index];
        else            return _packets[_num_packets + index];
    }

    CF_VEC_1D get_signal(const int& index) 
    {
        if (index < 0 or index >= _num_packets)
        {
            throw std::out_of_range("Signal index is out-of-range.");
        }
        return _signals[index];
    }

    CF_VEC_1D get_replica_chirp(const int& index)
    {
        if (index < 0 or index >= _num_packets)
        {
            throw std::out_of_range("Signal index is out-of-range.");
        }
        return _replica_chirps[index];
    }
};


typedef struct std::vector<Burst>              BURST_VEC_1D;
typedef struct std::vector<std::vector<Burst>> BURST_VEC_2D;
