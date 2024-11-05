#pragma once

#include "burst.h"
#include "packet.h"
#include "misc_types.h"


/* Class Representing a Swath of Bursts */
class Swath {
private:
    std::string  _swath;
    BURST_VEC_1D _bursts;
    int          _num_bursts;

    void _set_bursts(std::ifstream& data)
    {
        PACKET_VEC_2D bursts = L0Packet::get_packets_in_bursts(data, _swath);
    
        _num_bursts = bursts.size();

        for (int i = 0; i < _num_bursts; i++)
        {
            _bursts.push_back(Burst(bursts[i], _swath, i));
        }
    }

public:
    Swath(const std::string& filename, const std::string& swath)
    {
        _swath = swath;
        std::ifstream data = open_file(filename);
        _set_bursts(data);
    }

    Swath(std::ifstream& data, const std::string& swath) 
    {
        _swath = swath;
        _set_bursts(data);
    }

    ~Swath() {}

    int          get_num_bursts() { return _num_bursts; }
    std::string  get_swath()      { return _swath;      }
    BURST_VEC_1D get_bursts()     { return _bursts;     }

    int get_num_packets() 
    {
        int num_packets = 0;
        for (int i = 0; i < _num_bursts; i++)
        {
            num_packets += _bursts[i].get_num_packets();
        }
        return num_packets;
    }


    Burst get_burst(const int& burst_num) 
    { 
        if (burst_num < 0 or burst_num >= _num_bursts)
        {
            throw std::out_of_range("Burst number is out-of-range.");
        }
        return _bursts[burst_num];
    }


    CF_VEC_2D get_all_signals()
    {
        CF_VEC_2D signals;

        for (Burst& burst : _bursts)
        {
            for (CF_VEC_1D& signal : burst.get_signals())
            {
                signals.push_back(signal);
            }
        }
        return signals;
    }


    CF_VEC_2D get_all_replica_chirps()
    {
        CF_VEC_2D replica_chirps;

        for (Burst& burst : _bursts)
        {
            for (CF_VEC_1D& replica_chirp : burst.get_replica_chirps())
            {
                replica_chirps.push_back(replica_chirp);
            }
        }
        return replica_chirps;
    }
};


typedef struct std::vector<Swath>              SWATH_VEC_1D;
typedef struct std::vector<std::vector<Swath>> SWATH_VEC_2D;
