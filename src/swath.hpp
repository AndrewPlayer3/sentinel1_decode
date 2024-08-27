#pragma once

#include "burst.hpp"
#include "packet.hpp"


/* Class Representing a Swath of Bursts */
class Swath {
private:
    string        _swath;
    int           _num_bursts;
    vector<Burst> _bursts;

    void _set_bursts(ifstream& data)
    {
        vector<vector<L0Packet>> bursts = L0Packet::get_packets_in_bursts(data, _swath);

        _num_bursts = bursts.size();

        for (int i = 0; i < _num_bursts; i++)
        {
            _bursts.push_back(Burst(bursts[i], _swath, i));
        }
    }

public:
    Swath(const string& filename, const string& swath)
    {
        _swath = swath;
        ifstream data = open_file(filename);
        _set_bursts(data);
    }

    Swath(ifstream& data, const string& swath) 
    {
        _swath = swath;
        _set_bursts(data);
    }

    ~Swath() {}

    int           get_num_bursts() { return _num_bursts; }
    string        get_swath()      { return _swath;      }
    vector<Burst> get_bursts()     { return _bursts;     }


    Burst get_burst(const int& burst_num) 
    { 
        if (burst_num < 0 or burst_num >= _num_bursts)
        {
            throw out_of_range("Burst number is out-of-range.");
        }
        return _bursts[burst_num];
    }


    vector<vector<complex<float>>> get_all_signals()
    {
        vector<vector<complex<float>>> signals;

        for (Burst& burst : _bursts)
        {
            for (vector<complex<float>> signal : burst.get_signals())
            {
                signals.push_back(signal);
            }
        }
        return signals;
    }


    vector<vector<complex<float>>> get_all_replica_chirps()
    {
        vector<vector<complex<float>>> replica_chirps;

        for (Burst& burst : _bursts)
        {
            for (vector<complex<float>> replica_chirp : burst.get_replica_chirps())
            {
                replica_chirps.push_back(replica_chirp);
            }
        }
        return replica_chirps;
    }
};