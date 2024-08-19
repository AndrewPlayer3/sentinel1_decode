#!/bin/bash

mkdir bin 2>/dev/null

echo "Compiling main to bin/main."

if [[ "$OSTYPE" == "linux-gnu"* ]]; then
	g++ -fopenmp -std=c++20 -O3 src/main.cpp src/packet.cpp src/packet_decoding.cpp src/decoding_utils.cpp -o bin/main
elif [[ "$OSTYPE" == "darwin"* ]]; then
	/opt/homebrew/opt/llvm/bin/clang++ -fopenmp -std=c++20 -O3 src/main.cpp src/packet.cpp src/perf.cpp src/decoding_utils.cpp -o bin/main
else
	echo "$OSTYPE is not supported by this script, see the commands for what needs to be compiled."
fi


# See https://github.com/lava/matplotlib-cpp/blob/master/README.md for more information on matplotlib-cpp.
# Needs to be linked to Python.h and <numpy>, which in my case are both in my conda environment.

echo "Compiling plotting to bin/plot"

if [[ "$OSTYPE" == "linux-gnu"* ]]; then
	g++ -fopenmp -std=c++20 -O3 -w src/plotting.cpp src/packet.cpp src/decoding_utils.cpp -I/$HOME/miniforge3/envs/signal_processing/include/python3.12 -I$HOME/miniforge3/envs/signal_processing/lib/python3.12/site-packages/numpy/_core/include -lpython3.12 -o bin/plot
elif [[ "$OSTYPE" == "darwin"* ]]; then
	/opt/homebrew/opt/llvm/bin/clang++ -fopenmp -std=c++20 -O3 -w src/plotting.cpp src/packet.cpp src/decoding_utils.cpp -I/opt/homebrew/include -I/$HOME/miniforge3/envs/signal_processing/include/python3.12 -I$HOME/miniforge3/envs/signal_processing/lib/python3.12/site-packages/numpy/_core/include -L/opt/homebrew/lib -L/opt/homebrew/Frameworks/Python.framework/Versions/3.12/lib -lfftw3 -lm -lpython3.12 -o bin/plot
else
	echo "$OSTYPE is not supported by this script, see the commands for what needs to be compiled."
fi