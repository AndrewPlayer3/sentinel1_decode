#!/bin/bash

mkdir bin 2>/dev/null

echo "Compiling main to bin/main."

if [[ "$OSTYPE" == "linux-gnu"* ]]; then
	g++ -fopenmp -std=c++20 -O3 src/main.cpp src/packet.cpp src/perf.cpp src/decoding_utils.cpp src/aux_decoding.cpp -o bin/main
elif [[ "$OSTYPE" == "darwin"* ]]; then
	/opt/homebrew/opt/llvm/bin/clang++ -fopenmp -std=c++20 -O3 src/main.cpp src/packet.cpp src/perf.cpp src/decoding_utils.cpp src/aux_decoding.cpp -o bin/main
else
	echo "$OSTYPE is not supported by this script, see the commands for what needs to be compiled."
fi


# See https://github.com/lava/matplotlib-cpp/blob/master/README.md for more information on matplotlib-cpp.
# Needs to be linked to Python.h and <numpy>, which in my case are both in my conda environment.

echo "Compiling plotting to bin/plot."

export PYTHON_INCLUDE="$HOME/miniforge3/envs/signal_processing/include/python3.12"
export NUMPY_INCLUDE="$HOME/miniforge3/envs/signal_processing/lib/python3.12/site-packages/numpy/_core/include"

if [[ "$OSTYPE" == "linux-gnu"* ]]; then

	g++ -fopenmp -std=c++20 -O3 -w src/plotting.cpp src/packet.cpp src/decoding_utils.cpp -I$PYTHON_INCLUDE -I$NUMPY_INCLUDE -lfftw3 -lm -lpython3.12 -o bin/plot

elif [[ "$OSTYPE" == "darwin"* ]]; then

	export PYTHON_INCLUDE="$HOME/miniforge3/envs/signal_processing/include/python3.12"
	export NUMPY_INCLUDE="$HOME/miniforge3/envs/signal_processing/lib/python3.12/site-packages/numpy/_core/include"
	export HOMEBREW_INCLUDE="/opt/homebrew/include"
	export HOMEBREW_LIB="/opt/homebrew/lib"
	export HOMEBREW_PYTHON_INCLUDE="/opt/homebrew/Frameworks/Python.framework/Versions/3.12/lib"

	/opt/homebrew/opt/llvm/bin/clang++ -fopenmp -std=c++20 -O3 -w src/plotting.cpp src/packet.cpp src/decoding_utils.cpp -I$HOMEBREW_INCLUDE -I$PYTHON_INCLUDE -I$NUMPY_INCLUDE -L$HOMEBREW_LIB -L$HOMEBREW_PYTHON_INCLUDE -lfftw3 -lm -lpython3.12 -o bin/plot

else
	echo "$OSTYPE is not supported by this script, see the commands for what needs to be compiled."
fi