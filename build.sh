#!/bin/bash

mkdir bin 2>/dev/null

export MAIN_OPTIONS="-fopenmp -std=c++20 -O3"
export MAIN_SRC="src/main.cpp src/cli.cpp src/packet.cpp src/perf.cpp src/decoding_utils.cpp src/aux_decoding.cpp"
export MAIN_PATH="bin/main"

export WRITE_OPTIONS="-fopenmp -std=c++20 -O3"
export WRITE_SRC="src/write_main.cpp src/cli.cpp src/image_write.cpp src/image_formation.cpp src/packet.cpp src/decoding_utils.cpp src/aux_decoding.cpp src/signal_processing.cpp"
export WRITE_PATH="bin/write"
export WRITE_LIBS="-lfftw3f_omp -lfftw3f -lm -ltiff"

export PLOT_OPTIONS="-fopenmp -std=c++20 -O3"
export PLOT_SRC="src/plot_main.cpp src/cli.cpp src/plot.cpp src/image_formation.cpp src/packet.cpp src/decoding_utils.cpp src/signal_processing.cpp"
export PLOT_PATH="bin/plot"
export PLOT_LIBS="-lfftw3f_omp -lfftw3f -lm -lpython3.12"

export PYTHON_INCLUDE="$HOME/miniforge3/envs/signal_processing/include/python3.12"
export NUMPY_INCLUDE="$HOME/miniforge3/envs/signal_processing/lib/python3.12/site-packages/numpy/_core/include"

export HOMEBREW_INCLUDE="/opt/homebrew/include"
export HOMEBREW_LIB="/opt/homebrew/lib"
export HOMEBREW_PYTHON_INCLUDE="/opt/homebrew/Frameworks/Python.framework/Versions/3.12/lib"


echo "Compiling main to bin/main."
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
	g++ $MAIN_OPTIONS $MAIN_SRC -o $MAIN_PATH
elif [[ "$OSTYPE" == "darwin"* ]]; then
	/opt/homebrew/opt/llvm/bin/clang++ $MAIN_OPTIONS $MAIN_SRC -o $MAIN_PATH
else
	echo "$OSTYPE is not supported by this script, see the commands for what needs to be compiled."
fi


echo "Compiling image writing to bin/write"
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
	g++ -L/usr/lib/x86_64-linux-gnu/ $WRITE_OPTIONS $WRITE_SRC  $WRITE_LIBS -o $WRITE_PATH
elif [[ "$OSTYPE" == "darwin"* ]]; then
	/opt/homebrew/opt/llvm/bin/clang++ $WRITE_OPTIONS -w $WRITE_SRC -I$HOMEBREW_INCLUDE -I$PYTHON_INCLUDE \
	                                  -I$NUMPY_INCLUDE -L$HOMEBREW_LIB -L$HOMEBREW_PYTHON_INCLUDE $WRITE_LIBS \
									  -lpython3.12 -o $WRITE_PATH
else
	echo "$OSTYPE is not supported by this script, see the commands for what needs to be compiled."
fi


echo "Compiling plotting to bin/plot."
# # See https://github.com/lava/matplotlib-cpp/blob/master/README.md for more information on matplotlib-cpp.
# # Needs to be linked to Python.h and <numpy>, which in my case are both in my conda environment.
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
	g++ $PLOT_OPTIONS -w $PLOT_SRC -I$PYTHON_INCLUDE -I$NUMPY_INCLUDE $PLOT_LIBS -o $PLOT_PATH
elif [[ "$OSTYPE" == "darwin"* ]]; then
	/opt/homebrew/opt/llvm/bin/clang++ $PLOT_OPTIONS -w $PLOT_SRC -I$HOMEBREW_INCLUDE -I$PYTHON_INCLUDE \
	                                   -I$NUMPY_INCLUDE -L$HOMEBREW_LIB -L$HOMEBREW_PYTHON_INCLUDE $PLOT_LIBS \
									   -o $PLOT_PATH
else
	echo "$OSTYPE is not supported by this script, see the commands for what needs to be compiled."
fi
