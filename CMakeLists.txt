cmake_minimum_required(VERSION 3.10)

project(SignalProcessingProject)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O3")

find_package( OpenCV REQUIRED )

# Detect platform
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CLANG_COMPILER "/opt/homebrew/opt/llvm/bin/clang++")
else()
    set(CLANG_COMPILER "${CMAKE_CXX_COMPILER}")
endif()

# Set include directories and libraries for macOS
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(HOMEBREW_INCLUDE "/opt/homebrew/include")
    set(HOMEBREW_LIB "/opt/homebrew/lib")
    set(HOMEBREW_PYTHON_INCLUDE "/opt/homebrew/Frameworks/Python.framework/Versions/3.12/lib")
    set(PYTHON_LIBS "$ENV{HOME}/miniforge3/envs/signal_processing/include/python3.12/lib")
    set(PYTHON_INCLUDE "$ENV{HOME}/miniforge3/envs/signal_processing/include/python3.12")
    set(NUMPY_INCLUDE "$ENV{HOME}/miniforge3/envs/signal_processing/lib/python3.12/site-packages/numpy/_core/include")
else()
    set(PYTHON_INCLUDE "$ENV{HOME}/miniforge3/envs/signal_processing/include/python3.12")
    set(NUMPY_INCLUDE "$ENV{HOME}/miniforge3/envs/signal_processing/lib/python3.12/site-packages/numpy/_core/include")
endif()

include_directories(${PYTHON_INCLUDE} ${NUMPY_INCLUDE} ${HOMEBREW_INCLUDE} ${OpenCV_INCLUDE_DIRS})

link_directories(${PYTHON_LIBS} ${HOMEBREW_LIB} ${HOMEBREW_PYTHON_INCLUDE})

add_subdirectory("matplotplusplus")

# Executable: main
add_executable(s1_print src/main.cpp src/cli.cpp src/packet.cpp src/signal_processing.cpp src/image_formation.cpp src/perf.cpp src/decoding_utils.cpp src/aux_decoding.cpp src/state_vectors.cpp)
set_target_properties(s1_print PROPERTIES OUTPUT_NAME "../bin/s1_print.exe")
set_target_properties(s1_print PROPERTIES CXX_COMPILER "${CLANG_COMPILER}")
target_link_libraries(s1_print PRIVATE fftw3_omp fftw3 m)

# Executable: write
add_executable(s1_write_tiff src/write_main.cpp src/cli.cpp src/image_write.cpp src/image_formation.cpp src/packet.cpp src/decoding_utils.cpp src/aux_decoding.cpp src/signal_processing.cpp src/state_vectors.cpp)
set_target_properties(s1_write_tiff PROPERTIES OUTPUT_NAME "../bin/s1_write_tiff.exe")
target_link_libraries(s1_write_tiff PRIVATE fftw3_omp fftw3 m tiff )
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    link_directories(${HOMEBREW_LIB})
endif()

# Executable: plot
add_executable(s1_plot src/plot_main.cpp src/cli.cpp src/plot.cpp src/image_formation.cpp src/aux_decoding.cpp src/packet.cpp src/decoding_utils.cpp src/signal_processing.cpp src/state_vectors.cpp)
set_target_properties(s1_plot PROPERTIES OUTPUT_NAME "../bin/s1_plot.exe")
target_link_libraries(s1_plot PUBLIC matplot ${OpenCV_LIBS})
target_link_libraries(s1_plot PRIVATE python3.12)
target_link_libraries(s1_plot PRIVATE fftw3_omp fftw3 m)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    link_directories(${HOMEBREW_LIB} ${HOMEBREW_PYTHON_INCLUDE})
endif()
