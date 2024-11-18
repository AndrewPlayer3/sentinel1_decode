cmake_minimum_required(VERSION 3.10)

project(SignalProcessingProject)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O3")

find_package( OpenCV REQUIRED )

# Detect platform
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CLANG_COMPILER "/opt/homebrew/opt/llvm/bin/clang++")
else()
    set(CLANG_COMPILER "/bin/g++")
endif()

# Set include directories and libraries for macOS
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(HOMEBREW_INCLUDE "/opt/homebrew/include")
    set(HOMEBREW_LIB "/opt/homebrew/lib")
    include_directories(${HOMEBREW_INCLUDE} ${OpenCV_INCLUDE_DIRS})
    link_directories(${HOMEBREW_LIB})
else()
    include_directories(${OpenCV_INCLUDE_DIRS})
endif()

# Executable: main
add_executable(s1_print src/print_main.cpp src/cli.cpp src/packet.cpp src/signal_processing.cpp src/image_formation.cpp src/decoding_utils.cpp src/aux_decoding.cpp src/state_vectors.cpp src/s1_decoder.cpp)
set_target_properties(s1_print PROPERTIES OUTPUT_NAME "../bin/s1_print")
set_target_properties(s1_print PROPERTIES CXX_COMPILER "${CLANG_COMPILER}")
target_link_libraries(s1_print PRIVATE fftw3 m)

# Executable: write
add_executable(s1_write_tiff src/write_main.cpp src/cli.cpp src/image_write.cpp src/image_formation.cpp src/packet.cpp src/decoding_utils.cpp src/aux_decoding.cpp src/signal_processing.cpp src/state_vectors.cpp src/s1_decoder.cpp)
set_target_properties(s1_write_tiff PROPERTIES OUTPUT_NAME "../bin/s1_write")
target_link_libraries(s1_write_tiff PRIVATE fftw3 m tiff )

# Executable: plot
add_executable(s1_plot src/plot_main.cpp src/cli.cpp src/plot.cpp src/image_formation.cpp src/aux_decoding.cpp src/packet.cpp src/decoding_utils.cpp src/signal_processing.cpp src/state_vectors.cpp src/s1_decoder.cpp)
set_target_properties(s1_plot PROPERTIES OUTPUT_NAME "../bin/s1_plot")
target_link_libraries(s1_plot PUBLIC ${OpenCV_LIBS})
target_link_libraries(s1_plot PRIVATE fftw3 m)
