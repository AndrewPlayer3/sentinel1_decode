cmake_minimum_required(VERSION 3.10)

project(s1_decode)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O3")

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CLANG_COMPILER "/opt/homebrew/opt/llvm/bin/clang++")
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(CLANG_COMPILER "/bin/g++")
endif()
# Finding all only seems to matter on Windows
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    find_package( OpenCV CONFIG REQUIRED )
    find_package( FFTW3 CONFIG REQUIRED )
    find_package( TIFF CONFIG REQUIRED )
else()
    find_package( OpenCV REQUIRED )
endif()

# CMake attempts to find the dependencies automatically. If it can't, setting the environment
# variables manually prior to configuring CMake should allow things to work. 
#
# Here are some examples:
#
# export OpenCV_LIBS /opt/homebrew/include
# export OpenCV_INCLUDE_DIRS /opt/homebrew/lib
#
# export FFTW3_LIBS /usr/lib/x86_64-linux-gnu/
# export FFTW3_INCLUDE_DIRS /usr/include
#
# export TIFF_LIBS /usr/lib/x86_64-linux-gnu/
# export TIFF_INCLUDE_DIRS /usr/include

find_package( OpenCV )
find_package( FFTW3 )
find_package( TIFF )

if(NOT OpenCV_FOUND)
    set(OpenCV_LIBS "$ENV{OpenCV_LIBS}")
    set(OpenCV_INCLUDE_DIRS "$ENV{OpenCV_INCLUDE_DIRS}")
endif()

if(NOT FFTW3_FOUND)
    set(FFTW3_LIBS "$ENV{FFTW3_LIBS}")
    set(FFTW3_INCLUDE_DIRS "$ENV{FFTW3_INCLUDE_DIRS}")
endif()

if(NOT TIFF_FOUND)
    set(TIFF_LIBS "$ENV{TIFF_LIBS}")
    set(TIFF_INCLUDE_DIRS "$ENV{TIFF_INCLUDE_DIRS}")
endif()

include_directories(${OpenCV_INCLUDE_DIRS} ${FFTW3_INCLUDE_DIRS} ${TIFF_INCLUDE_DIRS})
link_libraries(${OpenCV_LIBS} ${FFTW3_LIBS} ${TIFF_LIBS})

# Executable: main
add_executable(
    s1_print
    src/print_main.cpp
    src/cli.cpp
    src/packet.cpp
    src/signal_processing.cpp
    src/image_formation.cpp
    src/decoding_utils.cpp
    src/aux_decoding.cpp
    src/state_vectors.cpp
    src/s1_decoder.cpp
)
set_target_properties(s1_print PROPERTIES OUTPUT_NAME "../bin/s1_print")
target_link_libraries(s1_print PRIVATE fftw3 m)

# Executable: write
add_executable(
    s1_write_tiff
    src/write_main.cpp
    src/cli.cpp
    src/image_write.cpp
    src/image_formation.cpp
    src/packet.cpp
    src/decoding_utils.cpp
    src/aux_decoding.cpp
    src/signal_processing.cpp
    src/state_vectors.cpp
    src/s1_decoder.cpp
)
set_target_properties(s1_write_tiff PROPERTIES OUTPUT_NAME "../bin/s1_write")
target_link_libraries(s1_write_tiff PRIVATE fftw3 tiff )

# Executable: plot
add_executable(
    s1_plot
    src/plot_main.cpp
    src/cli.cpp
    src/plot.cpp
    src/image_formation.cpp
    src/aux_decoding.cpp
    src/packet.cpp
    src/decoding_utils.cpp
    src/signal_processing.cpp
    src/state_vectors.cpp
    src/s1_decoder.cpp
)
set_target_properties(s1_plot PROPERTIES OUTPUT_NAME "../bin/s1_plot")
target_link_libraries(s1_plot PUBLIC ${OpenCV_LIBS})
target_link_libraries(s1_plot PRIVATE fftw3)
