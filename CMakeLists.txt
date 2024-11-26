cmake_minimum_required(VERSION 3.10)

project(s1_decode)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O3")

# Finding all only seems to matter on Windows
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    find_package( OpenCV CONFIG REQUIRED )
    find_package( FFTW3 CONFIG REQUIRED )
    find_package( TIFF CONFIG REQUIRED )
else()
    find_package( OpenCV REQUIRED )
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
target_link_libraries(s1_print PRIVATE fftw3)

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
