cmake_minimum_required(VERSION 3.10)
project(s1_decode CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Output directory
set(BIN_DIR ${CMAKE_BINARY_DIR}/bin)
file(MAKE_DIRECTORY ${BIN_DIR})

# Compiler flags
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang|GNU")
    add_compile_options(-O3)
endif()

# Add Homebrew and common prefixes
if(APPLE)
    list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew" "/usr/local")
    set(CMAKE_INSTALL_RPATH "/opt/homebrew/lib")
endif()

# OpenMP (portable detection)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "Found OpenMP")
endif()

# Dependencies via find_package
find_package(OpenCV QUIET)
# find_package(FFTW3 QUIET)
find_package(TIFF QUIET)

# Fallback: pkg-config
find_package(PkgConfig REQUIRED)

if(NOT OpenCV_FOUND)
    pkg_check_modules(OpenCV IMPORTED_TARGET opencv4)
endif()

pkg_check_modules(FFTW3 REQUIRED IMPORTED_TARGET fftw3)

if(NOT TIFF_FOUND)
    pkg_check_modules(TIFF IMPORTED_TARGET libtiff-4)
endif()

# Debug info
message(STATUS "FFTW3 include: ${FFTW3_INCLUDE_DIRS}")
message(STATUS "FFTW3 libs: ${FFTW3_LIBRARIES}")
message(STATUS "TIFF include: ${TIFF_INCLUDE_DIRS}")
message(STATUS "TIFF libs: ${TIFF_LIBRARIES}")
message(STATUS "OpenCV include: ${OpenCV_INCLUDE_DIRS}")
message(STATUS "OpenCV libs: ${OpenCV_LIBS}")

# Shared sources
set(SRC_COMMON
    src/cli.cpp
    src/packet.cpp
    src/decoding_utils.cpp
    src/signal_processing.cpp
    src/aux_decoding.cpp
    src/state_vectors.cpp
    src/image_formation.cpp
    src/s1_decoder.cpp
    src/doppler.cpp
)

# Function to define executables cleanly
function(add_s1_executable name main_file extra_srcs)
    add_executable(${name})
    target_sources(${name} PRIVATE ${main_file} ${SRC_COMMON} ${extra_srcs})
    set_target_properties(${name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})

    # Link OpenCV
    if(OpenCV_FOUND)
        target_include_directories(${name} PRIVATE ${OpenCV_INCLUDE_DIRS})
        target_link_libraries(${name} PRIVATE ${OpenCV_LIBS})
    elseif(TARGET PkgConfig::OpenCV)
        target_link_libraries(${name} PRIVATE PkgConfig::OpenCV)
    endif()

    # Link FFTW
    target_link_libraries(${name} PRIVATE PkgConfig::FFTW3)

    # Link TIFF
    if(TIFF_FOUND)
        target_include_directories(${name} PRIVATE ${TIFF_INCLUDE_DIRS})
        target_link_libraries(${name} PRIVATE ${TIFF_LIBRARIES})
    elseif(TARGET PkgConfig::TIFF)
        target_link_libraries(${name} PRIVATE PkgConfig::TIFF)
    endif()

    # Link OpenMP
    if(OpenMP_CXX_FOUND)
        target_link_libraries(${name} PRIVATE OpenMP::OpenMP_CXX)
    endif()

    # Link math
    target_link_libraries(${name} PRIVATE m)
endfunction()

# Executables
add_s1_executable(s1_print        src/print_main.cpp       "")
add_s1_executable(s1_write_tiff   src/image_write_main.cpp src/image_write.cpp)
add_s1_executable(s1_plot         src/plot_main.cpp        src/plot.cpp)

# Optional: install support
install(TARGETS s1_print s1_write_tiff s1_plot RUNTIME DESTINATION bin)
