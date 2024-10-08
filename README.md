# sentinel1_decode
For additional information on Level-0 product decoding, see the [SAR Space Packet Protocol Data Unit Specification](https://sentinels.copernicus.eu/documents/247904/2142675/Sentinel-1-SAR-Space-Packet-Protocol-Data-Unit.pdf) and the [Sentinel-1 Level-0 Data Decoding Package](https://sentinel.esa.int/documents/247904/0/Sentinel-1-Level-0-Data-Decoding-Package.pdf/a8742c59-4914-40c4-8309-c77515649f17).

## Table of Contents
* [Introduction](#introduction)
* [Commands and Examples](#commands)
  * [Image Writing](#writing-images)
  * [Plotting](#plotting)
  * [Packet Information and Performance Testing](#packet-information-and-performance-testing)
* [Results with Point Targets](#results-with-point-targets)
* [Compiling](#compiling)
   * [Build Script](#build-script)
   * [Dependencies](#dependencies)
## Introduction
**sentinel1_decode** is a C++ program/library for quickly decoding Level-0 Raw data from the Sentinel-1 satellite. I am also working on implementing most of the [Level-1 Algorithm](#results-with-point-targets). I am creating this as an education experience for myself, and because there isn't a public fast, simple, and robust program for decoding this data. Currently, I can decode the complex samples for all packets in a data file in approximately 30 seconds (on a Ryzen 5800). 

## Commands
### Writing Images
```bash
$ bin/write --help
burst [swath] [burst_num] [in_path] [out_path]
swath [swath] [in_path] [out_path]
burst_replica_chirps [swath] [burst_num] [in_path] [out_path]
swath_replica_chirps [swath] [in_path] [out_path]
range_compressed_burst [swath] [burst_num] [in_path] [out_path]
range_compressed_swath [swath] [in_path] [out_path]
Scaling Options: [--norm_log|--norm|--mag|--real|--imag]
```
#### Examples
The sample image *data/points/point.dat* is the VV data from [S1A_IW_RAW__0SDV_20240813T095440_20240813T095513_055193_06BA22_1119-RAW](https://search.asf.alaska.edu/#/?searchType=List%20Search&searchList=S1A_IW_RAW__0SDV_20240813T095440_20240813T095513_055193_06BA22_1119-RAW&resultsLoaded=true&granule=S1A_IW_RAW__0SDV_20240813T095440_20240813T095513_055193_06BA22_1119-RAW). The colors were added via *Singleband psuedocolor* with *Cumulative cut count* in QGIS.
```bash
$ bin/write swath IW2 data/points/point.dat IW2.tif --norm
```
![swath_write_example](imgs/swath_write_example.png)
```bash
$ bin/write range_compressed_swath IW2 data/points/point.dat point_targets.tif --norm
```
![swath_rc_write_example](imgs/swath_rc_write_examples.png)
```bash
$ bin/write burst_replica_chirps IW2 5 data/points/point.dat burst_replica_chirps.tif
```
![burst_replica_chirps_example](imgs/burst_replica_chirp_example.png)
### Plotting
```bash
$ bin/plot --help
signal [packet_index] [mode] [path]
swath [swath] [path]
burst [swath] [burst_num] [path]
fft [packet_index] [fft_size] [path] [--inverse]
fft2 [swath] [burst_num] [path] [fft_rows] [fft_cols] [--inverse]
fft_axis [swath] [burst_num] [axis] [fft_size] [path] [--inverse]
range_compressed_burst [swath] [burst_num] [path]
range_compressed_swath [swath] [path]
Scaling Options: [--norm_log|--norm|--mag|--real|--imag]
```
#### Examples
The sample image *data/sample/sample.dat* is the VV data from [S1A_IW_RAW__0SDV_20240806T135224_20240806T135256_055093_06B68A_AE41](https://search.asf.alaska.edu/#/?searchType=List%20Search&searchList=S1A_IW_RAW__0SDV_20240806T135224_20240806T135256_055093_06B68A_AE41&resultsLoaded=true&granule=S1A_IW_RAW__0SDV_20240806T135224_20240806T135256_055093_06B68A_AE41-RAW)
```
$ bin/plot signal 0 sample_data/sample.dat
```
![plot_command_example](imgs/complex_sample_plot_example.png)
```bash
$ bin/plot swath IW3 data/sample/sample.dat --norm
```
![plot_swath_example](imgs/iw3_swath.png)
```bash
$ bin/plot range_compressed_burst IW1 2 data/sample/sample.dat --norm
```
![range_compression_example](imgs/range_compressed.png)
```bash
$ bin/plot range_compressed_swath IW1 data/sample/sample.dat --norm
```
![range_compression_example](imgs/range_compressed_swath.png)

#### Packet Information and Performance Testing
```bash
$ bin/main --help
print_packet_info [packet_index] [path]
print_complex_samples [packet_index] [path]
print_swath_names [path]
print_index_records [path]
print_annotation_record [record_index] [path]
time [num_packets] [path]
thread_test [path]
omp_test [path]
```
#### Examples
The sample image *data/sample/sample.dat* is the VV data from [S1A_IW_RAW__0SDV_20240806T135224_20240806T135256_055093_06B68A_AE41](https://search.asf.alaska.edu/#/?searchType=List%20Search&searchList=S1A_IW_RAW__0SDV_20240806T135224_20240806T135256_055093_06B68A_AE41&resultsLoaded=true&granule=S1A_IW_RAW__0SDV_20240806T135224_20240806T135256_055093_06B68A_AE41-RAW)
```bash
$ bin/main print_packet_info 0 data/sample/sample.dat
Primary Header:
packet_version_number: 0
packet_type: 0
secondary_header_flag: 1
process_id: 65
process_category: 12
sequence_flags: 3
packet_sequence_count: 11157
packet_data_length: 18621

Secondary Header:
coarse_time: 1406987562
fine_time: 7757
sync_marker: 892270675
data_take_id: 225252800
ecc_number: 8
test_mode: 0
rx_channel_id: 0
instrument_configuration_id: 7
sc_data_word_index: 27
sc_data_word: 48835
space_packet_count: 240533
pri_count: 243384
error_flag: 0
baq_mode: 12
baq_block_length: 31
range_decimation: 8
rx_gain: 8
tx_ramp_rate_sign: 1
tx_ramp_rate_mag: 1605
pulse_start_frequency_sign: 0
pulse_start_frequency_mag: 12335
pulse_length: 1967
rank: 9
pri: 21859
swst: 3681
swl: 13979
ssb_flag: 0
polarisation: 7
temperature_compensation: 3
elevation_beam_address: 6
azimuth_beam_address: 385
calibration_mode: 0
tx_pulse_number: 6
signal_type: 0
swap: 1
swath_number: 10
num_quadratures: 11938

Operating Mode Info:
Data Format: D
BAQ Mode: fdbaq_mode_0
BAQ Block Length: 256
Test Mode: measurement_mode
Sensor Mode: interferomatric_wide_swath
Signal Type: echo
Error Status: nominal

Pulse Info:
Swath: IW1
RX Polarization: V
TX Polarization: V
Pulse Length: 52.4048
TX Ramp Rate (TXPRR): 1.07823
Start Frequency (TXPSF): -28.2515
PRI: 582.367
SWL: 372.428
SWST: 98.0692
RX Gain: -4
Range Decimation: 8
TX Pulse Number: 6
```

```bash
$ bin/main print_complex_samples 0 data/sample/sample.dat
# 50,000 lines later...
...
(-11.201,1.59988)
(1.59988,8.00034)
(-4.80058,-11.201)
(-4.80058,20.8031)
(4.80058,4.80058)
(17.6024,17.6024)
(4.80058,31.7194)
(-14.4017,8.00034)
(-11.201,20.8031)
(-4.80058,1.59988)
(1.59988,-4.80058)
(-8.00034,-4.80058)
(-4.80058,-14.4017)
(1.59988,11.201)
(1.59988,-4.80058)
```

## Results with Point Targets
Currently, basic range compression is the maximum level of processing that I have implemented. Although, I am working on implementing much of the [Level-1 Algorithm](https://sentinel.esa.int/documents/247904/1877131/Sentinel-1-Level-1-Detailed-Algorithm-Definition). Here is an example of range compression using some ships, outside of Shanghai, as point targets:
![shanghai_ships](imgs/shanghai_ships_v2.png)
```bash
$ bin/write range_compressed_swath IW2 data/points/point.dat point_targets.tif --norm
```
![shanghai_range_compression](imgs/shanghai_range_compression_v2.png)
The data is from the VV portion of this product: [S1A_IW_RAW__0SDV_20240813T095440_20240813T095513_055193_06BA22_1119-RAW](https://search.asf.alaska.edu/#/?searchType=List%20Search&searchList=S1A_IW_RAW__0SDV_20240813T095440_20240813T095513_055193_06BA22_1119-RAW&resultsLoaded=true&granule=S1A_IW_RAW__0SDV_20240813T095440_20240813T095513_055193_06BA22_1119-RAW).

## Compiling
### Build Script
At the moment, compiling is done via *build.sh*. You may have to edit include paths to match your system. To run it, run the following commands:
```bash
$ chmod +x build.sh
$ ./build.sh
```
Currently, there are three basic CLI utilities that get compiled:
 * **bin/main** for printing packet, index, and annotation info, and for performance timing.
 * **bin/write** for writing *tiff* images from raw and processed data.
 * **bin/plot** for plotting raw and processed data.

### Dependencies

My goal is to not involve too many external dependancies; however, some are necessary. For image writing and plotting, *[OpenMP](https://curc.readthedocs.io/en/latest/programming/OpenMP-C.html)* and *[FFTW3](https://www.fftw.org/)* are required. The image writing is done via *[libtiff](http://www.libtiff.org/)*. The plotting functionality also uses *[matplotlib-cpp](https://github.com/lava/matplotlib-cpp)*, which is included in the *include* directory. It is dependant on being linked to Python and Numpy. For me, they are located in my conda (miniforge3) environment. The *build.sh* file handles this for me on my MacOS and Linux setups, but may need to be modified to match your installs. 