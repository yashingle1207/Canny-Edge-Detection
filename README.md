Real-Time Canny Edge Detection (C → C++ → SystemC)

This repository demonstrates the iterative development and optimization of a Canny Edge Detector application, starting from a basic C prototype to a C++ implementation with static memory usage, and finally transforming into a SystemC specification model capable of processing video frames in near real-time.
Table of Contents

Table of Contents
  1.	Project Overview
  2.	Repository Structure
  3.	Prerequisites & Dependencies
  4.	Phases of Development
      o	Phase 1: C Prototype
      o	Phase 2: C++ Refined Model
      o	Phase 3: SystemC Modeling
      o	Phase 4: Pipelining & Parallelization
      o	Phase 5: Performance Profiling & Optimization
  5.	Building & Running the Project
      o	C/C++ Versions
      o	SystemC Version
  6.	Testing & Validation
  7.	Results
  8.	Future Improvements
  9.	References


__________________________________ ___________________________ _________________________________

Project Overview
The Canny Edge Detector is a classic algorithm used to detect edges in images by:
  1.	Applying Gaussian smoothing to reduce noise.
  2.	Computing gradients in x and y directions.
  3.	Calculating the magnitude and direction of edges.
  4.	Performing non-maximal suppression to thin out edges.
  5.	Using hysteresis thresholding with high and low thresholds to finalize edge pixels.


In this repository, you will find:
  •	The initial C code (canny_edge.c, hysteresis.c, pgm_io.c).
  •	 A C++ version with static memory allocation and compiler warning fixes.
  •	A SystemC pipeline model that handles high-resolution frames (e.g., 2704x1520) with near real-time performance using sc_fifo channels and parallelized blur stages.


______________________________________________________________________________

Repository Structure
    .
    ├── c_version
    │   ├── canny_edge.c
    │   ├── hysteresis.c
    │   ├── pgm_io.c
    │   └── Makefile
    ├── cpp_version
    │   ├── canny.cpp
    │   ├── Makefile
    │   └── README.md
    ├── systemc_version
    │   ├── canny.cpp
    │   ├── modules/
    │   │   ├── Gaussian_Smooth.cpp
    │   │   ├── Derivative_X_Y.cpp
    │   │   ├── Magnitude_X_Y.cpp
    │   │   ├── Non_Max_Supp.cpp
    │   │   ├── Apply_Hysteresis.cpp
    │   │   └── ...
    │   ├── top.cpp
    │   ├── Makefile
    │   └── README.md
    ├── images
    │   ├── sample_input.pgm
    │   ├── sample_output.pgm
    │   └── ...
    ├── docs
    │   └── project_report.pdf
    └── README.md

_________________________________________________________________________________

Prerequisites & Dependencies

  •	C/C++ Compiler: GCC or Clang for building the initial prototypes.
  •	Make: For convenient build scripts in all phases.
  •	SystemC Library (2.3+): Required for building and running the SystemC model.
  •	Linux Environment: (Optional but recommended) The examples are tested primarily on Linux-based systems.
  •	PGM Image Tools: Utilities like eog or display (optional) for viewing generated PGM files.

__________________________________________________________________________________

Phases of Development

Phase 1: C Prototype
  •	Files: canny_edge.c, hysteresis.c, pgm_io.c.
  •	Goal: Single-image edge detection using a straightforward approach.
  •	Highlights:
      o	Gaussian smoothing with a separable kernel.
      o	First derivatives in x and y with small finite filters.
      o	Non-maximal suppression to thin edges.
      o	Hysteresis with low/high thresholds (tlow, thigh).
      
Phase 2: C++ Refined Model
  •	Files: canny.cpp, minimal headers.
  •	Goal: Address dynamic memory and code warnings, use modern C++ structures.
  •	Highlights:
      o	Static memory for images (no malloc/free).
      o	No compiler warnings (e.g., -Wall, -O2).
      o	Hard-coded parameters (sigma, tlow, thigh) for embedded feasibility.
    
Phase 3: SystemC Modeling
  •	Files: canny.cpp plus separate modules in modules/.
  •	Goal: Convert the C++ code into a SystemC concurrency model.
  •	Highlights:
      o	Stimulus and Monitor modules for reading/writing PGM frames.
      o	Platform (DUT) that uses sc_fifo channels connecting:
          •	Gaussian_Smooth
          •	Derivative_X_Y
          •	Magnitude_X_Y
          •	Non_Max_Supp
          •	Apply_Hysteresis

      o	Each module reads/writes images or short-image buffers.
      
Phase 4: Pipelining & Parallelization
  •	Goal: Enable pipeline concurrency via sc_fifo channels of depth 1.
  •	Highlights:
    o	Break Gaussian smoothing into BlurX and BlurY, each parallelized in multiple SC_THREADs.
    o	Achieve higher throughput by overlapping stage operations on consecutive frames.
    
Phase 5: Performance Profiling & Optimization
  •	Goal: Attain near real-time performance for large frames (e.g., 2704x1520).
  •	Highlights:
    o	Profile code on Linux and Raspberry Pi.
    o	Back-annotate measured latencies (e.g., wait statements) into SystemC modules.
    o	Tweak pipeline stage delays, adjust slice-based parallelism.

_________________________________________________________________________________________


Building & Running the Project
C/C++ Versions

  1.Navigate to the respective folder (c_version or cpp_version).
  2.Build:
      make

Run:
    ./canny input_image.pgm sigma tlow thigh

    e.g. ./canny golfcart.pgm 0.6 0.3 0.8
    

___________________________________________________________________________________

SystemC Version

    1.Ensure SYSTEMC_HOME is set or Makefile is updated with correct SystemC include and lib paths.
    2.Navigate to systemc_version/.
    3.Build:
        make

    Run:
      ./canny

    •	The code may read multiple frames from video/ directory (e.g., Engineering001.pgm to Engineering030.pgm).
    •	Adjust stack size if needed (e.g., ulimit -s 128000).

___________________________________________________________________________________

Testing & Validation
  •	Unit Tests: Intermediate checks on gradient magnitudes, directions, or partial pipeline outputs.
  •	Visual Inspection: Compare generated *_edges.pgm images with reference outputs. Tools:
      o	eog output_edges.pgm
  •	Automated Diff: Compare final edges with reference using diff or ImageDiff.
_____________________________________________________________________________________

Results
  •	High Resolution: Successfully processes frames of size 2704x1520.
  •	Pipeline: Achieves overlapping concurrency, reducing total latency.
  •	Parallel Gaussian Blur: Splitting BlurX/BlurY across multiple threads yields significant throughput gains.
  •	Real-World Timing: Verified on Raspberry Pi for embedded feasibility; back-annotated times reflect realistic hardware behavior.

_______________________________________________________________________________________
  
Future Improvements
  •	GPU Acceleration: Potential to offload heavy computations (BlurX/BlurY) to GPU frameworks (CUDA/OpenCL).
  •	Dynamic Thresholding: Real-time adaptation of sigma, tlow, and thigh to varying scene conditions.
  •	Memory Footprint: Further reduce memory usage by reusing buffers or compressing intermediate representations.

________________________________________________________________________________________

References
  •	Mike Heath: Original canny code structure.
  •	SystemC 2.3: http://accellera.org/downloads/standards/systemc
  •	PBMplus: Tools for .pgm file manipulations.
________________________________________
Maintainer: Yash Ingle
Contact: yashingle1207@gmail.com


