# ReionYuga-Py
Python Code for simulation of EoR neutral hydrogen density field

This repository contains a Python project that integrates and utilizes C code, ReionYuga, using ctypes.
The original code is located at https://github.com/rajeshmondal18/ReionYuga
The primary purpose of this integration is to leverage the optimized C functions provided by the original codebase within a Python environment, combining the performance benefits of the C code with the flexibility and ease of use of Python.

For running this code, the outputs of N-body code (https://github.com/rajeshmondal18/N-body) and Friends-of-Friend Halo Finder(https://github.com/rajeshmondal18/FoF-Halo-finder) are needed.

# Required Library
For running this code, FFTW-3 library needs to be installed in the system. 
Download FFTW-3 from http://www.fftw.org/fftw3_doc/Installation-on-Unix.html#Installation-on-Unix
Then install it with the flags --enable-float and --enable-shared

# Instructions for Running the code
Download all the files in this repository. Before running the code, run the following command:
make shared_lib
This will create a shared object (reion.so). Once the shared object is created, the code can be run with the following command:
python3 ionz_main.py 13 11 10 9 8 7
Here the numbers after the ionz_main.py are values of redshifts.


