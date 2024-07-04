# ReionYuga-CFFI
Python Code for simulation of EoR neutral hydrogen density field. This reository uses the CFFI library for passing data between C and Python

# Required Libraries

CFFI: can be installed with pip install cffi

FFTW-3 library needs to be installed in the system. Download FFTW-3 from http://www.fftw.org/fftw3_doc/Installation-on-Unix.html#Installation-on-Unix Then install it with the flags --enable-float, --enable-shared and --enable-openmp

# Instructions for Running the Code

Download all the files in this repository as ZIP and extract.

If you have installed FFTW locally, edit the build_funcs.py file. Add the following arguments to the ffi.set_source() function:
<pre>
include_dirs=['/path/to/local/include/directory'],

library_dirs=['/path/to/local/lib/directory'],
</pre>
    
Then run the file build_funcs.py with 
<pre>
python3 build_funcs.py
</pre>
Then the main code can be run using:
<pre>
python3 ionz_main.py 7 8 9 10 11 13
</pre>

Here the numbers after the ionz_main.py are values of redshifts. Before running it, open the ionz_main.py file and update the nbody and halo-catalogue files path as directed in the comments in the code.

You can set the number of OpenMP threads too. Experiment with it and see for which value it is giving maximum speed-up. Making it too large for small array dimesnions causes large overhead, hence worse performance.

There is seperate flag and threads variable for running the CIC function. Making the flag True runs the CIC parallel and takes more memory. Also more CIC threads take more memory. Additional memory taken for running it parallel = (CIC_threads * N1 * N2 * N3 / sfac^3) *4/1024^3 GB. Decrease it if available memory is low.


To get 21-cm Brightness Temperature map, multiply the output of this code with $\bar T/\bar{\rho}$.
And to get 21-cm Power Spectrum, multiply the power spectrum with $(\bar T/\bar{\rho})^2$.
Where $\bar{\rho} = ({\mathrm Total \; number\; of\; DM\; particles}) / (N1\times N2 \times N3)$
