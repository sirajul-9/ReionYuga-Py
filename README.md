# ReionYuga-Py

This repository contains a Python project that integrates and utilizes the C code, ReionYuga. The original code is located at https://github.com/rajeshmondal18/ReionYuga. The primary purpose of this integration is to leverage the optimized C functions provided by the original codebase within a Python environment, combining the performance benefits of the C code with the flexibility and ease of use of Python. The C code uses shared memory parallelization by OpenMP for faster execution.

For running this code, the outputs of N-body code (https://github.com/rajeshmondal18/N-body) and Friends-of-Friend Halo Finder(https://github.com/rajeshmondal18/FoF-Halo-finder) are needed.



# Required Libraries

CFFI: can be installed with <pre>  pip install cffi </pre>
for local installation (without root access) you may use <pre> pip install --user cffi </pre> or <pre> pip3 install --user cffi </pre>
FFTW-3 library needs to be installed in the system. It can be downloaded with 
<pre> sudo apt-get install libfftw3-dev </pre>

If that does not work download FFTW-3 from http://www.fftw.org/fftw3_doc/Installation-on-Unix.html#Installation-on-Unix Then install it with the flags --enable-float, --enable-shared and --enable-openmp as written below

After extracting the downloaded file go to the directory and do the following.
<pre> 
./configure --enable-float --enable-shared --enable-openmp
 make
 sudo make install
</pre>

For installing locally without root access:
<pre> 
./configure --enable-float --enable-shared --enable-openmp --prefix=/path/to/your/local/install/directory
 make
 make install
</pre>
 

# Instructions for Running the Code

Download all the files in this repository (CFFI branch) as ZIP and extract. Or use the following command:
<pre>
 git clone https://github.com/sirajul-9/ReionYuga-Py
</pre>

If you have installed FFTW locally, edit the "build_funcs.py" file. Add the following arguments to the ffi.set_source() function:
<pre>
include_dirs=['/path/to/local/include/directory'],

library_dirs=['/path/to/local/lib/directory'],
</pre>

For example
<pre>
 ffi.set_source(
    "funcs_cffi",  
    """
    #include "source.h"
    """,
    sources=["source.c"],
    libraries=["m", "fftw3f", "fftw3f_omp"],
    extra_compile_args=['-w','-fopenmp', '-std=c99'],
    library_dirs=["/home/sirajul/local_install/lib"],
    include_dirs=["/home/sirajul/local_install/include/"]
)
</pre>
    
Then run the file build_funcs.py with 
<pre>
python3 build_funcs.py
</pre>
If you are running with locally installed libraries run the following command (or write it in your job submission script):
<pre> 
LD_LIBRARY_PATH=/path/to/local/lib/directory
export LD_LIBRARY_PATH
</pre>
Then the main code can be run using:
<pre>
python3 ionz_main.py 7 8 9 10 11 13
</pre>

Here the numbers after the ionz_main.py are values of redshifts. Before running it, open the ionz_main.py file and update the nbody and halo-catalogue files path as directed in the comments in the code.

You can set the number of OpenMP threads too. Experiment with it and see for which value it is giving maximum speed-up. Making it too large for small array dimesnions causes large overhead, hence worse performance.

To get 21-cm Brightness Temperature map, multiply the output of this code with $\bar T/\bar{\rho}$.

Where $$\bar{\rho} = \frac{\text{Total number of DM particles}}{N1 \times N2 \times N3 / \text{sfac}^3}$$ and $$\bar T = \frac{27}{\sqrt 10} \sqrt{1+z}$$  z = redshift, N1,N2 and N3 are grid dimensions of the N-body box.
