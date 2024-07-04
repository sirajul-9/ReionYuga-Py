from cffi import FFI
import os
ffi = FFI()

ffi.set_source(
    "funcs_cffi",  # Name of the generated module
    """
    #include "source.h"
    """,
    sources=["source.c"],
    libraries=["m", "fftw3f", "fftw3f_omp"],
    extra_compile_args=['-w','-fopenmp', '-std=c99']
)

ffi.cdef("""
void read_output(char *fname, int read_flag, long *box, float *rra, float *vva, float *aa);
void cic_vmass(float *ro_dum,float *data,long tot_particles, long N1,long N2,long N3,int xin, int yin, int zin, int min, int col);
void parallelize(int Nthreads);
void read_fof(char *fname, int read_flag, long *totcluster, float *halo, float *aa);
void cic_vmass_parallel(float *ro_dum, float *data, long tot_particles, long N1, long N2, long N3, int xin, int yin, int zin, int min, int col, int num);
void get_nhs(float *nh, float *ng, float *nx,long N1, long N2, long N3,float LL, float nion, float rmfp);
void density_2_mass(float *ro_dum,float *data,long MM, long N1,long N2,long N3,int xin, int yin, int zin, int col, float *mass);
void calpow_mom(float *ro,int Nbin,double* power, double* kmode, double* power_P2,double* power_P4, double *no, long N1, long N2, long N3, float LL);
""")

if __name__ == "__main__":
    ffi.compile(verbose=True)

