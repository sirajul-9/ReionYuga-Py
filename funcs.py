from ctypes import *
import numpy as np
import os
current_directory = os.getcwd()
lib=CDLL(current_directory+'/Reion.so')

read_output=lib.read_output
read_output.argtypes=[c_char_p,c_int,POINTER(c_long),POINTER(c_float),POINTER(c_float),POINTER(c_float)]
read_output.restype=None

read_fof=lib.read_fof
read_fof.argtypes=[c_char_p,c_int,POINTER(c_long),POINTER(c_float),POINTER(c_float)]
read_fof.restype=None

#function for openmp parallelization
parallel=lib.parallelize
parallel.argtypes=[c_int]
parallel.restype=None

def set_num_threads(nthreads):
    n=c_int(nthreads)
    parallel(n)


class density_params:
    def __init__(self, omega_m, omega_lam, omega_b):
        self.omega_m = omega_m
        self.omega_lam = omega_lam
        self.omega_b = omega_b

class grid_info():
    def __init__(self,dim,LL):
        self.dim=dim
        self.grid_spacing=LL

class nbody_info:
    def __init__(self,grid,params,vaa,tot_DM,DM_m):
        self.grid=grid
        self.params=params
        self.scale_factor=vaa
        self.tot_DM=tot_DM
        self.DM_mass=DM_m



#function for reading grid spacing,N1,N2,N3 and tot_DM
def get_box_info(fname):
    """
    Reads the grid dimensions, total dark matter, and other necessary parameters from a specified N-body simulation file and returns an instance of the `nbody_info` class.
    
    Parameters
    ----------
    fname : bytes
        The filename of the N-body data file.
    
    Returns
    -------
    nbody_info
        An object of the `nbody_info` class containing:
        - grid : grid_info
            The computational grid information(dimensions and spacing).
        - params : density_params
            The cosmological density parameters.
        - scale_factor : float
            The scale factor of cosmic expansion.
        - tot_DM : int
            The total number of dark matter particles.
        - DM_mass : float
            The mass of individual dark matter particles.
    
    Description
    -----------
    The `get_box_info` function reads data from a specified N-body output file. It extracts the grid dimensions (`N1, N2, N3`), grid-spacing, total dark matter particles (`tot_DM`), and other parameters as written above.
    
    Examples
    --------
    >>> filename = b"nbody_data"
    >>> nbody_data = get_box_info(filename)
    """
    box = (c_long * 4)()
    dummy = c_float()
    la = (c_float * 7)()
    read_output(fname, 1, box, byref(dummy), byref(dummy), la)
    p = density_params(np.float32(la[1]), np.float32(la[2]), np.float32(la[3]))
    g = grid_info(np.array([box[0], box[1], box[2]]).astype(np.int64), np.float32(la[0]))
    obj = nbody_info(g, p, np.float32(la[6]), np.int64(box[3]), np.float32(la[5]))
    return obj


#function for reading position,velocity and vaa
def read_nbody_output(fname):
    """
    Reads the position, velocity of dark matter particles from an N-body simulation output file.

    Parameters
    ----------
    fname : bytes
        The filename of the N-body output data file.

    Returns
    -------
    list
        A list containing the position and velocity arrays:
        - [position, velocity]
            - position : numpy.ndarray
                The rounded position array of shape (total_particles, 3).
            - velocity : numpy.ndarray
                The rounded velocity array of shape (total_particles, 3).

    Description
    -----------
    The `read_nbody_output` function reads the position and velocity data from the N-body simulation output file.

    Examples
    --------
    >>> filename = b"nbody_output.dat"
    >>> position, velocity = read_nbody_output(filename)
    >>> print(position.shape)
    (1000, 3)
    >>> print(velocity.shape)
    (1000, 3)
    """
    obj = get_box_info(fname)
    box=(c_long * 4)()
    tot_DM=obj.tot_DM
    rra = np.zeros([tot_DM,3],dtype=np.float32)
    vva = np.zeros([tot_DM,3],dtype=np.float32)
    la = (c_float * 7)()
    r_ptr=rra.ctypes.data_as(POINTER(c_float))
    v_ptr=vva.ctypes.data_as(POINTER(c_float))
    read_output(fname,2,box,r_ptr,v_ptr,la)
    return [rra,vva]  


#function which returns tot_cluster , halo and vaa
def read_halo_catalogue(filename):
    """
    Reads the halo catalogue data from a specified file.

    Parameters
    ----------
    filename : bytes
        The filename of the halo catalogue data file.

    Returns
    -------
    numpy.ndarray
        The halo catalogue array of shape (total_clusters, 7).
        column 1             ---> mass of a halo
        column 2 to column 4 ---> x,y and z co-ordinates
        column 5 to column 7 ---> x,y and z components of velocity

    Description
    -----------
    The `read_halo_catalogue` function reads the halo catalogue data from a specified file.

    Examples
    --------
    >>> filename = "halo_catalogue.dat"
    >>> halo_catalogue = read_halo_catalogue(filename)
    >>> print(halo_catalogue.shape)
    (1000, 7)
    """
    tot_cluster=c_long()
    dummy = c_float()
    vaa=c_float()
    read_fof(filename,1,byref(tot_cluster),byref(dummy),byref(vaa))
    halo = np.zeros([tot_cluster.value,7],dtype=np.float32)
    read_fof(filename,2,byref(tot_cluster),halo.ctypes.data_as(POINTER(c_float)),byref(vaa))
           
    return halo

#calculate hubble parameter
def hubble_param(scale_fac, omega):
    """
    Calculates Hubble parameter 

    Parameters
    ----------
    scale_fac: float
        scale factor of cosmic expansion
    omega: density_params
        Cosmological density parameters

    Returns
    ---------
    float
        The value of Hubble parameter
    """

    H = np.sqrt(np.round(omega.omega_m * scale_fac**(-3), 6) +
                np.round((1.0 - omega.omega_m - omega.omega_lam) * scale_fac**(-2), 6) +
                np.round(omega.omega_lam, 6))
    return H



#function for CIC
cic_vmass_=lib.cic_vmass
cic_vmass_.argtypes=[POINTER(c_float),POINTER(c_float),c_long,c_long,c_long,c_long,c_int,c_int,c_int,c_int,c_int]   
cic_vmass_.restype=None

cic_vmass_par=lib.cic_vmass_parallel
cic_vmass_par.argtypes=[POINTER(c_float),POINTER(c_float),c_long,c_long,c_long,c_long,c_int,c_int,c_int,c_int,c_int,c_int]   
cic_vmass_par.restype=None

def cic_vmass(data,out_arr_dims,xindex,yindex,zindex,mindex,run_parallel=False,num_threads=1):
    """
    Computes the Cloud-in-Cell (CIC) density at each grid point for a given field distribution.

    Parameters
    ----------
    data : numpy.ndarray
        The input dataset of particles.
    out_arr_dims : tuple
        The dimensions of the output array (N1, N2, N3).
    xindex, yindex, zindex, mindex : int
        Indices for accessing the particle positions and masses in the input dataset.
    run_parallel : bool, optional
        If True, the computation will run in parallel. Defaults to False.
    num_threads : int, optional
        Number of threads to use for parallel computation. Defaults to 1.

    Returns
    -------
    numpy.ndarray
        The output array containing the CIC density at each grid point.

    Description
    -----------
    The cic_vmass function computes the Cloud-in-Cell (CIC) Voronoi mass assignment for a given dataset. 
    It initializes an output array ro with dimensions specified by out_arr_dims. 
    The function then calculates the total number of particles and the number of columns in the input data. 
    If run_parallel is True, the computation is performed in parallel using the specified number of threads (num_threads). 
    Otherwise, the computation runs serially. The output array ro is rounded to 6 decimal places before being returned.

    Examples
    --------
    >>> data = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    >>> out_arr_dims = (3, 3, 3)
    >>> xindex, yindex, zindex, mindex = 0, 1, 2, 3
    >>> result = cic_vmass(data, out_arr_dims, xindex, yindex, zindex, mindex, run_parallel=True, num_threads=4)
    >>> print(result)
    """
    N1,N2,N3=out_arr_dims
    ro = np.zeros([N1,N2,N3],dtype=np.float32)
    tot_particles=data.shape[0]
    col=data.shape[1]
    if run_parallel:
        cic_vmass_par(ro.ctypes.data_as(POINTER(c_float)),data.ctypes.data_as(POINTER(c_float)),c_long(tot_particles),N1,N2,N3,xindex,yindex,zindex,mindex,c_int(col),
        c_int(num_threads))
        return np.round(ro,6)   
    cic_vmass_(ro.ctypes.data_as(POINTER(c_float)),data.ctypes.data_as(POINTER(c_float)),c_long(tot_particles),N1,N2,N3,xindex,yindex,zindex,mindex,c_int(col))
    return ro


#function for getting ionization field
get_ion_=lib.get_nhs
get_ion_.argtypes=[POINTER(c_float),POINTER(c_float),POINTER(c_float),c_long,c_long,c_long,c_float,c_float,c_float]
get_ion_.restype=None

def get_ion(nh,ngamma,grid_spacing,nion,rmfp):
    """
    Calculates the ionization density based on given parameters.

    Parameters
    ----------
    nh : numpy.ndarray
        The neutral hydrogen density array of shape (N1, N2, N3).
    ngamma : numpy.ndarray
        The ionizing photon density array of shape (N1, N2, N3).
    grid_spacing : float
        The grid spacing parameter.
    nion : float
        Parameter N_ion.
    rmfp : float
        Parameter R_mfp.

    Returns
    -------
    numpy.ndarray
        Array of shape (N1, N2, N3) containing ionization fraction at each grid point.

    Description
    -----------
    The get_ion function calculates the ionization density based on the given neutral hydrogen density (nh),
    ionizing photon density (ngamma)and the reionization model parameters (nion and rmfp). See documentation for more info.

    Examples
    --------
    >>> nh = np.random.rand(3, 3, 3) * 10  # Random neutral hydrogen density
    >>> ngamma = np.random.rand(3, 3, 3) * 5  # Random ionizing photon density
    >>> grid_spacing = 0.1
    >>> nion = 0.5
    >>> rmfp = 1.0
    >>> ion_density = get_ion(nh, ngamma, grid_spacing, nion, rmfp)
    >>> print(ion_density)
    """
    N1,N2,N3=nh.shape
    nxion=np.zeros([N1,N2,N3],dtype=np.float32)
    get_ion_(nh.ctypes.data_as(POINTER(c_float)),ngamma.ctypes.data_as(POINTER(c_float)),nxion.ctypes.data_as(POINTER(c_float)),N1,N2,N3,c_float(grid_spacing),
    c_float(nion),c_float(rmfp))
    return nxion

#function for density_2_mass
density_2_mass_=lib.density_2_mass
density_2_mass_.argtypes=[POINTER(c_float),POINTER(c_float),c_long,c_long,c_long,c_long,c_int,c_int,c_int,c_int, POINTER(c_float)]   
density_2_mass_.restype=None

def density_to_mass(ro,data,xindex,yindex,zindex):
    """
    Converts given density field(grid points) into masses at particle's positions'.

    Parameters
    ----------
    ro : numpy.ndarray
        The density array of shape (N1, N2, N3) with dtype float32.
    data : numpy.ndarray
        The input dataset of particles containing the positions.
    xindex, yindex, zindex : int
        Indices for accessing the particle positions in the input dataset.

    Returns
    -------
    numpy.ndarray
        Array containing masses at particle's positions'.

    Description
    -----------
    See documentation for understanding the process.

    Examples
    --------
    >>> ro = np.random.rand(3, 3, 3).astype(np.float32) * 10  # Random density array with dtype float32
    >>> data = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])  # Random input dataset
    >>> xindex, yindex, zindex = 0, 1, 2  # Particle position indices
    >>> mass_result = density_to_mass(ro, data, xindex, yindex, zindex)
    >>> print(mass_result)
    """

    N1,N2,N3=ro.shape
    tot_particles=data.shape[0]
    col=data.shape[1]
    mass=np.zeros(tot_particles,dtype=np.float32)
    density_2_mass_(ro.ctypes.data_as(POINTER(c_float)),data.ctypes.data_as(POINTER(c_float)),c_long(tot_particles),N1,N2,N3,xindex,yindex,zindex,c_int(col),
    mass.ctypes.data_as(POINTER(c_float)))
    return mass


#function for power spectrum
calpow=lib.calpow_mom
calpow.argtypes=[POINTER(c_float),c_int,POINTER(c_double),POINTER(c_double),POINTER(c_double),POINTER(c_double),POINTER(c_double),c_long,c_long,c_long]   
calpow.restype=None

def calculate_power_spec(ro, Nbin, grid_spacing, filename):
    """
    Calculates Power Spectrum multipoles of a given density field. Saves the power spectrum info into a file.
    Parameters
    ----------
    ro: numpy.ndarray
        The density array of shape (N1, N2, N3) with dtype float32.
    Nbin: int
        Number of bins for power spectrum
    grid_spacing: float
        grid spacing of computational grid
    filename: str
        filename for power spectrum output file

    """
    
    N1,N2,N3=ro.shape
    P0=np.zeros(Nbin, dtype=np.float64)
    P2=np.zeros(Nbin, dtype=np.float64)
    P4=np.zeros(Nbin, dtype=np.float64)
    no=np.zeros(Nbin, dtype=np.float64)
    kmode=np.zeros(Nbin, dtype=np.float64)
    calpow(ro.ctypes.data_as(POINTER(c_float)),c_int(Nbin),P0.ctypes.data_as(POINTER(c_double)),kmode.ctypes.data_as(POINTER(c_double)),
           P2.ctypes.data_as(POINTER(c_double)),P4.ctypes.data_as(POINTER(c_double)),no.ctypes.data_as(POINTER(c_double)), N1, N2, N3, c_float(grid_spacing))
    
    with open(filename, 'w') as out_file:
        for ii in range(Nbin):
            out_file.write(f"{kmode[ii]:e} {P0[ii]:e} {P2[ii]:e} {P4[ii]:e} {int(no[ii]):d}\n")




