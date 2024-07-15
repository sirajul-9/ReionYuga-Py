import numpy as np
from funcs_cffi import ffi, lib

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

#function for setting number of OpenMP threads
def set_num_threads(nthreads):
    """
    sets the number of OpenMP threads to nthreads for shared memory parllelism

    Parameters
    ----------
    nthreads: int
             Number of OpenMP threads
    """
    n=ffi.cast("int", nthreads)
    lib.parallelize(n)

#function for reading grid spacing,N1,N2,N3, total number of dark matter particles etc
def get_box_info(fname):
    """
    Reads the grid dimensions, total dark matter, and other necessary parameters from a specified N-body simulation file and returns an instance of the `nbody_info` class.
    
    Parameters
    ----------
    fname : str
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
    >>> filename = "output.nbody_8.000"
    >>> nbody_data = get_box_info(filename)
    """
    fname_c = ffi.new("char[]", fname.encode('utf-8'))
    box=np.zeros(4, dtype=np.int64)
    la = np.zeros(7, dtype=np.float32)
    dummy = ffi.new("float *")
    box_c = ffi.cast("long*", box.ctypes.data)
    #rra_c = ffi.cast("float*", rra.ctypes.data)
    #vva_c = ffi.cast("float*", vva.ctypes.data)
    la_c = ffi.cast("float*", la.ctypes.data)

    
    lib.read_output(fname_c, 1, box_c, dummy, dummy, la_c)
    
    p = density_params(np.float32(la[1]), np.float32(la[2]), np.float32(la[3]))
    g = grid_info(np.array([box[0], box[1], box[2]]).astype(np.int64), np.float32(la[0]))
    obj = nbody_info(g, p, np.float32(la[6]), np.int64(box[3]), np.float32(la[5]))
    return obj

#function for reading position and velocity of dark matter particles
def read_nbody_output(fname):
    """
    Reads the position, velocity of dark matter particles from an N-body simulation output file.

    Parameters
    ----------
    fname : str
        The filename of the N-body output data file.

    Returns
    -------
    list
        A list containing the position and velocity arrays:
        - [position, velocity]
            - position : numpy.ndarray
                The position array of shape (total_particles, 3).
            - velocity : numpy.ndarray
                The velocity array of shape (total_particles, 3).

    Description
    -----------
    The `read_nbody_output` function reads the position and velocity data from the N-body simulation output file.

    Examples
    --------
    >>> filename = "output.nbody_8.000"
    >>> position, velocity = read_nbody_output(filename)
    >>> print(position.shape)
    (1000, 3)
    >>> print(velocity.shape)
    (1000, 3)
    """
    fname_c = ffi.new("char[]", fname.encode('utf-8'))
    obj = get_box_info(fname)
    
    box=np.zeros(4, dtype=np.int64)
    tot_DM=obj.tot_DM
    rra = np.zeros([tot_DM,3],dtype=np.float32)
    vva = np.zeros([tot_DM,3],dtype=np.float32)
    la = np.zeros(7, dtype=np.float32)

    r_ptr=ffi.cast("float *", rra.ctypes.data)
    v_ptr=ffi.cast("float*", vva.ctypes.data)
    box_c = ffi.cast("long*", box.ctypes.data)
    la_c = ffi.cast("float*", la.ctypes.data)
    lib.read_output(fname_c,2,box_c,r_ptr,v_ptr,la_c)
    return [rra,vva]  

#function which returns halo mass, position and velocity
def read_halo_catalogue(fname):
    """
    Reads the halo catalogue data from a specified file.

    Parameters
    ----------
    filename : str
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
    >>> filename = "halo_catalogue_8.000"
    >>> halo_catalogue = read_halo_catalogue(filename)
    >>> print(halo_catalogue.shape)
    (1000, 7)
    """
    fname_c = ffi.new("char[]", fname.encode('utf-8'))
    tot_cluster = ffi.new("long *")
    dummy = ffi.new("float *")
    vaa = ffi.new("float *")
    
    lib.read_fof(fname_c, 1, tot_cluster, dummy, vaa)
    
    total_clusters = tot_cluster[0]
    halo = np.zeros([total_clusters, 7], dtype=np.float32)
    halo_ptr=ffi.cast("float *", halo.ctypes.data)    

    lib.read_fof(fname_c, 2, tot_cluster, halo_ptr, vaa)
    
    return halo


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
    return np.float32(H)


def cic_vmass(data,out_arr_dims,xindex,yindex,zindex,mindex):
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

    Returns
    -------
    numpy.ndarray
        The output array containing the CIC density at each grid point.

    Description
    -----------
    The cic_vmass function computes the Cloud-in-Cell (CIC) mass assignment for a given dataset. 
    
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

    ro_ptr=ffi.cast("float *", ro.ctypes.data)
    data_ptr=ffi.cast("float *", data.ctypes.data) 
    lib.cic_vmass(ro_ptr,data_ptr,tot_particles,N1,N2,N3,xindex,yindex,zindex,mindex,col)
    return ro


def get_ion(nh,ngamma,grid_spacing,nion,rmfp):
    """
    Calculates the ionization density based on given parameters. 

    Parameters
    ----------
    nh : numpy.ndarray
        The neutral hydrogen density array of shape (N1, N2, N3).
    ngamma : numpy.ndarray
        The halo density array of shape (N1, N2, N3).
    grid_spacing : float
        The grid spacing 
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
    ionizing photon density (ngamma)and the reionization model parameters (nion and rmfp).It implements the Excursion Set Formalism.
    See documentation for more info.

    Examples
    --------
    >>> nh = np.random.rand(3, 3, 3) * 10  # Random hydrogen density
    >>> ngamma = np.random.rand(3, 3, 3) * 5  # Random halo density
    >>> grid_spacing = 0.1
    >>> nion = 0.5
    >>> rmfp = 1.0
    >>> ion_density = get_ion(nh, ngamma, grid_spacing, nion, rmfp)
    >>> print(ion_density)
    """
    N1,N2,N3=nh.shape
    nxion=np.zeros([N1,N2,N3],dtype=np.float32)
    nh_ptr=ffi.cast("float *", nh.ctypes.data)
    ng_ptr=ffi.cast("float *", ngamma.ctypes.data)
    nxion_ptr=ffi.cast("float *", nxion.ctypes.data)
    lib.get_nhs(nh_ptr,ng_ptr,nxion_ptr,N1,N2,N3,grid_spacing,nion,rmfp)
    return nxion


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
    ro_ptr=ffi.cast("float *", ro.ctypes.data)
    data_ptr=ffi.cast("float *", data.ctypes.data)
    mass_ptr=ffi.cast("float *", mass.ctypes.data)
    lib.density_2_mass(ro_ptr,data_ptr,tot_particles,N1,N2,N3,xindex,yindex,zindex,col,mass_ptr)
    return mass







   
