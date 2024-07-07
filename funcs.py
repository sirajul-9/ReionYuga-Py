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

def set_num_threads(nthreads):
    n=ffi.cast("int", nthreads)
    lib.parallelize(n)

def get_box_info(fname):
    # Convert Python objects to CFFI objects
    fname_c = ffi.new("char[]", fname.encode('utf-8'))
    box=np.zeros(4, dtype=np.int64)
    la = np.zeros(7, dtype=np.float32)
    dummy = ffi.new("float *")
    box_c = ffi.cast("long*", box.ctypes.data)
    #rra_c = ffi.cast("float*", rra.ctypes.data)
    #vva_c = ffi.cast("float*", vva.ctypes.data)
    la_c = ffi.cast("float*", la.ctypes.data)

    # Call the C function
    lib.read_output(fname_c, 1, box_c, dummy, dummy, la_c)
    
    p = density_params(np.float32(la[1]), np.float32(la[2]), np.float32(la[3]))
    g = grid_info(np.array([box[0], box[1], box[2]]).astype(np.int64), np.float32(la[0]))
    obj = nbody_info(g, p, np.float32(la[6]), np.int64(box[3]), np.float32(la[5]))
    return obj


def read_nbody_output(fname):
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

#function which returns tot_cluster , halo and vaa
def read_halo_catalogue(fname):
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

#function for calculating power spectrum
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
 
    P0_ptr=ffi.cast("double *", P0.ctypes.data)
    P2_ptr=ffi.cast("double *", P2.ctypes.data)
    P4_ptr=ffi.cast("double *", P4.ctypes.data)
    no_ptr=ffi.cast("double *", no.ctypes.data)
    kmode_ptr=ffi.cast("double *", kmode.ctypes.data)
    ro_ptr=ffi.cast("float *", ro.ctypes.data)

    lib.calpow_mom(ro_ptr,Nbin,P0_ptr,kmode_ptr,P2_ptr,P4_ptr,no_ptr, N1, N2, N3, grid_spacing)
    
    with open(filename, 'w') as out_file:
        for ii in range(Nbin):
            out_file.write(f"{kmode[ii]:e} {P0[ii]:e} {P2[ii]:e} {P4[ii]:e} {int(no[ii]):d}\n")






   
