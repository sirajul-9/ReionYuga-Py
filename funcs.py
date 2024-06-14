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



#function for reading LL,N1,N2,N3 and tot_DM
def get_box_info(fname):
    """
    input parameter
    fname: str
    filename of the n-body
    retruns object of class nbody_info
    """
    box=(c_long * 4)()
    dummy = c_float()
    la=(c_float * 7)()
    read_output(fname,1,box,byref(dummy),byref(dummy),la)
    p=density_params(np.float32(la[1]),np.float32(la[2]),np.float32(la[3]))
    g=grid_info(np.array([box[0],box[1],box[2]]).astype(np.int64),np.float32(la[0]))
    obj=nbody_info(g,p,np.float32(la[6]),np.int64(box[3]),np.float32(la[5]))
    return obj

#function for reading position,velocity and vaa
def read_nbody_output(fname):
    obj = get_box_info(fname)
    box=(c_long * 4)()
    tot_DM=obj.tot_DM
    rra = np.zeros([tot_DM,3],dtype=np.float32)
    vva = np.zeros([tot_DM,3],dtype=np.float32)
    la = (c_float * 7)()
    r_ptr=rra.ctypes.data_as(POINTER(c_float))
    v_ptr=vva.ctypes.data_as(POINTER(c_float))
    read_output(fname,2,box,r_ptr,v_ptr,la)
    return [np.round(rra,6),np.round(vva,6)]  

#function which returns tot_cluster , halo and vaa
def read_halo_catalogue(filename):
    tot_cluster=c_long()
    dummy = c_float()
    vaa=c_float()
    read_fof(filename,1,byref(tot_cluster),byref(dummy),byref(vaa))
    halo = np.zeros([tot_cluster.value,7],dtype=np.float32)
    read_fof(filename,2,byref(tot_cluster),halo.ctypes.data_as(POINTER(c_float)),byref(vaa))
           
    return [np.int64(tot_cluster.value),np.round(halo,6)]

#function for Hubble parameter
get_Hf = lib.Hf
get_Hf.argtypes=[c_float]
get_Hf.restype=c_float   

def hubble_param(aa):
    aaa=c_float(aa)
    return np.round(get_Hf(aaa.value),6)  


#function for CIC
cic_vmass_=lib.cic_vmass
cic_vmass_.argtypes=[POINTER(c_float),POINTER(c_float),c_long,c_long,c_long,c_long,c_int,c_int,c_int,c_int]   
cic_vmass_.restype=None

def cic_vmass(data,out_arr_dims,xindex,yindex,zindex,mindex):
    N1,N2,N3=out_arr_dims
    ro = np.zeros([N1,N2,N3],dtype=np.float32)
    tot_particles=data.shape[0]
    cic_vmass_(ro.ctypes.data_as(POINTER(c_float)),data.ctypes.data_as(POINTER(c_float)),tot_particles,N1,N2,N3,xindex,yindex,zindex,mindex)
    return np.round(ro,6)


#function for smoothing
get_ion_=lib.get_nhs
get_ion_.argtypes=[POINTER(c_float),POINTER(c_float),POINTER(c_float),c_long,c_long,c_long,c_float,c_float,c_float]
get_ion_.restype=None

def get_ion(nh,ngamma,grid_spacing,nion,rmfp):
    N1,N2,N3=nh.shape
    nxion=np.zeros([N1,N2,N3],dtype=np.float32)
    get_ion_(nh.ctypes.data_as(POINTER(c_float)),ngamma.ctypes.data_as(POINTER(c_float)),nxion.ctypes.data_as(POINTER(c_float)),N1,N2,N3,c_float(grid_spacing),
    c_float(nion),c_float(rmfp))
    return np.round(nxion,6)

#function for density_2_mass
density_2_mass_=lib.density_2_mass
density_2_mass_.argtypes=[POINTER(c_float),POINTER(c_float),c_long,c_long,c_long,c_long,c_int,c_int,c_int,c_int]   
density_2_mass_.restype=None

def density_2_mass(ro,data,out_arr_dims,xindex,yindex,zindex,mindex):
    N1,N2,N3=out_arr_dims
    tot_particles=data.shape[0]
    density_2_mass_(ro.ctypes.data_as(POINTER(c_float)),data.ctypes.data_as(POINTER(c_float)),tot_particles,N1,N2,N3,xindex,yindex,zindex,mindex)

