import numpy as np
from funcs import *
import time
import gc
import sys

#Initialise required variables
#**********************************************************#
nbody_path="/home/sirajul/packages/N-body"                 #write the n-body output directory path here, don't put slash at the end
halo_cat_path='/home/sirajul/packages/FoF-Halo-finder'     #write the halo_catalogue directory path here, don't put slash at the end
sfac = 2                                                   #factor for rescaling the  grid
Nthreads=2                                                 #set the number of openmp threads
cic_parallel_flag=True                                     #run CIC in parallel or not
cic_threads=2                                              #number of threads for CIC

#parameters of the reionization model
nion=23.21                                                 #parameter Nion, fid. value=23.21
rmfp=20.00                                                 #parameter Rmfp, fid. value=20.00
mmin=10                                                    #parameter Mmin(in units of DM_mass), fid. value=1.09 * 10**9 M_sun
#**********************************************************#



#read redshifts from command line args
#**********************************************************#
if len(sys.argv)<2:
    print("Please enter redshift values")
    sys.exit(0)

redshifts=np.array(sys.argv[1:]).astype(np.float32)
print("\nGenerating HI maps at {} redshifts".format(len(redshifts)))
#**********************************************************#


start=time.time()
set_num_threads(Nthreads)                                  #setting openmp parallelsation with Nthreads threads
if not os.path.exists("ionz_out"):
        os.makedirs("ionz_out")
print("****************************************************")


#starting redshift loop
for z in redshifts:
    print("Redshift = {:.3f}".format(z))

    #reading nbody_outputs
    #******************************************************#
    print("Reading N-body output file")
    filename1 = nbody_path+"/output.nbody_{:.3f}".format(z)      
    filename1=filename1.encode('utf-8')

    box = get_box_info(filename1)
    N1,N2,N3=box.grid.dim
    LL=box.grid.grid_spacing
    tot_DM=box.tot_DM
    scale_fac=box.scale_factor

    DM_data,DM_vel = read_nbody_output(filename1)          #storing position and velocity of DM particles
    print("N1 = {} N2 = {} N3 = {} LL={:.4f} and total DM = {}".format(N1,N2,N3,LL,tot_DM))
    #******************************************************#


    #redefining grid
    N1=int(N1/sfac)  
    N2=int(N2/sfac) 
    N3=int(N3/sfac)  
    LL=LL*sfac; 


    robar=tot_DM/(1.*N1*N2*N3)                             #mean number density (grid)^{-3}
    vfac =1./(hubble_param(scale_fac)*scale_fac*scale_fac) #for redshift space distortion



    #storing all the required DM data in single 2D array
    #******************************************************#
    #storing redshift space coordinates along LoS(Z axis) in a single array and scaling the data
    vel_z = DM_vel[:, 2]                                   
    del DM_vel
    DM_data = np.hstack((DM_data, np.zeros((tot_DM, 2), dtype=np.float32)))    #two more columns for RS co-ordinates and mass
    DM_data[:, 3] = (DM_data[:, 2] + vfac * vel_z)/sfac    #applying RSD along LoS (Z direction)
    del vel_z
    DM_data[:,0:3]/=sfac                                   #scaling all the four co-ordinates
    gc.collect()                                           #collect garbages

    #applying periodic boundary condition to the RS space co-ordinates 
    DM_data[:,3] += N3
    DM_data[:,3]=DM_data[:,3] - 1.0 * N3 * np.trunc(np.floor(DM_data[:,3])/N3)  
  
    DM_data[:,4]=1                                         #taking unit mass of all DM particles
    DM_data=np.round(DM_data,6)
    #******************************************************#

    #reading halo_catalogues
    #******************************************************#
    print("Reading halo catalogues")
    filename2 =halo_cat_path+"/halo_catalogue_{:.3f}".format(z)      
    filename2=filename2.encode('utf-8')
    halo_data = read_halo_catalogue(filename2)             #getting halo position, velocity and mass data
    halo_data=halo_data[halo_data[:,0]>=mmin]              #keeping halos with mass >= mmin only
    gc.collect()
    tot_cluster=halo_data.shape[0]
    print("Total number of halos with mass greater than {} = {}".format(mmin,tot_cluster))
    print("Minimum Halo mass = ", halo_data[:,0].min())

    #redefining halo
    halo_data[:,1:4] /= sfac
    halo_data=np.round(halo_data,6)
    #******************************************************#    


    #cic on DM_data and halo_data to get nh and ngamma respectively
    #******************************************************#
    dimensions=[N1,N2,N3]

    ngamma = cic_vmass(halo_data,dimensions,xindex=1,yindex=2,zindex=3,mindex=0)
    del halo_data                                          
    gc.collect()

    start_cic=time.time()
    nh= cic_vmass(DM_data,dimensions,xindex=0,yindex=1,zindex=2,mindex=4,run_parallel=cic_parallel_flag,num_threads=cic_threads)
    print("CIC done in {:.2f} sec".format(time.time()-start_cic))
    #******************************************************#


    #getting ionisation field using excursion set formalism
    #******************************************************#
    nxion=get_ion(nh=nh,ngamma=ngamma,grid_spacing=LL,nion=nion,rmfp=rmfp)
    #******************************************************#


    #getting neutral hydrogen density field in real space
    #******************************************************#
    np.maximum(0, 1 - nxion, out=nxion)                    #neutral hydrogen frac = 1-nxion; 0 if it comes out to be negative
    nh *= nxion                                            #neutral hydrogen density in grid^{-3}
    print("volume averaged neutral hydrgen fraction(for real-space data) = {:.6f}".format(np.average(nxion)))
    print("Mass averaged Neutral Hydrogen Fraction (for data in real space)= {:.6f}".format(np.average(nh)/robar))
    
    #storing real-space data into file
    filename3="ionz_out/HI_map_{:.3f}".format(z)
    with open(filename3, 'wb') as file:
        file.write(np.array([N1, N2, N3], dtype=np.int32).tobytes())
        file.write(nh.tobytes())
    #******************************************************#
    

    #getting neutral hydrogen density field in redshift space
    #******************************************************#
    print("mapping to redshift space")
    DM_data[:,4]=density_to_mass(nxion, DM_data, xindex=0, yindex=1, zindex=2)  #getting HI mass at particle pos. in RS space
    nh=cic_vmass(DM_data,dimensions,xindex=0,yindex=1,zindex=2,mindex=4,run_parallel=cic_parallel_flag,num_threads=cic_threads)
                                                                                #geiing HI density at grids 

    #storing redshift-space data into file
    filename4="ionz_out/HI_maprs_{:.3f}".format(z)
    with open(filename4, 'wb') as file:
        file.write(np.array([N1, N2, N3], dtype=np.int32).tobytes())
        file.write(nh.tobytes())

    
    print("Mass averaged Neutral Hydrogen Fraction (for data in redshift space) = {:.6f}".format(np.average(nh)/robar))
    #******************************************************#
    print("-------------------------------------------------")
    gc.collect()
print(*****************************************************)
tot_time=time.time()-start
hr=int(tot_time/3600)
minute=int((tot_time-int(tot_time/3600)*3600)/60)
sec=tot_time-int(tot_time/60)*60
print("Total time taken={} hr {} min {} sec".format(hr,minute,int(sec)))

