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
Nthreads=2                                                 #set the number of openmp threads, see for which valule you are getting maximum speedup

cic_parallel_flag=True                                     #run CIC in parallel or not, parallel will require more memory, set it to False if low memory
cic_threads=2                                              #number of threads for CIC, more threads for CIC takes more memory, so keep it <= 2

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
if not os.path.exists("ionz_out"):                         #making directory for storing the outputs
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

    box = get_box_info(filename1)                          #getting info regarding the nbody box
    N1,N2,N3=box.grid.dim                                  #grid diemnsions
    LL=box.grid.grid_spacing                               #grid spacing
    tot_DM=box.tot_DM                                      #total number of dark matter particles
    scale_fac=box.scale_factor                             #scale factor of cosmological expansion

    DM_data,DM_vel = read_nbody_output(filename1)          #storing position and velocity of DM particles
    print("N1 = {} N2 = {} N3 = {} LL={:.4f} and total DM = {}".format(N1,N2,N3,LL,tot_DM))
    #******************************************************#


    #redefining grid
    N1=int(N1/sfac)                                        #reducing number of grid-points along X
    N2=int(N2/sfac)                                        #reducing number of grid-points along Y 
    N3=int(N3/sfac)                                        #reducing number of grid-points along Z 
    LL=LL*sfac;                                            #scaling grid-spacing to retain the same volume of the box


    robar=tot_DM/(1.*N1*N2*N3)                             #mean number density (grid)^{-3}
    vfac =1./(hubble_param(scale_fac)*scale_fac*scale_fac) #for redshift space distortion



    #storing all the required DM data in single 2D array
    #******************************************************#
    #storing redshift space coordinates along LoS(Z axis) in the same array and scaling the data
    vel_z = DM_vel[:, 2]                                   #getting Z component of peculiar velocity
    del DM_vel
    DM_data = np.hstack((DM_data, np.zeros((tot_DM, 2), dtype=np.float32)))    #two more columns for RS co-ordinates and mass
    DM_data[:, 3] = (DM_data[:, 2] + vfac * vel_z)/sfac    #applying RSD along LoS (Z direction) and scaling
    del vel_z
    DM_data[:,0:3]/=sfac                                   #scaling all the three real-space co-ordinates
    gc.collect()                                           #collect garbages

    #applying periodic boundary condition to the RS space co-ordinates so that any particle is not out of the box
    DM_data[:,3] += N3
    DM_data[:,3]=DM_data[:,3] - 1.0 * N3 * np.trunc(np.floor(DM_data[:,3])/N3)  
  
    DM_data[:,4]=1                                         #taking unit mass of each DM particle
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
    halo_data[:,1:4] /= sfac                              #scaling the co-ordinates
    halo_data=np.round(halo_data,6)
    #******************************************************#    


    #cic on DM_data and halo_data to get nh and ngamma respectively
    #******************************************************#
    dimensions=[N1,N2,N3]

    ngamma = cic_vmass(halo_data,dimensions,xindex=1,yindex=2,zindex=3,mindex=0) #density of ionizing sources
    del halo_data                                          
    gc.collect()

    start_cic=time.time()
    nh= cic_vmass(DM_data,dimensions,xindex=0,yindex=1,zindex=2,mindex=4,run_parallel=cic_parallel_flag,num_threads=cic_threads)
                                                                                          #density of hydrogen
    print("CIC done in {:.2f} sec".format(time.time()-start_cic))
    #******************************************************#


    #getting ionisation field using excursion set formalism
    #******************************************************#
    nxion=get_ion(nh=nh,ngamma=ngamma,grid_spacing=LL,nion=nion,rmfp=rmfp)           #nxion is ionization fraction at each grid point
    #******************************************************#


    #getting neutral hydrogen density field in real space
    #******************************************************#
    np.maximum(0, 1 - nxion, out=nxion)                    #neutral hydrogen frac = 1-nxion; 0 if it comes out to be negative; storing it in nxion
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
    DM_data[:,4]=density_to_mass(nxion, DM_data, xindex=0, yindex=1, zindex=2)  #getting HI mass at particle pos; storing the mass back to data[:,4]
    nh=cic_vmass(DM_data,dimensions,xindex=0,yindex=1,zindex=3,mindex=4,run_parallel=cic_parallel_flag,num_threads=cic_threads)
                                                                                #getting HI density at grids taking particles in RS space

    #storing redshift-space data into file
    filename4="ionz_out/HI_maprs_{:.3f}".format(z)
    with open(filename4, 'wb') as file:
        file.write(np.array([N1, N2, N3], dtype=np.int32).tobytes())
        file.write(nh.tobytes())

    
    print("Mass averaged Neutral Hydrogen Fraction (for data in redshift space) = {:.6f}".format(np.average(nh)/robar))
    #******************************************************#
    print("-------------------------------------------------")
    gc.collect()
print("*****************************************************")

#calculating total time for execution of the code
#***********************************************************#
tot_time=time.time()-start
hr=int(tot_time/3600)
minute=int((tot_time-int(tot_time/3600)*3600)/60)
sec=tot_time-int(tot_time/60)*60
print("Total time taken={} hr {} min {} sec".format(hr,minute,int(sec)))
#****************************END****************************#

