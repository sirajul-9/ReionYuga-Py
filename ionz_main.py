import numpy as np
from funcs import *
import time
import gc
import sys


nbody_path="/home/sirajul/packages/N-body"                 #write the n-body output file path here, don't put slash at the end
halo_cat_path='/home/sirajul/packages/FoF-Halo-finder'              #write the halo_catalogue file path here, don't put slash at the end
sfac = 2                                                   #factor for rescaling the  grid
Nthreads=2                                                 #set the number of openmp threads

nion=23.21
rmfp=20.00

if len(sys.argv)<2:
    print("Please enter redshift values")
    sys.exit(0)

redshifts=np.array(sys.argv[1:]).astype(np.float32)
print("\nGenerating HI maps at {} redshifts".format(len(redshifts)))
print("----------------------------------------------")
start=time.time()
set_num_threads(Nthreads)
if not os.path.exists("ionz_out"):
        os.makedirs("ionz_out")

for z in redshifts:
    print("Redshift = {:.3f}".format(z))

    #reading nbody_outputs
    print("Reading N-body output file")
    filename1 = nbody_path+"/output.nbody_{:.3f}".format(z)      
    filename1=filename1.encode('utf-8')
    box = get_box_info(filename1)
    N1,N2,N3=box.grid.dim
    LL=box.grid.grid_spacing
    tot_DM=box.tot_DM
    vaa=box.scale_factor
    pos,vel = read_nbody_output(filename1)
    print("N1 = {} N2 = {} N3 = {} LL={:.4f} and total DM = {}".format(N1,N2,N3,LL,tot_DM))
   
    #redefining grid
    N1=int(N1/sfac)  
    N2=int(N2/sfac) 
    N3=int(N3/sfac)  
    LL=LL*sfac; 
    robar=tot_DM/(1.*N1*N2*N3) #mean number density (grid)^{-3}
    vfac =1./(hubble_param(vaa)*vaa*vaa) #for redshift space distortion

    #storing scaled position and velocity of dark matter particles in data array
    data = np.zeros((tot_DM,5),dtype=np.float32)
    data[:,0:3]=pos/sfac

    #applying redshift space distortion along Z(LoS)
    data[:,3] = (pos[:,2]+vfac*vel[:,2])/sfac

    #applying periodic boundary condition
    data[:,3] += N3
    data[:,3]=data[:,3] - 1.0 * N3 * np.trunc(np.floor(data[:,3])/N3)    
    data[:,4]=1 
    data=np.round(data,6)

    #releasing memory which is not needed anymore
    del pos, vel
    gc.collect()

    #reading halo_catalogues
    print("Reading halo catalogues")
    filename2 =halo_cat_path+"/halo_catalogue_{:.3f}".format(z)      
    filename2=filename2.encode('utf-8')
    tot_cluster,halo= read_halo_catalogue(filename2)
    print("Total number of halos = {}".format(tot_cluster))
    #redefining halo
    halo[:,1:4] /= sfac
    halo=np.round(halo,6)

    
    #cic on data and halo
    dimensions=[N1,N2,N3]
    nh= cic_vmass(data,dimensions,0,1,2,4)
    ngamma = cic_vmass(halo,dimensions,1,2,3,0)
    print("CIC done")

    #get ionisation field using excursion set formalism
    nxion=get_ion(nh,ngamma,LL,nion,rmfp)

    #getting neutral hydrogen density field
    np.maximum(0, 1 - nxion, out=nxion)
    nh *= nxion
    print("volume averaged neutral hydrgen fraction(for real-space data) = {:.6f}".format(np.average(nxion)))
    print("Mass averaged Neutral Hydrogen Fraction (for data in real space)= {:.6f}".format(np.average(nh)/robar))
    
    #storing real-space data into file
    filename3="ionz_out/HI_map_{:.3f}".format(z)
    with open(filename3, 'wb') as file:
        file.write(np.array([N1, N2, N3], dtype=np.int32).tobytes())
        file.write(nh.tobytes())
 
    
    #getting neutral hydrogen density field in redshift space
    print("mapping to redshift space")
    density_2_mass(nxion, data,dimensions, 0, 1, 2, 4)
    nh=cic_vmass(data,dimensions, 0, 1, 3, 4)

    #storing redshift-space data into file
    filename4="ionz_out/HI_maprs_{:.3f}".format(z)
    with open(filename4, 'wb') as file:
        file.write(np.array([N1, N2, N3], dtype=np.int32).tobytes())
        file.write(nh.tobytes())

    
    print("Mass averaged Neutral Hydrogen Fraction (for data in redshift space) = {:.6f}".format(np.average(nh)/robar))
    print("-------------------------------------------------")
tot_time=time.time()-start
hr=int(tot_time/3600)
minute=int((tot_time-int(tot_time/3600)*3600)/60)
sec=tot_time-int(tot_time/60)*60
print("Total time taken={} hr {} min {} sec".format(hr,minute,int(sec)))

