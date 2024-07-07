#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fftw3.h>
#include<omp.h>
#include"source.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static float ***rosp; // sphere for smoothing
static fftwf_plan p_rosp; // for FFT
static float ***nhs, ***ngammas;
int nthreads;
typedef struct 
  {
    long      npart[6];
    double   mass[6];
    double   time;
    double   redshift;
    int     flag_sfr;
    int     flag_feedback;
    long    npartTotal[6];
    int     flag_cooling;
    int     num_files;
    double   BoxSize;
    double   Omega0;
    double   OmegaLambda;
    double   HubbleParam; 
    double   Omegab; 
    double   sigma_8_present;
    long  Nx;
    long  Ny;
    long  Nz;
    float LL;
    int output_flag;
    int in_flag;
    long int seed;
    char  fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 6*8 -6*4 -8];  // fills to 256 Bytes 
  }  io_header;

float  vhh, // Hubble parameter 
  vomegam, // Omega_matter; total matter density (baryons+CDM) parameter
  vomegalam, // Cosmological Constant 
  vomegab, //Omega_baryon
  tcmb, //  CMB temperature
  sigma_8_present ,//  Last updated value of sigma_8 (Presently WMAP)
  vnn; // Spectral index of primordial Power spectrum
long N1,N2,N3;
io_header header1;
long    MM; // Number of particles 
float   LL; // grid spacing in Mpc
float   DM_m,// Darm matter mass of simulation particle (calculated) 
  norm, // normalize Pk
  pi=M_PI;

void parallelize(int Nthreads){
  omp_set_num_threads(Nthreads);
  int e=fftwf_init_threads(); 
  if (e==0) printf("Error in FFT multi-threading\n");
  fftwf_plan_with_nthreads(Nthreads);
  printf("Number of openmp threads = %d\n",Nthreads);
  nthreads=Nthreads;
}


void read_output(char *fname, int read_flag,long *box,float *rra,float *vva,float *aa)
{
  FILE *fp1;
  long ii;
  float vaa;
  fp1=fopen(fname,"r");
  int output_flag,in_flag;
  long int seed;
  int kk,dummy;
  
  // header reading
  fread(&dummy,sizeof(dummy),1,fp1);
  fread(&header1,sizeof(io_header),1,fp1);
  fread(&dummy,sizeof(dummy),1,fp1);
  
  
  vaa=(float)header1.time;     //scale factor of  nbody output
  aa[6]=vaa;
  
  MM=(long)header1.npartTotal[1]; //Total DM particles in this simulation
  DM_m=(float)header1.mass[1]; //DM particle mass in units of 10^10 M_sun/h
  
  vomegam=(float)header1.Omega0;    //Omega_
  vomegalam=(float)header1.OmegaLambda;     //OmegaLambda
  vhh=(float)header1.HubbleParam;     //HubbleParam
  vomegab=(float)header1.Omegab;    //Omega_b
  sigma_8_present=(float)header1.sigma_8_present;
  N1=(long)header1.Nx;
  N2=(long)header1.Ny;
  N3=(long)header1.Nz;
  LL=(float)header1.LL;
  output_flag=header1.output_flag; // input units ? (!=1) => kp/h; else grid
  in_flag=header1.in_flag; //  input file generated by ? 1 => zel  else nbody
  seed=header1.seed;

  aa[0]=LL;
  box[0]=N1;
  box[1]=N2;
  box[2]=N3;
  box[3]=MM;
  aa[1]=vomegam;
  aa[2]=vomegalam;
  aa[3]=vomegab;
  aa[4]=vhh;
  aa[5]=DM_m;
  if(read_flag!=1)
    {
      fread(&dummy,sizeof(dummy),1,fp1);
      for(ii=0;ii<MM;ii++)
	{
	  fread(rra+ii*3,sizeof(float),3,fp1);
	  
	  if(output_flag!=1)
	    {
	      *(rra+ii*3)=*(rra+ii*3)/(LL*1000.*vhh);//coordinates in grid
	      *(rra+ii*3+1)=*(rra+ii*3+1)/(LL*1000.*vhh);
	      *(rra+ii*3+2)=*(rra+ii*3+2)/(LL*1000.*vhh);
	      
	      rra[ii*3] = rra[ii*3]-1.0*N1*(long)(floor(rra[ii*3])/(1.*N1));
	      rra[ii*3+1] = rra[ii*3+1]-1.0*N2*(long)(floor(rra[ii*3+1])/(1.*N2));  // imposing periodic boundary condition
	      rra[ii*3+2] = rra[ii*3+2]-1.0*N3*(long)(floor(rra[ii*3+2])/(1.*N3));
	    }
	}
      fread(&dummy,sizeof(dummy),1,fp1);
      
      fread(&dummy,sizeof(dummy),1,fp1);
      for(ii=0;ii<MM;ii++)
	{
	  fread(&vva[ii*3],sizeof(float),3,fp1);
	  
	  if (output_flag!=1)
	    {
	      vva[ii*3]=vva[ii*3]/(LL*vhh*100./vaa);//velocities in 
	      vva[ii*3+1]=vva[ii*3+1]/(LL*vhh*100./vaa);
	      vva[ii*3+2]=vva[ii*3+2]/(LL*vhh*100./vaa);
	    }
	}
      fread(&dummy,sizeof(dummy),1,fp1);
    }
  
  fclose(fp1);
}


void read_fof(char *fname, int read_flag, long *totcluster, float *halo, float *aa)
{
  FILE *fp1;
  long ii,t;
  float vaa;
  fp1=fopen(fname,"r");
  
  int kk,dummy;
  int output_flag;
  
  //***********************************************************************************
  //                               header reading
  //***********************************************************************************
  
  fread(&dummy,sizeof(dummy),1,fp1);
  fread(&header1,sizeof(io_header),1,fp1);
  fread(&dummy,sizeof(dummy),1,fp1);
  

  
  vaa=(float)header1.time;     //scale factor of  nbody output
  *aa=vaa;
  
  MM=(long)header1.npartTotal[1]; //Total DM particles in this simulation
  DM_m=(float)header1.mass[1]; //DM particle mass in units of 10^10 M_sun/h
  
  vomegam=(float)header1.Omega0;    //Omega_
  vomegalam=(float)header1.OmegaLambda;     //OmegaLambda
  vhh=(float)header1.HubbleParam;     //HubbleParam
  vomegab=(float)header1.Omegab;    //Omega_b
  sigma_8_present=(float)header1.sigma_8_present;
  N1=(long)header1.Nx;
  N2=(long)header1.Ny;
  N3=(long)header1.Nz;
  LL=(float)header1.LL;
  output_flag=header1.output_flag; // input units ? (!=1) => kp/h; else grid
  //*in_flag=header1.in_flag; //  input file generated by ? 1 => zel  else nbody
  
  //***********************************************************************************
  //                        reading total number of haloes
  //***********************************************************************************
  fread(&dummy,sizeof(dummy),1,fp1);
  fread(&t,sizeof(long),1,fp1);
  fread(&dummy,sizeof(dummy),1,fp1);
  *totcluster=t;
  
  //***********************************************************************************
  //                        reading halo position and velocity
  //***********************************************************************************
  
  if(read_flag!=1)
    {
      fread(&dummy,sizeof(dummy),1,fp1);
      for(ii=0; ii<t; ii++)
	  fread(&halo[ii*7],sizeof(float),7,fp1); 
      fread(&dummy,sizeof(dummy),1,fp1);
      
      if(output_flag!=1)
	  for(ii=0; ii<t; ii++)
	  {
	    halo[ii*7]=halo[ii*7]/DM_m;

	    halo[ii*7+1]=halo[ii*7+1]/(LL*1000.*vhh);   //coordinates in grid
	    halo[ii*7+2]=halo[ii*7+2]/(LL*1000.*vhh);
	    halo[ii*7+3]=halo[ii*7+3]/(LL*1000.*vhh);
	    halo[ii*7+1] = halo[ii*7+1]-1.0*N1*(long)(floor(halo[ii*7+1])/(1.*N1));
	    halo[ii*7+2] = halo[ii*7+2]-1.0*N2*(long)(floor(halo[ii*7+2])/(1.*N2));   // imposing periodic boundary condition
	    halo[ii*7+3] = halo[ii*7+3]-1.0*N3*(long)(floor(halo[ii*7+3])/(1.*N3));
	    
	    halo[ii*7+4]=halo[ii*7+4]/(LL*vhh*100./vaa);  //velocities
	    halo[ii*7+5]=halo[ii*7+5]/(LL*vhh*100./vaa);
	    halo[ii*7+6]=halo[ii*7+6]/(LL*vhh*100./vaa);
	  } 
    }
  
  fclose(fp1);
}


//function for CIC
void cic_vmass(float *ro_dum,float *data,long tot_particles, long N1,long N2,long N3,int xin, int yin, int zin, int min, int col)
/* This uses Cloud in Cell for calculating density given posns.*/
/* The i/p is simply the array containing the posns. */
{
  int ii, jj, kk, ix, jy, kz;
  long i, j, k, a[2], b[2], c[2], pin;
  float xx,yy,zz,delx,dely,delz,wx,wy,wz,W;

  
  /* Clear out the array ro. ******/
  for (int iii=0;iii<N1;iii++){
    for (int jjj=0;jjj<N2;jjj++)
      for (int kkk=0;kkk<N3;kkk++)
           ro_dum[iii*N2*N3 + jjj*N3 + kkk] = 0.0;
    }

   for(pin=0;pin<tot_particles;pin++)
    { // begin particle index loop 
      // (a/b/c)[0] or (a/b/c)[1] can never be greater than (N1/N2/N3) 
      
      a[0]=floor(data[pin*col+xin]);
      b[0]=floor(data[pin*col+yin] );
      c[0]=floor(data[pin*col+zin]);
      
      a[1]=(a[0]+1);
      b[1]=(b[0]+1);
      c[1]=(c[0]+1);
      
      xx=data[pin*col+xin];
      yy=data[pin*col+yin];
      zz=data[pin*col+zin];
      
      // for each of the 8 corner points 
      for(ii=0;ii<=1;ii++)
	    for(jj=0;jj<=1;jj++)
	     for(kk=0;kk<=1;kk++)
	     { // begin 8 corners loop 
	      delx = xx - a[ii];
	      dely = yy - b[jj];
	      delz = zz - c[kk];
	      
	      ix=a[ii]%N1;
	      jy=b[jj]%N2;
	      kz=c[kk]%N3;
	      
	      // assigning of weights to the corners 
	      wx=1.0-fabs(delx);
	      wy=1.0-fabs(dely);
	      wz=1.0-fabs(delz);
	      W=wx*wy*wz*data[pin*col+min]; //multiplying the product of weights with mass of halo
	      
	      ro_dum[ix*N2*N3+jy*N3+kz]+= W;
	    } // end of <8 grid corners loop>	
    } //end of each particle loop      
        
} // end of function cic_vmass






void smooth(float ***ro_dum, float Radii, long N1, long N2, long N3)
{
  long index, i, j, k; 
  float tempre,tempim;
  double tot;
  
  fftwf_complex *A;
  fftwf_complex *B;
  
  //generating the filtering function
  #pragma omp parallel for private(i,j,k) default(shared)
  for(i=0;i<N1;i++)
    for(j=0;j<N2;j++)
      for(k=0;k<N3;k++)
  	rosp[i][j][k]=0.0;
  
  //generating a sphere at the centre of the box
  
  tot=0.;
  #pragma omp parallel for private(i,j,k) reduction(+:tot)
  for(i=0;i<N1;i++)
    for(j=0;j<N2;j++)
      for(k=0;k<N3;k++)
  	{
	  if((float)((N1/2-i)*(N1/2-i)+(N2/2-j)*(N2/2-j)+(N3/2-k)*(N3/2-k))<=Radii*Radii)
	    rosp[i][j][k]=1.0;//centre N1/2,N2/2,N3/2
	  
	  tot += (double)rosp[i][j][k];
	}
  
  
  //Sphere generation complete 
  //Doing Fourier Transform of the sphere
  
  fftwf_execute(p_rosp);
  B=(fftwf_complex*)&(rosp[0][0][0]);
  
  //We will multiply the factor powf((-1.),(i+j+k)) with FT of the sphere to shift it to one corner of the box from box centre while applying boundary condition below
  //----------------------------------------------------------------------
  
  //Doing Fourier Transform of the density field
  
  
  fftwf_plan p_ro_dum=fftwf_plan_dft_r2c_3d(N1, N2, N3, &(ro_dum[0][0][0]), (fftwf_complex*)&(ro_dum[0][0][0]),FFTW_ESTIMATE);
  
  fftwf_plan q_ro_dum=fftwf_plan_dft_c2r_3d(N1, N2, N3, (fftwf_complex*)&(ro_dum[0][0][0]), &(ro_dum[0][0][0]), FFTW_ESTIMATE);
  
  
  fftwf_execute(p_ro_dum);
  A=(fftwf_complex*)&(ro_dum[0][0][0]);
  
  #pragma omp parallel for private(i,j,k,tempre,tempim,index)
  for(i=0;i<N1;i++)
    for(j=0;j<N2;j++)
      for(k=0;k<(N3/2+1);k++)
	{ 
	  
	  index = i*N2*(N3/2 +1) + j*(N3/2 +1) + k;
	  
	  tempre=(A[index][0]*B[index][0]-A[index][1]*B[index][1])*powf((-1.),1.*(i+j+k))/(1.*tot);
	  tempim=(A[index][1]*B[index][0]+A[index][0]*B[index][1])*powf((-1.),1.*(i+j+k))/(1.*tot);
	  //multiplying the factor powf((-1.),(i+j+k)) with FT of the sphere to shift it to one corner of the box from box centre	  
	  A[index][0]=tempre;
	  A[index][1]=tempim;
	}
  
  fftwf_execute(q_ro_dum);
  #pragma omp parallel for private(i,j,k)
  for(i=0;i<N1;i++)
    for(j=0;j<N2;j++)
      for(k=0;k<=N3;k++)
  	    ro_dum[i][j][k] = ro_dum[i][j][k]/(1.*N1*N2*N3);
  
  fftwf_destroy_plan(p_ro_dum);
  fftwf_destroy_plan(q_ro_dum); 
}

float  ***allocate_fftwf_3d(long N1,long N2,long N3)
{
  long ii,jj;
  long asize,index;
  float ***phia, *phi;

  phia=(float ***) fftwf_malloc (N1 *  sizeof(float **));


  for(ii=0;ii<N1;++ii)
      phia[ii]=(float **) fftwf_malloc (N2 *  sizeof(float *));

  asize = N1*N2;
  asize = asize*N3;

  if(!(phi = (float *) calloc(asize,sizeof(float))))
    {
      printf("error in allocate_fftwf_3d");
      exit(0);
    }

  for(ii=0;ii<N1;++ii)
    for(jj=0;jj<N2;++jj)
      {
	index = N2*N3;
	index = index*ii + N3*jj;
	phia[ii][jj]=phi+ index;
      }
  return(phia);
}


void deallocate_3d_double(float ***array, long N1, long N2) {
    if (array == NULL) {
        return; // Check for NULL pointer to avoid issues with double freeing
    }
     if (array[0][0] != NULL) {
            
            free(array[0][0]); // Free the pointer to 1D array phi
        }
    for (long ii = 0; ii < N1; ++ii) {
        if (array[ii] != NULL) {
            
            fftwf_free(array[ii]); // Free the 1D array of pointers 
        }
    }
    fftwf_free(array); // Free the top-level array of pointers
}

void get_nhs(float *nh, float *ng, float *nx,long N1, long N2, long N3,float LL, float nion, float rmfp){
    
    long ii,jj,kk;
    rosp = allocate_fftwf_3d(N1,N2,N3+2);
    nhs=allocate_fftwf_3d(N1,N2,N3+2);
    ngammas=allocate_fftwf_3d(N1,N2,N3+2);
    p_rosp = fftwf_plan_dft_r2c_3d(N1, N2, N3, &(rosp[0][0][0]), (fftwf_complex*)&(rosp[0][0][0]), FFTW_ESTIMATE);
    #pragma omp parallel for private(ii,jj,kk)
    for(ii=0;ii<N1;ii++)
	  for(jj=0;jj<N2;jj++)
	  for(kk=0;kk<N3;kk++)
	    {
          if(nh[ii*N2*N3+jj*N3+kk]>nion*ng[ii*N2*N3+jj*N3+kk]) // checking ionization condition
		{
		  nx[ii*N2*N3+jj*N3+kk]=nion*ng[ii*N2*N3+jj*N3+kk]/nh[ii*N2*N3+jj*N3+kk];
		}
	      
	      else
		{
		  nx[ii*N2*N3+jj*N3+kk]=1.;
		}
	    }
     float r_max=rmfp/LL; // Mpc/LL in grid unit

     float Radii=1;
     float dr;
      while(Radii < r_max)
	{
      #pragma omp parallel for private(ii,jj,kk)
	  for(ii=0;ii<N1;ii++)
	    for(jj=0;jj<N2;jj++)
	      for(kk=0;kk<N3;kk++)
		{
		  nhs[ii][jj][kk]=nh[ii*N2*N3+jj*N3+kk];
		  ngammas[ii][jj][kk]=ng[ii*N2*N3+jj*N3+kk];
		}
	  //printf("starting smoothing for radius of size %e\n",Radii);
	  
	  smooth(nhs,Radii,N1,N2,N3);
	  
	  smooth(ngammas,Radii,N1,N2,N3);
	  
	  #pragma omp parallel for private(ii,jj,kk)
	  for(ii=0;ii<N1;ii++)
	    for(jj=0;jj<N2;jj++)
	      for(kk=0;kk<N3;kk++)
		{
		  if(nhs[ii][jj][kk]<=nion*ngammas[ii][jj][kk])  // checking ionization condition
		    nx[ii*N2*N3+jj*N3+kk]=1.;
		}
	  
	  dr=(Radii*0.1) < 2.0 ? (Radii*0.1) : 2.0; //increment of the smoothing radius
	  Radii += dr;
	}
    deallocate_3d_double(rosp,N1,N2);
    deallocate_3d_double(nhs,N1,N2);
    deallocate_3d_double(ngammas,N1,N2); 
    fftwf_destroy_plan(p_rosp);
}


void density_2_mass(float *ro_dum,float *data,long MM, long N1,long N2,long N3,int xin, int yin, int zin, int col, float *mass)
{
  long pin;
  #pragma omp parallel for private(pin)
  for(pin=0;pin<MM;pin++)
    mass[pin]=0;


  #pragma omp parallel for private(pin)
  for(pin=0;pin<MM;pin++) /* begin particle index loop */
    {
      //rr=data[pin];
      //rr[min]=0.0;
      long ii,jj,kk;
      float delx,dely,delz, wx,wy,wz,W,*rr;
      long a[2],b[2],c[2],ix,jy,kz;
      float xx,yy,zz;
      /* left most corner of the cube enclosing the particle */
      
      a[0]=(long)floor(data[pin*col+xin]);
      b[0]=(long)floor(data[pin*col+yin]);
      c[0]=(long)floor(data[pin*col+zin]);                  
      
      /* right most corner of the cube enclosing the particle */
      
      a[1]= a[0]+1;
      b[1]= b[0]+1;
      c[1]= c[0]+1;
      
      /* particle co-ordinates itself */
      xx=data[pin*col+xin];
      yy=data[pin*col+yin];
      zz=data[pin*col+zin];
      
      
      for(ii=0;ii<=1;ii++)
        for(jj=0;jj<=1;jj++)
          for(kk=0;kk<=1;kk++)
            { /* begin 8 corners loop */
	          ix = a[ii]%N1;
              jy = b[jj]%N2;
              kz = c[kk]%N3;
	      /* ix,jy,kz are the current co-ordinates of the cube vertex point */
	      
              /* calculating the difference from the respective corner */   
              delx = xx - a[ii];
              dely = yy - b[jj];
              delz = zz - c[kk];
	      
              /* assigning of weights to the points acc to distance from pts*/
              wx=1.0-fabs(delx);
              wy=1.0-fabs(dely);
              wz=1.0-fabs(delz);
	      
              W = wx*wy*wz;
	      
	      //data[pin*col+min]+=W*ro_dum[ix*N2*N3+jy*N3+kz];
           mass[pin]+=W*ro_dum[ix*N2*N3+jy*N3+kz];
	    } /* End of 8 corners loop */
    } /* end particle index loop */
} /* end density_2_mass */


//Function for calculating power spectrum
void calpow_mom(float *ro,int Nbin,double* power, double* kmode, double* power_P2,double* power_P4, double *no, long N1, long N2, long N3, float LL)
{ 
  
  /******************** TO FIND POWER SPECTRUM **************************/
  
  long i, j, k, a, b, c, d;
  long index, index1, index2;
  
  double *no2, *no4;
  fftwf_complex *comp_ro;
  
  float fac1, fac2, fac3, m, mu, P2, P4, scale;
  double norml;
  
  norml=1./(1.*N1*N2*N3);
  
  float ***ro_dum= allocate_fftwf_3d(N1,N2,N3+2);
  
  #pragma omp parallel for private(j,k)
  for(i=0;i<N1;i++)
    for(j=0;j<N2;j++)
      for(k=0;k<N3;k++)
        ro_dum[i][j][k]=ro[i*N2*N3 + j*N3 + k];
  
  /*************** TAKING FOURIER TRANSFORM OF RO. **************/
  
  fftwf_plan p_ro_dum=fftwf_plan_dft_r2c_3d(N1, N2, N3, &(ro_dum[0][0][0]), (fftwf_complex*)&(ro_dum[0][0][0]),FFTW_ESTIMATE);
  
  fftwf_plan q_ro_dum=fftwf_plan_dft_c2r_3d(N1, N2, N3, (fftwf_complex*)&(ro_dum[0][0][0]), &(ro_dum[0][0][0]),FFTW_ESTIMATE);
  
  
  fftwf_execute(p_ro_dum);
  
  comp_ro=(fftwf_complex*)&(ro_dum[0][0][0]);
  
  /*********** TO FIND POWER SPECTRUM OF RO. **************/
  
  no2=calloc((size_t)Nbin,sizeof(double));
  no4=calloc((size_t)Nbin,sizeof(double));
  
  fac1=1./(1.*N1*N1);
  fac2=1./(1.*N2*N2);
  fac3=1./(1.*N3*N3);
  
  /**************** BINNING POWER SPECTRA **********************/
  
  for(i=0;i<Nbin;++i)
    {
      kmode[i]=0.0;
      
      power[i]=0.0;
      power_P2[i]=0.0;
      power_P4[i]=0.0;
      
      no[i]=0.0;
      no2[i]=0.0;
      no4[i]=0.0;
    }
  
  //***********************************************************************
  
  scale=log10(0.5*N1)/Nbin;
  
  //***********************************************************************

  /*---------------------- BINNING POWER SPECTRA -------------------*/
  
  /*-------------------------- half lines ------------------------- */
  
  for(i=1;i<=N1/2;i++)
    for(j=0;j<=N2/2;j=j+N2/2)
      for(k=0;k<=N3/2;k=k+N3/2)
	{
	  index = i*N2*(N3/2+1) + j*(N3/2+1) + k;
	  
	  m = sqrt(fac1*i*i + fac2*j*j + fac3*k*k);	      
	  
	  mu = (1.0*k)/(N3*m);
	  
	  P2 = 0.5*(3.*mu*mu - 1.);
	  P4 = 0.125*(35.0*powf(mu,4.0)-30.0*mu*mu + 3);
	  
	  d=(int)floorf(log10(m*N1)/scale);  // logarithmic bins
	  
	  if(d>=0 && d<Nbin)
	    {
	      power[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))); 
	      power_P2[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))*P2);
	      power_P4[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))*P4);
	      
	      kmode[d]+=(double)(m);
	      no[d]= no[d] + (double)(1.0);
	      no2[d]=no2[d] + (double)(P2*P2);
	      no4[d]=no4[d] + (double)(P4*P4);
	    }
	}	  
  
  /*----------------------- half planes -----------------------*/
  
  for(i=0;i<N1;i++)
    {
      a=(i>N1/2)? N1-i: i;
      index1 = i*N2*(N3/2+1) ;
      
      for(j=1;j<N2/2;j++) 
	{
	  b=j; 
	  index2 = index1 + j*(N3/2+1) ;
	  
	  for(k=0;k<=N3/2;k=k+N3/2)
	    {
	      c=k;
	      index = index2 + k;
	      
	      m = sqrt(fac1*a*a + fac2*b*b + fac3*c*c);	      
	      
	      mu = (1.0*c)/(N3*m);
	      
	      P2 = 0.5*(3.*mu*mu - 1.);
	      P4 = 0.125*(35.0*powf(mu,4.0)-30.0*mu*mu + 3);
	      
	      d=(int)floorf(log10(m*N1)/scale);//logarithmic bins
	      
	      if(d>=0 && d<Nbin)
		{
		  power[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))); 
		  power_P2[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))*P2);
		  power_P4[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))*P4);
		  
		  kmode[d]+=(double)(m);
		  no[d]= no[d] + (double)(1.0);
		  no2[d]=no2[d] + (double)(P2*P2);
		  no4[d]=no4[d] + (double)(P4*P4);
		}
	      
	    }	  
	}
    }
  
  /**************** half cube **********************/
  
  for(i=0;i<N1;i++)
    {
      a=(i>N1/2)? N1-i: i;
      index1 = i*N2*(N3/2+1) ;
      
      for(j=0;j<N2;j++)
	{
	  b=(j>N2/2)? N2-j: j;
	  index2 = index1 + j*(N3/2+1) ;
	  
	  for(k=1;k<N3/2;k++)
	    {
	      c=k;	  	      
	      index = index2 + k;
	      
	      m = sqrt(fac1*a*a + fac2*b*b + fac3*c*c);	      
	      
	      /* m*(2 * pi/LL) is |k| */
	      /* m=1/2 corresponds to kmode[Nbin-1] i.e. Nyquits */
	      
	      mu = (1.0*c)/(N3*m);
	      
	      P2 = 0.5*(3.*mu*mu - 1.);
	      P4 = 0.125*(35.0*powf(mu,4.0)-30.0*mu*mu + 3);
	      
	      d=(int)floorf(log10(m*N1)/scale);//logarithmic bins
	      
	      if(d>=0 && d<Nbin)
		{
		  power[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))); 
		  power_P2[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))*P2);
		  power_P4[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))*P4);
		  
		  kmode[d]+=(double)(m);
		  no[d]= no[d] + (double)(1.0);
		  no2[d]=no2[d] + (double)(P2*P2);
		  no4[d]=no4[d] + (double)(P4*P4);
		}
	    } /* end k for loop */
	  
	}
    }
  
  //***********************************************************************
  //***********************************************************************
  
  for(i=0;i<Nbin;i++)
    {
      if (no[i]>0.0) 
	{
	  power[i] = pow(LL,3.)*norml*power[i]/(1.0*no[i]);	  
	  kmode[i]=2.*pi*kmode[i]/(no[i]*LL);
	}
      if(no2[i]>0.0)
	power_P2[i] = pow(LL,3.)*norml*power_P2[i]/(1.0*no2[i]);
      
      if(no4[i]>0.0)
	power_P4[i] = pow(LL,3.)*norml*power_P4[i]/(1.0*no4[i]);
    }
  
  //***********************************************************************
  
  for(i=0;i<N1;i++)
    {
      index1 = i*N2*(N3/2+1) ;	  
      for(j=0;j<N2;j++)
	{
	  index2=index1 + j*(N3/2+1) ;
	  for(k=0;k<(N3/2+1);k++)
	    {
	      index=index2 + k;
	      comp_ro[index][0]=comp_ro[index][0]*norml;
	      comp_ro[index][1]=comp_ro[index][1]*norml;
	    }
	}
    }
  
  /*  now convert the array back to real space */
  
  //***********************************************************************
  
  fftwf_execute(q_ro_dum);
  fftwf_destroy_plan(p_ro_dum);
  fftwf_destroy_plan(q_ro_dum); 
  deallocate_3d_double(ro_dum,N1,N2);
} /* end function */

