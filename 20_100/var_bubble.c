#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<complex.h> 
#include <time.h>
#include <float.h>
#include"funcwritepm3d.c"
#include<fftw3.h> 
#include"funcwriteoned.c"
//#include"ran2.c"



int main()
{

//variable declaration start
FILE *fpr,*fpw;
char ch[50],NAME[50];
int i,j,k,l,i1,i2,i3,ijk;
double li,li2,maxeta;
int iw,js,ie,jn,kn,ks;
double temp,etatemp,setatemp,laplacian_n,laplacian_nn;
double rxtemp;
int xtemp,ytemp,ztemp,xx,yy,zz;
int *n_nx,*feta,*maxidx;
double *ipeta,*ineta,*jpeta,*jneta,*kneta,*kpeta;

double *s_eta,*sig_eta,*M_var;
int *nap;


double min,max,tmpelement,AA;
int ns,a;
int *nsize;
int nx,ny,nz,no,nmax,istep,ostep,nstep;
double dx,dy,dz,dt;
double gama,cutoff;
double realc;
fftw_complex *g,*comp;
fftw_complex *comp_partx,*comp_party,*comp_partz;
fftw_plan planF, planB;


if( (fpr = fopen("system_size","r")) == NULL){
printf("Unable to open system_size. Exiting\n");
exit(0);
}
else{
fpr = fopen("system_size","r");
}
fscanf(fpr,"%d%d%d%le%le%le",&nx,&ny,&nz,&dx,&dy,&dz);
fclose(fpr);


if( (fpr = fopen("system_constants","r")) == NULL){
printf("Unable to open system_constants. Exiting\n");
exit(0);
}
else{
fpr = fopen("system_constants","r");
}
fscanf(fpr,"%le%d%d%le%d%d",&gama,&no,&nmax,&dt,&nstep,&ostep);

fclose(fpr);





double *eta[nmax],*eta_t[nmax];
int *num[nmax];

//variable allocation start


n_nx = (int *) malloc( nx*ny*nz* sizeof(int));
maxidx = (int *) malloc( nx*ny*nz* sizeof(int));


feta = (int *) malloc((size_t) no* sizeof(int));

ipeta= (double *) malloc( no* sizeof(double));
ineta= (double *) malloc( no* sizeof(double));
jpeta= (double *) malloc( no* sizeof(double));
jneta= (double *) malloc( no* sizeof(double));
kpeta= (double *) malloc( no* sizeof(double));
kneta= (double *) malloc( no* sizeof(double));

nsize= (int *) malloc((size_t) no* sizeof(int));
nsize= (int *) malloc((size_t) no* sizeof(int));

s_eta = (double *) malloc((size_t) nx*ny*nz* sizeof(double));
sig_eta = (double *) malloc((size_t) nx*ny*nz* sizeof(double));
M_var = (double *) malloc((size_t) nx*ny*nz* sizeof(double));



for (l = 1; l<nmax; ++l)
  {
eta[l] = (double *) malloc((size_t) nx*ny*nz* sizeof(double));
eta_t[l]= (double *) malloc((size_t) nx*ny*nz* sizeof(double));
num[l]= (int *) malloc((size_t) nx*ny*nz* sizeof(int));
}

nap= (int *) malloc((size_t) nx*ny*nz* sizeof(int));


g = fftw_malloc(nx*ny*nz*sizeof(fftw_complex));
comp = fftw_malloc(nx*ny*nz*sizeof(fftw_complex));
comp_partx = fftw_malloc(nx*ny*nz*sizeof(fftw_complex));
comp_party = fftw_malloc(nx*ny*nz*sizeof(fftw_complex));
comp_partz = fftw_malloc(nx*ny*nz*sizeof(fftw_complex));


double mobility;

int half_nx,half_ny,half_nz;
double kx,delta_kx,kx2,delta_ky,ky,ky2,delta_kz,kz,kz2;
double k2, k4,k2y,k4y;
double inv_denom;
double inv_denom_eta;
double inv_nx,inv_ny, inv_nz;

 
double alpha, beta, gamma, kappa, deltat; 
int un1, un2;
long seed;
double g_bub,g_int,beta_s,M,L;
int temp1,temp2,temp3;
double temp4,temp5;

//variable allocation ends

//system("rm -rf output/dataslice/*.dat");
//system("rm -rf output/data/*.dat");
system("rm -rf output/grain/*.dat");
system("rm -rf output/*.dat");
//initializing constants for every eta starts



for (l = 1; l<nmax+1; ++l)
 {
feta[l] = 0.0;
ipeta[l] = 0.0;
ineta[l] = 0.0;
jpeta[l] = 0.0;
jneta[l] = 0.0;
kpeta[l] = 0.0;
kneta[l] = 0.0;
/*
ipjpeta[l] = 0.0;
ipjneta[l] = 0.0;
injpeta[l] = 0.0;
injneta[l] = 0.0;
*/

}

 for(i=0; i < nx; ++i)
		{
		for(j=0; j < ny; ++j)
			{
		for(k=0; k < nz; ++k)
			{
			ijk=k+i*ny*nz+j*nz;
	  for( l=1;l<nmax;++l)
     		{
	eta[l][ijk] = 0.0;
	eta_t[l][ijk] = 0.0;
	num[l][ijk] = 0;
	
	}
	n_nx[ijk] = 0;
	nap[ijk] = 0;
	comp[ijk] = 0.0;


	//printf("%d",i);
}}}



 gamma= gama;
 alpha = gamma;
 beta  = alpha;
 beta_s = gama;
 kappa = 2.0*gamma;
 M = 0.01;
inv_nx = 1.0/(double)(nx);
inv_ny = 1.0/(double)(ny);
inv_nz = 1.0/(double)(nz);

fftw_plan_with_nthreads(48);


half_nx = (int) nx/2;

delta_kx = (2.0*M_PI)/(nx*dx);

half_ny = (int) ny/2;

delta_ky = (2.0*M_PI)/(ny*dy);

half_nz = (int) nz/2;

delta_kz = (2.0*M_PI)/(nz*dz);

//initializing constants for every eta starts

planF = fftw_plan_dft_3d(nx,ny,nz,g,g,FFTW_FORWARD,FFTW_PATIENT);
planB = fftw_plan_dft_3d(nx,ny,nz,g,g,FFTW_BACKWARD,FFTW_PATIENT);

//initial profile starts

seed =-11160;
//srand((unsigned)time(NULL));
 //rxtemp = 40.0;



	sprintf(ch,"1000grain_id.txt");
 if( (fpr = fopen(ch,"r")) == NULL){
	printf("Unable to open file %s",ch);
	printf("Exiting\n");
	exit(0);
	}
else{
	fpw = fopen(NAME,"r");
	}

 for(i=0; i < nx; ++i)
		{
		for(j=0; j < ny; ++j)
{
fscanf(fpr,"%d\t%d\t%d\t%lf\t%lf\n", &temp1,&temp2,&temp3,&temp4,&temp5);
for(k=16; k < 32; ++k)
				{
	
	ijk= k + j*nz + ny*nz*i;
n_nx[ijk] = 1;
     eta[1][ijk] = 1.0;
     num[1][ijk] = temp3;

}
}
}
fclose(fpr);

 for(i=0; i < nx; ++i)
		{
		for(j=0; j < ny; ++j)
			{
for(k=0; k < nz; ++k)
			{
			
xx = (i+nx)%nx;
               yy = (j+ny)%ny;
	ijk= k + yy*nz + ny*nz*xx;
  if(k<16 || k>31)
		{
			n_nx[ijk] = 1;
     			eta[1][ijk] = 0.0;
     			num[1][ijk] = 0;
		        comp[ijk] = 1.0 +_Complex_I*0.0;
		}

    }
   
   }//j
	
   }//i   

  for(i=0; i < nx; ++i)
		{
		for(j=0; j < ny; ++j)
			{
for(k=0; k < nz; ++k)
			{
			
			ijk=k+i*ny*nz+j*nz;
	  for( l=1;l<n_nx[ijk]+1;++l)
     		{
	eta[l][ijk] = fabs(eta[l][ijk]);
        s_eta[ijk] = s_eta[ijk]+eta[l][ijk]*eta[l][ijk];
        sig_eta[ijk] = sig_eta[ijk]+eta[l][ijk];
	
        if(s_eta[ijk] > 1.0) {s_eta[ijk] = 1.0;}
    		}
	maxidx[ijk] = num[1][ijk];

   }
  } 
}
		sprintf(NAME,"output/grain/3d_grain0.dat");
		writepm3d_gg(comp,maxidx,s_eta,sig_eta,NAME,nx,ny,nz);
		sprintf(NAME,"output/grain/2d_grain0.dat");
		writeonedyprofile(comp,maxidx,s_eta,sig_eta,NAME,nx,ny,nz);


//printing initial profile 


//sprintf(NAME,"output/dataslice/comp0.dat");
//writeonedxprofile(comp,NAME,nx,ny,nz);
//sprintf(NAME,"output/data/comp0.dat");
//writepm3d(comp,NAME,nx,ny,nz);


//time loop starts
for(istep=1; istep < nstep; ++istep)
{
	printf("%d\n",istep);

if (istep < 2000)

{ L = 2.0;
dt = 0.05;}
else { L = 2.0;
dt = 0.05;}
for(i=0; i < nx; ++i)
		{
		for(j=0; j < ny; ++j)
			{
for(k=0; k < nz; ++k)
			{
			
			ijk=k+i*ny*nz+j*nz;
  s_eta[ijk] = 0.0;
  sig_eta[ijk] = 0.0;

}}}


  for(i=0; i < nx; ++i)
		{
		for(j=0; j < ny; ++j)
			{
for(k=0; k < nz; ++k)
			{
			
			ijk=k+i*ny*nz+j*nz;
	  for( l=1;l<n_nx[ijk]+1;++l)
     		{
	eta[l][ijk] = fabs(eta[l][ijk]);
        s_eta[ijk] = s_eta[ijk]+eta[l][ijk]*eta[l][ijk];
        sig_eta[ijk] = sig_eta[ijk]+eta[l][ijk];
	
        if(s_eta[ijk] > 1.0) {s_eta[ijk] = 1.0;}
    		}

   }
  } 
}


for(k=0; k < nz; ++k)
{
   kn = k+1;
   if(k==nz-1){kn=0;}
   ks = k-1;
   if(k==0){ks =nz-1;}
  for(j=0; j < ny; ++j)
{
   jn = j+1;
   if(j==ny-1){jn=0;}
   js = j-1;
   if(j==0){js =ny-1;}
   
   for(i=0; i < nx; ++i)
{
    ie = i+1;
    if(i==nx-1){ie=0;}
    iw = i-1;
    if(i==0){iw =nx-1;}
	ijk=k+i*ny*nz+j*nz;
    for( l=1;l<n_nx[ijk]+1;++l)
    {
     feta[num[l][ijk]] = 1;
    }
  
    nap[ijk] = n_nx[ijk];
/*    
    for (l=1;l<n_nx[k+iw*ny*nz+js*nz]+1;++l)
	{
     if (feta[num[l][k+iw*ny*nz+j*nz]] == 0 && feta[num[l][k+i*ny*nz+js*nz]] == 0 && feta[num[l][k+iw*ny*nz+js*nz]]== 0) 
	{
      feta[num[l][k+iw*ny*nz+js*nz]] = 1;
      nap[ijk] = nap[ijk] + 1;

      num[nap[ijk]][ijk] = num[l][k+iw*ny*nz+js*nz];

     }
    }    
    

for (l=1;l<n_nx[k+iw*ny*nz+jn*nz]+1;++l)
   {  if (feta[num[l][k+iw*ny*nz+j*nz]] == 0 && feta[num[l][k+i*ny*nz+jn*nz]] == 0 && feta[num[l][k+iw*ny*nz+jn*nz]]== 0) {
      feta[num[l][k+iw*ny*nz+jn*nz]] = 1;
      nap[ijk] = nap[ijk] + 1;

      num[nap[ijk]][ijk] = num[l][k+iw*ny*nz+jn*nz];

     }
    } 


   for (l=1;l<n_nx[k+ie*ny*nz+js*nz]+1;++l)
    { if (feta[num[l][k+ie*ny*nz+j*nz]] == 0 && feta[num[l][k+i*ny*nz+js*nz]] == 0 && feta[num[l][k+ie*ny*nz+js*nz]]== 0) {
      feta[num[l][k+ie*ny*nz+js*nz]] = 1;
      nap[ijk] = nap[ijk] + 1;

      num[nap[ijk]][ijk] = num[l][k+ie*ny*nz+js*nz];

     }
    }   

 for (l=1;l<n_nx[k+ie*ny*nz+jn*nz]+1;++l)
    { if (feta[num[l][k+ie*ny*nz+j*nz]] == 0 && feta[num[l][k+i*ny*nz+jn*nz]] == 0 && feta[num[l][k+ie*ny*nz+jn*nz]]== 0) {
      feta[num[l][k+ie*ny*nz+jn*nz]] = 1;
      nap[ijk] = nap[ijk] + 1;

      num[nap[ijk]][ijk] = num[l][k+ie*ny*nz+jn*nz];

     }
    }


*/  
    for (l=1;l<n_nx[k+iw*ny*nz+j*nz]+1;++l)
    { if (feta[num[l][k+iw*ny*nz+j*nz]] == 0 ) {
      feta[num[l][k+iw*ny*nz+j*nz]] = 1;
      nap[ijk] = nap[ijk] + 1;

      num[nap[ijk]][ijk] = num[l][k+iw*ny*nz+j*nz];
      eta_t[nap[ijk]][ijk] = 0.0;
     }
    } 


 for (l=1;l<n_nx[k+ie*ny*nz+j*nz]+1;++l)
    { if (feta[num[l][k+ie*ny*nz+j*nz]] == 0 ) {
      feta[num[l][k+ie*ny*nz+j*nz]] = 1;
      nap[ijk] = nap[ijk] + 1;

      num[nap[ijk]][ijk] = num[l][k+ie*ny*nz+j*nz];
      eta_t[nap[ijk]][ijk] = 0.0;
     }
    }
    
    
   for (l=1;l<n_nx[k+i*ny*nz+js*nz]+1;++l)
    { if (feta[num[l][k+i*ny*nz+js*nz]] == 0 ) {
      feta[num[l][k+i*ny*nz+js*nz]] = 1;
      nap[ijk] = nap[k+i*ny*nz+j*nz] + 1;

      num[nap[ijk]][ijk] = num[l][k+i*ny*nz+js*nz];
      eta_t[nap[ijk]][ijk] = 0.0;
     }
    }   
   
        

   for (l=1;l<n_nx[k+i*ny*nz+jn*nz]+1;++l)
    { if (feta[num[l][k+i*ny*nz+jn*nz]] == 0 ) {
      feta[num[l][k+i*ny*nz+jn*nz]] = 1;
      nap[ijk] = nap[k+i*ny*nz+j*nz] + 1;
      num[nap[ijk]][ijk] = num[l][k+i*ny*nz+jn*nz];
      eta_t[nap[ijk]][ijk] = 0.0;
     }
    }    

for (l=1;l<n_nx[ks+i*ny*nz+j*nz]+1;++l)
    { if (feta[num[l][ks+i*ny*nz+j*nz]] == 0 ) {
      feta[num[l][ks+i*ny*nz+j*nz]] = 1;
      nap[ijk] = nap[k+i*ny*nz+j*nz] + 1;

      num[nap[ijk]][ijk] = num[l][ks+i*ny*nz+j*nz];
      eta_t[nap[ijk]][ijk] = 0.0;
     }
    }   
   
        

   for (l=1;l<n_nx[kn+i*ny*nz+j*nz]+1;++l)
    { if (feta[num[l][kn+i*ny*nz+j*nz]] == 0 ) {
      feta[num[l][kn+i*ny*nz+j*nz]] = 1;
      nap[ijk] = nap[k+i*ny*nz+j*nz] + 1;
      num[nap[ijk]][ijk] = num[l][kn+i*ny*nz+j*nz];
      eta_t[nap[ijk]][ijk] = 0.0;
     }
    } 
  
   
     for (l=1;l<n_nx[k+iw*ny*nz+j*nz]+1;++l){
     ipeta[num[l][k+iw*ny*nz+j*nz]] = eta[l][k+iw*ny*nz+j*nz];
    }
     for (l=1;l<n_nx[k+ie*ny*nz+j*nz]+1;++l){
     ineta[num[l][k+ie*ny*nz+j*nz]] = eta[l][k+ie*ny*nz+j*nz];
    }
     for (l=1;l<n_nx[k+i*ny*nz+js*nz]+1;++l){
     jpeta[num[l][k+i*ny*nz+js*nz]] = eta[l][k+i*ny*nz+js*nz];
    }
     for (l=1;l<n_nx[k+i*ny*nz+jn*nz]+1;++l){
     jneta[num[l][k+i*ny*nz+jn*nz]] = eta[l][k+i*ny*nz+jn*nz];
    }
     for (l=1;l<n_nx[ks+i*ny*nz+j*nz]+1;++l){
     kpeta[num[l][ks+i*ny*nz+j*nz]] = eta[l][ks+i*ny*nz+j*nz];
    }
     for (l=1;l<n_nx[kn+i*ny*nz+j*nz]+1;++l){
     kneta[num[l][kn+i*ny*nz+j*nz]] = eta[l][kn+i*ny*nz+j*nz];
    }
 /*    for (l=1;l<n_nx[k+iw*ny*nz+js*nz]+1;++l){
     ipjpeta[num[l][k+iw*ny*nz+js*nz]] = eta[l][k+iw*ny*nz+js*nz];
    }
     for (l=1;l<n_nx[k+iw*ny*nz+jn*nz]+1;++l){
     ipjneta[num[l][k+iw*ny*nz+jn*nz]] = eta[l][k+iw*ny*nz+jn*nz];
    }
     for (l=1;l<n_nx[k+ie*ny*nz+js*nz]+1;++l){
     injpeta[num[l][k+ie*ny*nz+js*nz]] = eta[l][k+ie*ny*nz+js*nz];
    }
     for (l=1;l<n_nx[k+ie*ny*nz+jn*nz]+1;++l){
     injneta[num[l][k+ie*ny*nz+jn*nz]] = eta[l][k+ie*ny*nz+jn*nz];
    }
*/
    for( l=1;l<nap[ijk]+1;++l)
     		{
     li = ipeta[num[l][ijk]] + ineta[num[l][ijk]];
    
     li = li+jpeta[num[l][ijk]] + jneta[num[l][ijk]];
     li = li + kpeta[num[l][ijk]] + kneta[num[l][ijk]];
     li = li-6.0*eta[l][ijk];
  /*   
     li2 = ipjpeta[num[l][ijk]]+ ipjneta[num[l][ijk]];
     li2 = li2 + injpeta[num[l][ijk]] + injneta[num[l][ijk]];
     li2 = li2-4.0*eta[l][ijk];
    */ 
     laplacian_n = li;
    // laplacian_nn = li2;
    
 //To Protect Overflow
     etatemp = eta[l][ijk];
     setatemp = s_eta[ijk];

 // Allen-Cahn equation
     temp = etatemp*alpha-beta*etatemp*etatemp*etatemp-2*gamma*etatemp*(setatemp-etatemp*etatemp) - 2.0*beta_s*etatemp*creal(comp[ijk])*creal(comp[ijk]);
  
     eta_t[l][ijk] = etatemp+dt*L*(temp+kappa*laplacian_n/(dx*dx));//+kappa*laplacian_nn/(4.0*dx*dx));
    }
    
    
    for( l=1;l<nap[ijk]+1;++l)
     		{
     feta[num[l][ijk]] = 0;    
    } 
    
    for( l=1;l<n_nx[k+iw*ny*nz+j*nz]+1;++l)
     		{
     ipeta[num[l][k+iw*ny*nz+j*nz]] = 0;    
    } 
    
    for( l=1;l<n_nx[k+ie*ny*nz+j*nz]+1;++l)
     		{
     ineta[num[l][k+ie*ny*nz+j*nz]] = 0;    
    } 
    

    for( l=1;l<n_nx[k+i*ny*nz+js*nz]+1;++l)
     		{
     jpeta[num[l][k+i*ny*nz+js*nz]] = 0;    
    } 
    
    for( l=1;l<n_nx[k+i*ny*nz+jn*nz]+1;++l)
     		{
     jneta[num[l][k+i*ny*nz+jn*nz]] = 0;    
    } 

for( l=1;l<n_nx[ks+i*ny*nz+j*nz]+1;++l)
     		{
     kpeta[num[l][ks+i*ny*nz+j*nz]] = 0;    
    } 
    
    for( l=1;l<n_nx[kn+i*ny*nz+j*nz]+1;++l)
     		{
     kneta[num[l][kn+i*ny*nz+j*nz]] = 0;    
    } 
    //clear
   
    
 /*   
    
    for( l=1;l<n_nx[k+iw*ny*nz+js*nz]+1;++l)
     		{
     ipjpeta[num[l][k+iw*ny*nz+js*nz]] = 0;    
    }

 for( l=1;l<n_nx[k+iw*ny*nz+jn*nz]+1;++l)
     		{
     ipjneta[num[l][k+iw*ny*nz+jn*nz]] = 0;    //
    } 

	for( l=1;l<n_nx[k+ie*ny*nz+js*nz]+1;++l)
     		{
     injpeta[num[l][k+ie*ny*nz+js*nz]] = 0;    //
    } 
    	
	for( l=1;l<n_nx[k+ie*ny*nz+jn*nz]+1;++l)
     		{
     injneta[num[l][k+ie*ny*nz+jn*nz]] = 0;    //
    }  
*/
 }//i 
 }//j
}//k  
if (istep > 0)
{


 if(istep < 100){
	mobility = 0.1;}
	else{
	mobility = 100.0;
	}
	for(i1=0; i1<nx; ++i1)
		{
		for(i2=0; i2<ny; ++i2)
		{
			for(i3=0; i3 < nz; ++i3)
			{
			ijk=i3+i1*ny*nz+i2*nz;

			realc = creal(comp[ijk]);
		

/*calculating g*/
        g_bub=(realc*realc*realc-realc);
        g_int=2*beta_s*realc*s_eta[ijk];
        g[ijk]=(g_bub+g_int)+_Complex_I*0.0;

		M_var[ijk] = M  + mobility*16.0*creal(comp[ijk])*creal(comp[ijk])*(1-creal(comp[ijk]))*(1-creal(comp[ijk]));

		
		}
		}
		}


	max = M_var[0];
	for(i1=0; i1<nx; ++i1)
		{
		for(i2=0; i2<ny; ++i2)
		{
			for(i3=0; i3 < nz; ++i3)
			{
			ijk=i3+i1*ny*nz+i2*nz;
			tmpelement = M_var[ijk];
			if(tmpelement>max)
			max = tmpelement;

		}
		}
		}

	min = M_var[0];
	for(i1=0; i1<nx; ++i1)
		{
		for(i2=0; i2<ny; ++i2)
		{
			for(i3=0; i3 < nz; ++i3)
			{
			ijk=i3+i1*ny*nz+i2*nz;
			tmpelement = M_var[ijk];
		
			if(tmpelement<min)
			min = tmpelement;


		}
		}
		}

AA = (min + max)/2.0;

//printf("%d\n",a);

  //evolution of c starts 
   fftw_execute_dft(planF,comp,comp);
    fftw_execute_dft(planF,g,g);
    fftw_execute_dft(planF,comp_partx,comp_partx);
    fftw_execute_dft(planF,comp_party,comp_party);
    fftw_execute_dft(planF,comp_partz,comp_partz);


for(i1=0; i1 < nx; ++i1)
    {
    if(i1 < half_nx) kx = i1*delta_kx;
    else kx = (i1-nx)*delta_kx;
    kx2 = kx*kx;
    for(i2=0; i2 < ny; ++i2)
    {
    if(i2 < half_ny) ky = i2*delta_ky;
    else ky = (i2-ny)*delta_ky;
    ky2 = ky*ky;
    for(i3=0; i3 < nz; ++i3)
    {
    if(i3 < half_nz) kz = i3*delta_kz;
    else kz = (i3-nz)*delta_kz;
    kz2 = kz*kz;
	
    k2 = kx2+ky2+kz2;
    k4 = k2*k2;
		ijk=i3+i1*ny*nz+i2*nz;

	comp_partx[ijk] =  _Complex_I*kx*((g[ijk]) + kappa*k2*comp[ijk]);
	comp_party[ijk] =  _Complex_I*ky*((g[ijk]) + kappa*k2*comp[ijk]);
	comp_partz[ijk] =  _Complex_I*kz*((g[ijk]) + kappa*k2*comp[ijk]);
	

	}
	}
	}


 fftw_execute_dft(planB,comp_partx,comp_partx);
    fftw_execute_dft(planB,comp_party,comp_party);
    fftw_execute_dft(planB,comp_partz,comp_partz);
    
    
	for(i1=0; i1 < nx; ++i1)
        {
	for(i2=0;i2<ny; ++i2)
	{
	for(i3=0;i3<nz; ++i3)
	{
	
	ijk=i3+i1*ny*nz+i2*nz;
        comp_partx[ijk] = (comp_partx[ijk])*(inv_nx*inv_ny*inv_nz) ;
        comp_party[ijk] = (comp_party[ijk])*(inv_nx*inv_ny*inv_nz) ;
        comp_partz[ijk] = (comp_partz[ijk])*(inv_nx*inv_ny*inv_nz) ;
	      
	}
	}
	}
	

	for(i1=0; i1 < nx; ++i1)
   	 {
    
   	for(i2=0; i2 < ny; ++i2)
   	 {
        for(i3=0; i3 < nz; ++i3)
   	 {
   		 ijk=i3+i1*ny*nz+i2*nz;
		
		
		comp_partx[ijk] = M_var[ijk]*comp_partx[ijk];
		comp_party[ijk] = M_var[ijk]*comp_party[ijk];
		comp_partz[ijk] = M_var[ijk]*comp_partz[ijk];
 	  }
	  }
   	  }


    fftw_execute_dft(planF,comp_partx,comp_partx);
    fftw_execute_dft(planF,comp_party,comp_party);
    fftw_execute_dft(planF,comp_partz,comp_partz);

for(i1=0; i1 < nx; ++i1)
    {
    if(i1 < half_nx) kx = i1*delta_kx;
    else kx = (i1-nx)*delta_kx;
    kx2 = kx*kx;
    for(i2=0; i2 < ny; ++i2)
    {
    if(i2 < half_ny) ky = i2*delta_ky;
    else ky = (i2-ny)*delta_ky;
    ky2 = ky*ky;
    for(i3=0; i3 < nz; ++i3)
    {
    if(i3 < half_nz) kz = i3*delta_kz;
    else kz = (i3-nz)*delta_kz;
    kz2 = kz*kz;
	
    k2 = kx2+ky2+kz2;
    k4 = k2*k2;
		ijk=i3+i1*ny*nz+i2*nz;

		inv_denom = 1.0 + AA*kappa*k4*dt;
		inv_denom = 1.0/inv_denom;
		comp[ijk] = inv_denom*((1.0 + AA*kappa*k4*dt)*comp[ijk] + 	_Complex_I*dt*(kx*comp_partx[ijk] + ky*comp_party[ijk] + kz*comp_partz[ijk]) );
	
	}
	}
	}
  // evolution of c ends


    fftw_execute_dft(planB,comp,comp);

for(i1=0; i1 < nx; ++i1)
        {
	for(i2=0;i2<ny; ++i2)
	{
	for(i3=0;i3<nz; ++i3)
	{
	
	ijk=i3+i1*ny*nz+i2*nz;
        comp[ijk] = creal(comp[ijk]*(inv_nx*inv_ny*inv_nz)) +_Complex_I*0.0;
        
        }
	}
	}
}

  for(i=0; i < nx; ++i)
		{
		for(j=0; j < ny; ++j)
			{
for(k=0; k < nz; ++k)
			{
			
			ijk=k+i*ny*nz+j*nz;
    n_nx[ijk] = 0;
    for( l=1;l<nap[ijk]+1;++l)
     		{
     if (fabs(eta_t[l][ijk]) > 0.001) 
      {n_nx[ijk] = n_nx[ijk]+1;
      eta[n_nx[ijk]][ijk]= eta_t[l][ijk];
      num[n_nx[ijk]][ijk] = num[l][ijk];
     }
     } //l
    } //i
   } //j
  
}
 
 if(/*istep > 5080 &&*/ istep%ostep ==0)
		{
		printf("%d\n",(int)(istep));
 /*
   sprintf(NAME,"output/data/comp%d.dat",(int) (istep));
		writepm3d(comp,NAME,nx,ny,nz);
		sprintf(NAME,"output/dataslice/comp%d.dat",(int) (istep));
		writeonedxprofile(comp,NAME,nx,ny,nz);
*/	


  for(i=0; i < nx; ++i)
		{
		for(j=0; j < ny; ++j)
			{
for(k=0; k < nz; ++k)
			{
			
			ijk=k+i*ny*nz+j*nz;
  s_eta[ijk] = 0.0;
  sig_eta[ijk] = 0.0;
  maxidx[ijk] = 0;
}}}

  for(i=0; i < nx; ++i)
		{
		for(j=0; j < ny; ++j)
			{
for(k=0; k < nz; ++k)
			{
			
			ijk=k+i*ny*nz+j*nz;
	  for( l=1;l<n_nx[ijk]+1;++l)
     		{
	eta[l][ijk] = fabs(eta[l][ijk]);
        s_eta[ijk] = s_eta[ijk]+eta[l][ijk]*eta[l][ijk];
        sig_eta[ijk] = sig_eta[ijk]+eta[l][ijk];
        if(s_eta[ijk] > 1.0) {s_eta[ijk] = 1.0;}
    		}
   }
  } 

}
  

if (istep%50000==0)
{
sprintf(NAME,"output/grain/eta_3d%d.dat",(int) (istep));
                if( (fpw = fopen(NAME,"w")) == NULL){
	printf("Unable to open file %s",NAME);
	printf("Exiting\n");
	exit(0);
	}
else{
	fpw = fopen(NAME,"w");
	}
 
}
for(i=0; i < nx; ++i)
		{
		for(j=0; j < ny; ++j)
			{
for(k=0; k < nz; ++k)
			{
			
			ijk=k+i*ny*nz+j*nz;
maxeta = fabs(eta[1][ijk]);		 
		for( l=2;l<n_nx[ijk]+1;++l)
     		{
		 if (fabs(eta[l][ijk])>maxeta)
	{        maxeta=fabs(eta[l][ijk]);
        	} 
		}

     for( l=1;l<nap[ijk]+1;++l)
     		{
		if (istep%50000==0)
		{
		fprintf(fpw,"%d\t%d\t%d\t%d\t%lf\t%lf\n",i,j,k,num[l][ijk],eta[l][ijk],s_eta[ijk]);
		} 


     if(maxeta==fabs(eta[l][ijk]))
       {maxidx[ijk] = num[l][ijk];

      }
     }//l

}}}


 		sprintf(NAME,"output/grain/3d_grain%d.dat",(int) (istep));
		writepm3d_gg(comp,maxidx,s_eta,sig_eta,NAME,nx,ny,nz);
		sprintf(NAME,"output/grain/2d_grain%d.dat",(int) (istep));
		writeonedyprofile(comp,maxidx,s_eta,sig_eta,NAME,nx,ny,nz);


   

for(i=0; i < nx; ++i)
		{
		for(j=0; j < ny; ++j)
			{
for(k=0; k < nz; ++k)
			{
			
			ijk=k+i*ny*nz+j*nz;
for( l=1;l<n_nx[ijk]+1;++l)
     		{

ns = 0;
nsize[l]=0;

}}}}
for(i=0; i < nx; ++i)
		{
		for(j=0; j < ny; ++j)
			{
for(k=0; k < nz; ++k)
			{
			
			ijk=k+i*ny*nz+j*nz;

		maxeta = fabs(eta[1][ijk]);		 
		for( l=1;l<n_nx[ijk]+1;++l)
     		{
		 if (fabs(eta[l][ijk])>maxeta)
	{        maxeta=fabs(eta[l][ijk]);
        	} 
		}
    for( l=1;l<n_nx[ijk]+1;++l)
     		{
       
         
          //if (maxeta==eta[l][ijk]) {

            if (maxeta > 0.95) {
             nsize[num[l][ijk]]=nsize[num[l][ijk]]+1;
             if (nsize[num[l][ijk]]==1) {
              ns = ns+1;
             }
            }
          // }
         }//l 
      
}}}
printf("%d\n",ns);
sprintf(NAME,"output/grn%d.dat",istep);
 if( (fpw = fopen(NAME,"w")) == NULL){
	printf("Unable to open file %s",NAME);
	printf("Exiting\n");
	exit(0);
	}
else{
	fpw = fopen(NAME,"w");
	}
 for( l=1;l<no;++l)
     	{
if (nsize[l]!=0)
		{
		fprintf(fpw,"%d\t%d\n",l,nsize[l]);
		}
		}	
	
	

 fclose(fpw);



 
}//mod_nstep/istep if



}// time loop ends istep 
 


free(n_nx); 
 free(feta);
 free(ipeta);
 free(ineta);
 free(jpeta);
 free(jneta);
 free(kpeta);
 free(kneta);
/*
 free(ipjpeta);
 free(ipjneta);
 free(injpeta);
 free(injneta);
 */
 free(nsize);

 free(s_eta);
free(sig_eta);
 free(maxidx);
 free(M_var);
for( l=1;l<nmax;++l)
     		{
 free(eta[l]);
 free(eta_t[l]);
 free(num[l]);
}
free(nap);
fftw_free(g);
fftw_free(comp);
fftw_free(comp_partx);
fftw_free(comp_partz);
fftw_free(comp_party);


fftw_destroy_plan(planF);
fftw_destroy_plan(planB);



 return(0);
} //main


