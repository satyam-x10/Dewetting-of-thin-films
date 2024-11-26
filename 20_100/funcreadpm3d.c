#include <stdio.h>
#include <stdlib.h>
#include<complex.h>
#include<fftw3.h>
#include "headers/functions.h"
#include "headers/nrutil.h"
void readpm3d(fftw_complex *array,char NAME[],int n_x,int n_y,int n_z)
{int i,j,k,ijk;
 int l,m,n;
 double tmpelement,junk;
 FILE *fpr;
 if( (fpr = fopen(NAME,"r")) == NULL){
	printf("Unable to open file %s",NAME);
	printf("Exiting\n");
	exit(0);
	}
else{
	fpr = fopen(NAME,"r");
	}
 	for(i=0;i<n_x;i++)
		{
	for(j=0;j<n_y;j++)
		{
 	for(k=0;k<n_z;k++)
			{
				ijk=k+j*n_z+n_y*n_z*i;
				tmpelement=0.0;
				fscanf(fpr,"%d%d%d%lf",&l,&m,&n,&tmpelement);
 				array[ijk]=tmpelement+_Complex_I*0.0;
			}
		
		}
		//fscanf(fpr,"%lf",&junk);
		}

 fclose(fpr);
}



