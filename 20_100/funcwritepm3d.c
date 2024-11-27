#include <stdio.h>
#include <stdlib.h>
#include<complex.h>
#include<fftw3.h>

void writepm3d(fftw_complex *array,char NAME[],int n_x,int n_y,int n_z)
{int i,j,k,ijk;
 double tmpelement;
 FILE *fpw;
 if( (fpw = fopen(NAME,"w")) == NULL){
	printf("Unable to open file %s",NAME);
	printf("Exiting\n");
	exit(0);
	}
else{
	fpw = fopen(NAME,"w");
	}
 	for(i=0;i<n_x;i++)
		{
	for(j=0;j<n_y;j++)
		{
 	for(k=0;k<n_z;k++)
			{
				ijk=k+j*n_z+n_y*n_z*i;
 				tmpelement=creal(array[ijk]);
				fprintf(fpw,"%d\t%d\t%d\t%lf\n",i,j,k,tmpelement);
			}
		
		}
		fprintf(fpw,"\n");
		}

 fclose(fpw);
}

void writepm3d_gg(fftw_complex *array3, int *array,double *array1,double *array2,char NAME[],int n_x,int n_y,int n_z)
{int i,j,ij,k,ijk;
 double tmpelement;
 FILE *fpw;
 if( (fpw = fopen(NAME,"w")) == NULL){
	printf("Unable to open file %s",NAME);
	printf("Exiting\n");
	exit(0);
	}
 for(i=0;i<n_x;i++)
	{
 	for(j=0;j<n_y;j++)
		{
		for(k=n_z/2;k<n_z;k++)
			{
				ij=k+j*n_z+n_y*n_z*i;
				tmpelement=creal(array3[ij]);
		fprintf(fpw,"%d\t%d\t%d\t%lf\t%d\t%lf\t%lf\n",i,j,k,tmpelement,array[ij],array1[ij],array2[ij]);
		}
		fprintf(fpw,"\n");
	}
	}
 fclose(fpw);
}



