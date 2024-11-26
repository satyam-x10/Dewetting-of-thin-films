#include <stdio.h>
#include <stdlib.h>
#include<complex.h>
#include<fftw3.h>

void writeonedxprofile(fftw_complex *array,char NAME[],int n_x,int n_y,int n_z)
{int i,j,ijk,k;
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
j=n_y/2;
 for(i=0;i<n_x;i++)
	{
 	//for(j=0;j<n_y;j++)
	//	{
		for(k=0;k<n_z;k++)
		{
		ijk=k+j*n_z+n_z*n_y*i;
 		tmpelement=creal(array[ijk]);
		//if(j==n_y/2)
		fprintf(fpw,"%d\t%d\t%lf\n",i,k,tmpelement);
		}fprintf(fpw,"\n");
		}
	//}
 fclose(fpw);
}

void writeonedyprofile(fftw_complex *array3, int *array,double *array1,double *array2,char NAME[],int n_x,int n_y,int n_z)
{int i,j,ijk,k;
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
		k=n_z/2;
		
		ijk=k+j*n_z+n_z*n_y*i;
		tmpelement=creal(array3[ijk]);
		fprintf(fpw,"%d\t%d\t%lf\t%d\t%lf\t%lf\n",i,j,tmpelement,array[ijk],array1[ijk],array2[ijk]);

		}		fprintf(fpw,"\n");
	}
 fclose(fpw);
}
void writeonedprofile(fftw_complex *array,char NAME[],int n_x)
{int i,j;
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

 		tmpelement=creal(array[i]);
		fprintf(fpw,"%d\t%lf\n",i,tmpelement);
		
	}

fclose(fpw);
}
