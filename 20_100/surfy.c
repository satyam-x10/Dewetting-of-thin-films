#include <stdio.h>

#include <stdlib.h>

#include <math.h>

int main()
{


int ijk,i,i1,i2,i3,nx=512,ny=512,nz=48,list=330050,interval=10000;
FILE *fpr,*fpw;
double tmpelement;
char ch[50],NAME[50];

double *comp = (double *)malloc(nx*ny*nz* sizeof(double));
double *eta0 = (double *)malloc(nx*ny*nz* sizeof(double));
double *eta1 = (double *)malloc(nx*ny*nz* sizeof(double));
double *eta2 = (double *)malloc(nx*ny*nz* sizeof(double));
double *eta3 = (double *)malloc(nx*ny*nz* sizeof(double));
double *z = (double *)malloc(nx*ny*nz* sizeof(double));
double *x = (double *)malloc(nx*ny*nz* sizeof(double));
double *y = (double *)malloc(nx*ny*nz* sizeof(double));



i=0;
while(i<list)
{
	sprintf(ch,"output/grain/3d_grain%d.dat",i);

 if( (fpr = fopen(ch,"r")) == NULL){
	printf("Unable to open file %s",ch);
	printf("Exiting\n");
	exit(0);
	}
else{
	fpw = fopen(NAME,"r");
	}

for(i1=0; i1<nx; ++i1)
{
for(i2=0; i2<ny; ++i2)
{
for(i3=0; i3<nz; ++i3)
{
ijk = i3 + i2*nz + i1*nz*ny;
fscanf(fpr,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &x[ijk],&x[ijk],&y[ijk],&x[ijk],&comp[ijk],&x[ijk]);
//printf("%lf",comp[ijk]);
}
}
}
fclose(fpr);

	sprintf(NAME,"ydata/2d_grain%d.dat",i);
 if( (fpw = fopen(NAME,"w")) == NULL){
	printf("Unable to open file %s",NAME);
	printf("Exiting\n");
	exit(0);
	}
else{
	fpw = fopen(NAME,"w");
	}


//for(i2=0; i2<ny; ++i2)
//{
for(i1=0; i1<nx; ++i1)
{  
for(i2=0; i2<ny; ++i2)
{
i3 = nz/2 + 3;

		ijk = i3 + i2*nz + i1*nz*ny;
 		tmpelement=comp[ijk];
		fprintf(fpw,"%d\t%d\t%lf\n",i1,i2,tmpelement);
		}fprintf(fpw,"\n");
		//}
	}
 fclose(fpw);

/*
sprintf(ch,"../data/eta0_%d.dat",i);
 if( (fpr = fopen(ch,"r")) == NULL){
	printf("Unable to open file %s",ch);
	printf("Exiting\n");
	exit(0);
	}
else{
	fpw = fopen(NAME,"r");
	}

for(i1=0; i1<nx; ++i1)
{
for(i2=0; i2<ny; ++i2)
{
for(i3=0; i3<nz; ++i3)
{
ijk = i3 + i2*nz + i1*nz*ny;
fscanf(fpr,"%lf\t%lf\t%lf\t%lf\n", &x[ijk],&y[ijk],&z[ijk],&comp[ijk]);
//printf("%lf",comp[ijk]);
}
}
}
fclose(fpr);

	sprintf(NAME,"../slicey/eta0_%d.dat",i);
 if( (fpw = fopen(NAME,"w")) == NULL){
	printf("Unable to open file %s",NAME);
	printf("Exiting\n");
	exit(0);
	}
else{
	fpw = fopen(NAME,"w");
	}

for(i1=0; i1<nx; ++i1)
{
for(i2=0; i2<ny; ++i2)
{
for(i3=0; i3<nz; ++i3)
{
		ijk = i3 + i2*nz + i1*nz*ny;
 		tmpelement=comp[ijk];
		if(i3==63&&i2==116)
		fprintf(fpw,"%d\t%lf\n",i1,tmpelement);
		}
		}
	}
 fclose(fpw);


sprintf(ch,"../data/eta1_%d.dat",i);
 if( (fpr = fopen(ch,"r")) == NULL){
	printf("Unable to open file %s",ch);
	printf("Exiting\n");
	exit(0);
	}
else{
	fpw = fopen(NAME,"r");
	}

for(i1=0; i1<nx; ++i1)
{
for(i2=0; i2<ny; ++i2)
{
for(i3=0; i3<nz; ++i3)
{
ijk = i3 + i2*nz + i1*nz*ny;
fscanf(fpr,"%lf\t%lf\t%lf\t%lf\n", &x[ijk],&y[ijk],&z[ijk],&comp[ijk]);
//printf("%lf",comp[ijk]);
}
}
}
fclose(fpr);

	sprintf(NAME,"../slicey/eta1_%d.dat",i);
 if( (fpw = fopen(NAME,"w")) == NULL){
	printf("Unable to open file %s",NAME);
	printf("Exiting\n");
	exit(0);
	}
else{
	fpw = fopen(NAME,"w");
	}

for(i1=0; i1<nx; ++i1)
{
for(i2=0; i2<ny; ++i2)
{
for(i3=0; i3<nz; ++i3)
{
		ijk = i3 + i2*nz + i1*nz*ny;
 		tmpelement=comp[ijk];
		if(i3==63&&i2==116)
		fprintf(fpw,"%d\t%lf\n",i1,tmpelement);
		}
		}
	}
 fclose(fpw);


sprintf(ch,"../data/eta2_%d.dat",i);
 if( (fpr = fopen(ch,"r")) == NULL){
	printf("Unable to open file %s",ch);
	printf("Exiting\n");
	exit(0);
	}
else{
	fpw = fopen(NAME,"r");
	}

for(i1=0; i1<nx; ++i1)
{
for(i2=0; i2<ny; ++i2)
{
for(i3=0; i3<nz; ++i3)
{
ijk = i3 + i2*nz + i1*nz*ny;
fscanf(fpr,"%lf\t%lf\t%lf\t%lf\n", &x[ijk],&y[ijk],&z[ijk],&comp[ijk]);
//printf("%lf",comp[ijk]);
}
}
}
fclose(fpr);

	sprintf(NAME,"../slicey/eta2_%d.dat",i);
 if( (fpw = fopen(NAME,"w")) == NULL){
	printf("Unable to open file %s",NAME);
	printf("Exiting\n");
	exit(0);
	}
else{
	fpw = fopen(NAME,"w");
	}

for(i1=0; i1<nx; ++i1)
{
for(i2=0; i2<ny; ++i2)
{
for(i3=0; i3<nz; ++i3)
{
		ijk = i3 + i2*nz + i1*nz*ny;
 		tmpelement=comp[ijk];
		if(i3==63&&i2==116)
		fprintf(fpw,"%d\t%lf\n",i1,tmpelement);
		}
		}
	}
 fclose(fpw);

sprintf(ch,"../data/eta3_%d.dat",i);
 if( (fpr = fopen(ch,"r")) == NULL){
	printf("Unable to open file %s",ch);
	printf("Exiting\n");
	exit(0);
	}
else{
	fpw = fopen(NAME,"r");
	}

for(i1=0; i1<nx; ++i1)
{
for(i2=0; i2<ny; ++i2)
{
for(i3=0; i3<nz; ++i3)
{
ijk = i3 + i2*nz + i1*nz*ny;
fscanf(fpr,"%lf\t%lf\t%lf\t%lf\n", &x[ijk],&y[ijk],&z[ijk],&comp[ijk]);
//printf("%lf",comp[ijk]);
}
}
}
fclose(fpr);

	sprintf(NAME,"../slicey/eta3_%d.dat",i);
 if( (fpw = fopen(NAME,"w")) == NULL){
	printf("Unable to open file %s",NAME);
	printf("Exiting\n");
	exit(0);
	}
else{
	fpw = fopen(NAME,"w");
	}

for(i1=0; i1<nx; ++i1)
{
for(i2=0; i2<ny; ++i2)
{
for(i3=0; i3<nz; ++i3)
{
		ijk = i3 + i2*nz + i1*nz*ny;
 		tmpelement=comp[ijk];
		if(i3==63&&i2==116)
		fprintf(fpw,"%d\t%lf\n",i1,tmpelement);
		}
		}
	}
 fclose(fpw);

*/
i=i+interval;

}

}
