#include <stdio.h>

#include <stdlib.h>

#include <math.h>
//#include"funcreadpm3d.c" 


int main()
{

int ijk,i,j,i1,i2,nx=512,ny=512,list=220000,interval=2000;
FILE *fpr,*fpw;
double tmpelement;
char ch[50],NAME[50];
double v_frac[list];
int volume[200000];




double radius[50000];
int *comp = (int *)malloc(nx*ny* sizeof(int));
int tmp;
int x,y;
double m,n;
int a,p,b1,min;
p = 5001;





double r_avg;
double delt=0.01;
double delx = 1.0;
i=2000;
while(i<list)
{
	printf("%d\n",i);
sprintf(ch,"output/grain/2d_grain%d.dat",i);
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

ijk =  i2 + i1*ny;


fscanf(fpr,"%d\t%d\t%lf\t%d\t%lf\t%lf\n", &x,&y,&y,&tmp,&m,&n);
comp[ijk] = tmp;
//printf("%d\n",comp[ijk]);

}
}


for(a=0; a<p; ++a)
{
volume[a] = 0;
}


for(a=0; a<p; ++a)
{
for(i1=0; i1<nx; ++i1)
{
for(i2=0; i2<ny; ++i2)
{

ijk =  i2 + i1*ny;
if (comp[ijk] == a)
{
volume[a] = volume[a] + 1;

}



}
}


}


sprintf(NAME,"output/volume/volume%d.dat",i);
 if( (fpw = fopen(NAME,"w")) == NULL){
	printf("Unable to open file %s",NAME);
	printf("Exiting\n");
	exit(0);
	}
else{
	fpw = fopen(NAME,"w");
	}
b1=0;
r_avg = 0.0;
for(a=0; a<p; ++a)
{

		if(volume[a] >0)
		{b1=b1+1;
		tmpelement = (volume[a]/(M_PI));
			radius[a] =delx* pow(tmpelement,1/2.0);
			r_avg = r_avg + radius[a];
		fprintf(fpw,"%d\t%d\t%d\t%lf\n",b1,a,volume[a],radius[a]);
		}
}
	

 fclose(fpw);

r_avg = 512*512/(b1*M_PI);
r_avg =delx* pow(r_avg,1/2.0);
sprintf(NAME,"avg_rad.dat");
 if( (fpw = fopen(NAME,"a")) == NULL){







	printf("Unable to open file %s",NAME);
	printf("Exiting\n");
	exit(0);
	}
else{
	fpw = fopen(NAME,"a");
	}

		fprintf(fpw,"%lf\t%lf\t%lf\t%d\n",i*delt,r_avg,r_avg*r_avg*r_avg*r_avg,b1);
		
	

 fclose(fpw);
	
//printf("%d",p);
i=i+interval;
}
//I LOOP ENDS
free(comp);




}
//PROGRAM ENDS



