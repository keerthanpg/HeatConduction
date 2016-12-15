/*1. initialise positions to forced coordinates
  2. initialise velocities to uniform random with mean 0, and std_dev sqrt(k*T/m)
  3. calculate accelerations and potentials at time 0, using the positions at time 0
  4. start the time loop(so, you have x(0), v(0) and a(0))
  5. calculate x(t)=x(t-1)+v(t-1)+1/2*a(t-1)*timestep*timestep
  6. apply periodic boundary conditions ie if (x>=l) x=x-l; else if(x<0) x=x+l;
  7. calculate a(t) using x(t)
  8. calculate v(t) using verlet algo ie v(t)=v(t-1)+1/2*(a(t)+a(t-1))
  7. store a(t) in place of a(t-1)
  */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
//#include "ran_uniform.h"

double k=1.0;//boltzmann constant,in SI units
double m=1.0 ;//mass of hydrogen molecule, in kg
double Tl=2.0;//left bath temperature
double Tr=2.0;//right
double sigma=1.0;//for hydrogen
double epsilon=1.0;//for hydrogen
#define dis_cut 2.5*sigma
//#define dis_far_nbour 2.0*dis_cut
#define dis_far_nbour 1000.00
double pi=3.14;
int no_moles=1000; 
double L=100.000;//Number density = 0.728
double timestep=0.1, time=0.0;
int boundingl=0, boundingr=0;

void initialize_pos(FILE* pos, double x[])
{

int i=0;
double x1;

fprintf(pos, "time=0\n");

for(i=0;i<=no_moles;i++)
{
	x[i]=/*L/(2.0*(double)no_moles)*/+ L/((double)no_moles+1)*(double)i;
	fprintf(pos, "\n%d %lf ", i, x[i]);
				
}

fprintf(pos, "\n");
return;

}

void initialize(FILE* pos, FILE* vel, double x[], double v[], double *K) 
{
int i=0, j;
double velocity;

fprintf(vel, "\ntime=0\n");
*K=0.0;
initialize_pos(pos, x);
FILE* Inivel=fopen("RandVel2_0.txt", "r");
for(i=0;i<no_moles;i++)
if(fscanf(Inivel, "%lf ", &velocity)!=EOF)
{
	v[i]=velocity;
	fprintf(vel, "\n%d %lf ", i, v[i]);
	*K+=0.5*m*v[i]*v[i];	
}
fprintf(vel, "\n");
fclose(Inivel);
return;

}
double LeftRecursiveX(FILE* RandVel1_0, double* x, double* v, double t)
{
	double bct, t1, t2;//before collision time
	double randvelocity;
	boundingl++;
	*x-=(*v)*(t);
	
	printf("\n Bounding at left!!! No. %d\n", boundingl);
	bct=fabs((*x)/(*v));
	
	if(fscanf(RandVel1_0, "%lf ", &randvelocity)!=EOF)
		{
			*v=randvelocity; 
			if (*v<0)
			*v=fabs(*v);				
		}
	else
		exit(0);
	*x=0.0;
	*x+=(*v)*(t-bct);
	if(boundingl==4)
	{
	printf("\n %lf, %lf\n ", *x, *v);
	//exit(0);
	
	}
	
	return(t-bct);
}

double RightRecursiveX(FILE* RandVel2_0, double *x, double *v, double t)
{
	double bct, t1, t2;//before collision time
	double randvelocity;
	boundingr++;
	*x-=(*v)*(t);
	printf("\n Bounding at right!!! No. %d\n", boundingr);
	bct=fabs((L-*x)/(*v));
	
	if(fscanf(RandVel2_0, "%lf ", &randvelocity)!=EOF)
		{
			*v=randvelocity; 
			if (*v>0)
			*v=-*v;				
		}
	else
		exit(0);
	*x=L;
	*x+=(*v)*(t-bct);
	printf(" %lf", *x);
	//exit(0);
	return(t-bct);
}

void position(FILE*RandVelT1, FILE*RandVelT2, FILE* pos, double x[no_moles], double v[no_moles], double time) 
{
	int i, j=0;
	double randvelocity;
	fprintf(pos, "\nTime=%lf\n", time);
	double t;
	
	for(i=0;i<no_moles;i++)
    {         
		x[i]+=v[i]*timestep;
		if(x[i]<=0)
		{
			
			boundingl++;
			printf("\nBounding at left. Number:%d\n", boundingl);
			x[i]-=v[i]*timestep;
			t=fabs(x[i]/v[i]);
			if(fscanf(RandVelT1, "%lf ", &randvelocity)!=EOF)
			{
				v[i]=randvelocity; 
				if (v[i]<0)
				v[i]=-v[i];				
			}
			else
			exit(0);
			x[i]=0.0+v[i]*(timestep-t);
			
			/*
			t=timestep;
			do
			{
			t=LeftRecursiveX(RandVelT1,&x[i], &v[i], t);
			}while(x[i]<=0);*/
			
		}
		else if(x[i]>=L)
		{
			boundingr++;
			printf("\nBounding at right. Number: %d\n", boundingr);
			x[i]-=v[i]*timestep;
			t=fabs((L-x[i])/v[i]);
			if(fscanf(RandVelT2, "%lf ", &randvelocity)!=EOF)
			{
				v[i]=randvelocity; 
				if (v[i]>0)
				v[i]=-v[i];				
			}
			x[i]=L+v[i]*(timestep-t);
			/*t=timestep;			
			do
			{
			t=RightRecursiveX(RandVelT2, &x[i], &v[i], t);
			}while(x[i]>=L);	*/
		}		
		
	}
	double temp;
	/*for(i=0;i<no_moles;i++)
	{
		for(j=0;j<no_moles;j++)if(i-1>=0)
		{
			if(x[i]<=x[i-1]&& v[i]*v[i-1]>=0)
			{
				temp=v[i];
				v[i]=v[i-1];
				v[i-1]=temp;
			}
			if(x[i]<=x[i-1]&& v[i]*v[i-1]<0)
			{
				temp=v[i-1];
				v[i-1]=-v[i];
				v[i]=temp;
				temp=x[i-1];
				x[i-1]=x[i];
				x[i]=temp;
			}//cases of collision		
			 
			fprintf(pos, "\n%d %lf ", i, x[i]);
		}
	}*/
	double x_temp;
	for(i=1; i<no_moles; i++) 
	{ 
		x_temp = x[i]; 
		j = i-1; 
		while(x_temp<x[j] && j>=0) 
		{ 
			x[j+1] = x[j];
			j = j-1; 
		} 
		x[j+1] = x_temp;
	} 	
	for(i=0;i<no_moles;i++)
	fprintf(pos, "%d %lf ", i, x[i]);
	return;
}

void velocity(FILE* vel, double v[no_moles], double *K)
{
	int i;	
	*K = 0.0;
	
	fprintf(vel, "\ntime=%lf\n ", time);
	for(i=0;i<no_moles;i++)
    {   
		*K +=0.5*m*v[i]*v[i];			
		fprintf(vel, "%d %lf ", i, v[i]);
	}
	
	return;
}

void TempProfile(FILE*tempprof, double x[no_moles], double v[no_moles])
{
	int i=0, no_bins=10, bin;
	int count[no_bins];
	double unit=L/no_bins, T[no_bins];
	for(i=0;i<no_moles;i++)
	{
		printf(" %f,", x[i]);
	}
	i=0;	
	for(bin=0;bin<no_bins;bin++)
	{
		count[bin]=0;
		T[bin]=0;
		while(x[i]>=((double)bin)*unit && x[i]<((double)(bin+1)*unit))
		{
			T[bin]+=m*v[i]*v[i]/k;
			count[bin]++;
			i++;
		}
		printf("%f\n", T[bin]);
		T[bin]=T[bin]/count[bin];
		fprintf(tempprof, "%d %f, ", bin, T[bin]); 
	}
	
	return;
}
	

void main()
{

double v[no_moles], x[no_moles], K, Temp;
int i, j;


/*file naming*/
FILE *pos=fopen("position.txt", "w");
FILE *vel=fopen("velocity.txt", "w");
FILE *ene=fopen("energy.txt", "w");
FILE *temp=fopen("temperature.txt", "w");
FILE* RandVelT1=fopen("RandVel2_0.txt", "r");
FILE* RandVelT2=fopen("RandVel4_0.txt", "r");
/*file naming end*/


initialize(pos, vel, x, v, &K); 
for(i=0;i<no_moles;i++)
	{
		printf(" %f,", v[i]);
	}

for (time=timestep; time<=1000000*timestep; time=time+timestep)
{	
	position(RandVelT1, RandVelT2, pos, x, v, time);	
	velocity(vel, v, &K);
	
	Temp=K/(0.5*(double)no_moles*k);
   	printf("Time: %lf K: %lf ",time, K);
   	printf("Temp=%lf K\n", K/(0.5*k*(double)no_moles));
   	fprintf(temp,"\n %lf, %lf", time, K/(0.5*k*(double)no_moles));
	fprintf(ene, "\nTime: %lf Temp: %lf K K: %lf",  time, Temp, K);
	
}

printf("\nboundingl=%d, boundingr=%d \n", boundingl, boundingr);

fclose(pos);
fclose(vel);
fclose(ene);
fclose(temp);
fclose(RandVelT1);
fclose(RandVelT2);

FILE* tempprof=fopen("tempprofile.txt", "w");
TempProfile(tempprof, x, v);
fclose(tempprof);


printf("\ndone\n");

}
