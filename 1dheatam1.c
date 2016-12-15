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
#include<ctype.h>
//#include "ran_uniform.h"

double k=1.0;//boltzmann constant,in SI units
double mass=1.0 ;//mass of hydrogen molecule, in kg
double Tl=2.0;//left bath temperature
double Tr=2.0;//right
double sigma=1.0;//for hydrogen
double epsilon=1.0;//for hydrogen
#define radius 0.000001
double pi=3.14;
int no_moles=1000; 
double L=1000.000;
double timestep=0.0001, time=0.0;
int boundingl=0, boundingr=0, collision=0;
int flag =0;

void initialize_pos(FILE* pos, double x[])
{

int i=0;
double x1;

fprintf(pos, "time=0\n");

for(i=0;i<=no_moles;i++)
{
	x[i]=L/(2.0*(double)no_moles)+ L/((double)no_moles)*(double)i;
	fprintf(pos, "\n%d %lf ", i, x[i]);
				
}

fprintf(pos, "\n");
return;

}

void initialize(FILE* pos, FILE* vel, double m[], double x[], double v[], double *K) 
{
int i=0, j;
double velocity;

fprintf(vel, "\ntime=0\n");
*K=0.0;
for(i=0;i<no_moles;i++)
{
	m[i]=mass;
}
initialize_pos(pos, x);
FILE* Inivel=fopen("RandVel1_5.txt", "r");
for(i=0;i<no_moles;i++)
if(fscanf(Inivel, "%lf ", &velocity)!=EOF)
{
	v[i]=velocity;
	fprintf(vel, "\n%d %lf ", i, v[i]);
	*K+=0.5*m[i]*v[i]*v[i];	
}
fprintf(vel, "\n");
fclose(Inivel);
return;

}
void calc_timestep(double x[no_moles], double v[no_moles])
{
	int i=0; 
	double t[no_moles+1], t_min=1000.0;
	if(v[0]<0)
	{
		t[0]=fabs(x[0]/v[0]);
		if (t_min>t[no_moles])
		t_min=t[no_moles];
	}
	else t[0]=1000.0;//any large number
	if(v[no_moles-1]>0)
	{
		t[no_moles]=fabs(x[no_moles-1]/v[no_moles-1]);
		if (t_min>t[no_moles])
		t_min=t[no_moles];
	}
	else t[no_moles]=1000.0;
	for(i=1;i<no_moles;i++)
	{
		if(v[i-1]>v[i]&&v[i]*v[i-1]>0)
		{
			t[i]=fabs((x[i]-x[i-1])/(v[i]-v[i-1]));
		}
		
		else if(v[i-1]>0&&v[i]<0)
		{
			t[i]=fabs((x[i]-x[i-1])/(v[i]+v[i-1]));
		}
		
		else if(v[i-1]<=0&&v[i]>=0)
		t[i]=1000.0;
		else t[i]=1000.0;
		
		if (t[i]<t_min)
		t_min=t[i];
	}
	timestep=t_min;
	return;
}
	
	


void position(FILE*RandVelT1, FILE*RandVelT2, FILE* pos, double m[no_moles], double x[no_moles], double v[no_moles], double time) 
{
	int i, j=0;
	double randvelocity;
	fprintf(pos, "\nTime=%lf\n", time);
	double t, u[no_moles];
	for(i=0;i<no_moles;i++)
	{
		u[i]=v[i];
	}
	
	for(i=0;i<no_moles;i++)
    {         
		x[i]+=u[i]*timestep;
		if(x[i]+2*radius>=L)
		{
			printf("\n Bastards trying to cross!!!\n");	
			flag++;
		}	
	}
	if(x[0]-radius<=0)
	{
		boundingl++;
		printf("\nBounding at left. Number:%d\n", boundingl);
		x[0]-=v[0]*timestep;
		if(v[0]!=0)
		{
			t=fabs(x[0]/v[0]);
			printf("\nBefore collision time =%lf, x[0]=%lf\n ", t, x[0]);
			if(fscanf(RandVelT1, "%lf ", &randvelocity)!=EOF)
			{
			v[0]=randvelocity; 
			if (v[0]<0)
			v[0]=-v[0];				
			}
			else
			exit(0);
			x[0]=0.0+v[0]*(timestep-t);
		}
		else printf("\n v[i]=0");
						
	}
	if(x[no_moles-1]+radius>=L)
	{
		boundingr++;
		printf("\nBounding at right. Number: %d\n", boundingr);
		x[no_moles-1]-=v[no_moles-1]*timestep;
		t=fabs((L-x[no_moles-1])/v[no_moles-1]);
		printf("\nBefore collision time =%lf\n x[no_moles-1]=%lf\n", t, x[no_moles-1]);
		if(fscanf(RandVelT2, "%lf ", &randvelocity)!=EOF)
		{
			v[no_moles-1]=randvelocity; 
			if (v[no_moles-1]>0)
			v[no_moles-1]=-v[no_moles-1];				
		}
		x[no_moles-1]=L+v[no_moles-1]*(timestep-t);
	
	}		
	double temp;
	
	double x_temp,v_temp;
	for(i=0; i<no_moles-1; i++) 
	{ 
		if(x[i+1]-radius<=x[i]+radius)
		{
		collision++;
		x[i]-=u[i]*timestep;
		x[i+1]-=u[i]*timestep;
		t=(2*radius+x[i+1]-x[i])/(u[i]-u[i+1]);//time before collision
		x[i]=x[i]+u[i]*t;
		x[i+1]=x[i];
		v[i]=(m[i]-m[i+1])/(m[i]+m[i+1])*u[i]+2*m[i+1]/(m[i]+m[i+1])*u[i+1];
		v[i+1]=2*m[i]/(m[i]+m[i+1])*u[i]-(m[i]-m[i+1])/(m[i]+m[i+1])*u[i+1];
		x[i]+=v[i]*(timestep-t);
		x[i+1]+=v[i+1]*(timestep-t);
		}
	} 	
	for(i=0;i<no_moles;i++)
	fprintf(pos, "%d %lf ", i, x[i]);
	return;
}

void velocity(FILE* vel, double m[no_moles], double v[no_moles], double *K)
{
	int i;	
	*K = 0.0;
	
	fprintf(vel, "\ntime=%lf\n ", time);
	for(i=0;i<no_moles;i++)
    {   
		*K +=0.5*m[i]*v[i]*v[i];			
		fprintf(vel, "%d %lf ", i, v[i]);
	}
	
	return;
}

void TempProfile(FILE*tempprof, double m[no_moles], double x[no_moles], double v[no_moles])
{
	int i=0, no_bins=10, bin;
	int count[no_bins];
	double unit=L/no_bins, T[no_bins];
	for(i=0;i<no_moles;i++)
	{
		printf(" %f,", x[i]);
	}
	i=0;	
	printf("\n Temperature in bins");
	for(bin=0;bin<no_bins;bin++)
	{
		count[bin]=0;
		T[bin]=0;
		
		if(x[i]<=L&&x[i]>=0)
		{
			while(x[i]>=((double)bin)*unit && x[i]<((double)(bin+1)*unit))
			{
				T[bin]+=m[i]*v[i]*v[i]/k;
				count[bin]++;
				i++; 
			}
			
		}
		else i++;		
		
		printf("\n%f", T[bin]);
		T[bin]=T[bin]/count[bin];
		fprintf(tempprof, "%d %f, ", bin, T[bin]); 
	}
	
	return;
}
	

void main()
{

double m[no_moles], v[no_moles], x[no_moles], K, Temp;
int i, j;


/*file naming*/
FILE *pos=fopen("position.txt", "w");
FILE *vel=fopen("velocity.txt", "w");
FILE *ene=fopen("energy.txt", "w");
FILE *temp=fopen("temperature.txt", "w");
FILE* RandVelT1=fopen("RandVel1_0.txt", "r");
FILE* RandVelT2=fopen("RandVel3_0.txt", "r");
/*file naming end*/


initialize(pos, vel, m, x, v, &K); 

calc_timestep(x, v);
for (time=timestep; time<=10*timestep; time=time+timestep)
{	
	position(RandVelT1, RandVelT2, pos, m, x, v, time);	
	velocity(vel, m, v, &K);
	
	Temp=K/(0.5*(double)no_moles*k);
   	printf("Time: %lf K: %lf ",time, K);
   	printf("Temp=%lf K\n", K/(0.5*k*(double)no_moles));
   	fprintf(temp,"\n %lf, %lf", time, K/(0.5*k*(double)no_moles));
	fprintf(ene, "\nTime: %lf Temp: %lf K K: %lf",  time, Temp, K);
	
}

printf("\nboundingl=%d, boundingr=%d collision=%d\n", boundingl, boundingr, collision);

fclose(pos);
fclose(vel);
fclose(ene);
fclose(temp);
fclose(RandVelT1);
fclose(RandVelT2);

FILE* tempprof=fopen("tempprofile.txt", "w");
TempProfile(tempprof, m, x, v);
fclose(tempprof);

printf("\ndone\n flag=%d", flag);
printf("\ntimestep=%lf", timestep);

}
