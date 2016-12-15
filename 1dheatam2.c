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
#define radius 0.001
double pi=3.14;
int no_moles=3; 
double L=100.000;
double timestep=0.001, time=0.0;
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
	if(i%2==0) m[i]=1.0;
	else m[i]=1.078;
}
initialize_pos(pos, x);
FILE* Inivel=fopen("RandVel1_5.txt", "r");
for(i=0;i<no_moles;i++)
{
	v[i]=0.0;
	if(fscanf(Inivel, "%lf ", &velocity)!=EOF)
	{
	v[i]=velocity;
	fprintf(vel, "\n%d %lf ", i, v[i]);
	*K+=0.5*m[i]*v[i]*v[i];	
	}
}
v[0]=2.0;
v[no_moles-1]=-4.0;
fprintf(vel, "\n");
fclose(Inivel);
return;

}

void calc_timestep(double m[no_moles], double x[no_moles], double v[no_moles])
{
	int i=0; 
	double t[no_moles+1], t_min=1000.0;
	if(v[0]<0)
	{
		if(x[0]>radius)
		{
			t[0]=fabs(x[0]/v[0]);
			if (t_min>t[0])
			t_min=t[0];
		}
		else
		{
			printf("collision");
			//exit(0);
		}
	}
	else t[0]=1000.0;//any large number
	if(v[no_moles-1]>0)
	{
		if(L-x[no_moles-1]>radius)
		{		
			t[no_moles]=fabs((L-x[no_moles-1])/v[no_moles-1]);
			if (t_min>t[no_moles])
			t_min=t[no_moles];
		}
		else
		{
			printf("collision");
			//exit(0);
		}
	}
	else t[no_moles]=1000.0;
	for(i=1;i<no_moles;i++)
	{
		int j=i-1;
		if(fabs(x[j+1]-x[j])<2*radius)
		{
			double u1=v[j];
			double u2=v[j+1];
			v[j]=(m[j]-m[j+1])/(m[j]+m[j+1])*u1+2*m[j+1]/(m[j]+m[j+1])*u2;
			v[i+1]=2*m[j]/(m[j]+m[j+1])*u1-(m[j]-m[j+1])/(m[j]+m[j+1])*u2;
		}
		if(v[i-1]>v[i]&&v[i]*v[i-1]>0)
		{
			t[i]=(fabs(x[i]-x[i-1])-2*radius)/fabs(v[i]-v[i-1]);
		}
		
		else if(v[i-1]>0&&v[i]<0)
		{
			t[i]=fabs(x[i]-x[i-1]-2*radius)/fabs(-v[i]+v[i-1]);
		}
		
		else if(v[i-1]<=0&&v[i]>=0)
		t[i]=1000.0;
		else t[i]=1000.0;
		
		if (t[i]<t_min&&t[i]>0.00001)
		t_min=t[i];
	}
	timestep=t_min;
	return;
}


double position(FILE*RandVelT1, FILE*RandVelT2, FILE* pos, double m[no_moles], double x[no_moles], double v[no_moles], double *time) 
{
	int i, j=0, flag=0;
	double randvelocity;
	fprintf(pos, "\nTime=%lf\n", *time);
	double t, u[no_moles];
	flag=0;
	for(i=0;i<no_moles;i++)
	{
		u[i]=v[i];
	}

/*	
	for(i=0;i<no_moles;i++)
    {         
		x[i]+=u[i]*timestep;
		if(x[i]>150.0||x[i]<-50.0)
		exit(0);
		if(x[i]+radius>=L)
		{
			printf("\n Trying to cross!!!\n");	
			flag++;			
		}	
	}
*/
	
	double dt = 100.0,dt_temp = 100.0;
	int i_temp=0,j_temp=0;

	if(v[0] < 0.0)		dt = -1.0*(x[0])/v[0];
	if(v[no_moles - 1] > 0.0)	dt_temp = (L - x[no_moles-1])/v[no_moles-1];
	
	if(dt_temp < dt) 
	{
		dt = dt_temp;
		flag=-2;//collision with right wall has least time
	}
	else flag=-1;
	for(i=0;i<no_moles-1;i++)
	{
		 if(v[i] > v[i+1])	 dt_temp = -1.0*(x[i+1]-x[i])/(v[i+1] - v[i]);
		 if(dt_temp > 0.0)
		 {
			 if(dt_temp < dt) 
			 {
				 dt = dt_temp;
				 i_temp = i;
				 j_temp = i+1;
				 flag=1;
			 }
		 }
	}

	//printf(" bEFORE: i_temp: %d j_temp:%d x[i_temp]=%lf x[j_temp]=%lf, v[i_temp]=%lf, v[j_temp]=%lf dt = %lf\n", i_temp, j_temp, x[i_temp], x[j_temp], v[i_temp],v[j_temp], dt);

	
	for(i=0;i<no_moles;i++)
	{
		if(flag>0)
		{
			x[i] += u[i]*dt;
			if(i == i_temp)
			{
				v[i]=(m[i]-m[i+1])/(m[i]+m[i+1])*u[i]+2.0*m[i+1]/(m[i]+m[i+1])*u[i+1];
				v[i+1]=2.0*m[i]/(m[i]+m[i+1])*u[i]-(m[i]-m[i+1])/(m[i]+m[i+1])*u[i+1];
			}
		}
		else x[i] += u[i]*dt;
	}
	
	//printf(" After: i_temp:%d j_temp:%d x[i_temp]=%lf x[j_temp]=%lf, v[i_temp]=%lf, v[j_temp]=%lf dt = %lf\n", i_temp, j_temp, x[i_temp], x[j_temp], v[i_temp],v[j_temp], dt);
	
/*	
	for(i=0; i<no_moles-1; i++) 
	{ 
		if(x[i+1]-radius<=x[i]+radius)
		{
		collision++;
		x[i]-=u[i]*timestep;
		x[i+1]-=u[i+1]*timestep;
		t=fabs((2.0*radius+x[i+1]-x[i])/(u[i]-u[i+1]));//time before collision
		x[i]=x[i]+u[i]*t;
		x[i+1]=x[i];
		v[i]=(m[i]-m[i+1])/(m[i]+m[i+1])*u[i]+2.0*m[i+1]/(m[i]+m[i+1])*u[i+1];
		v[i+1]=2.0*m[i]/(m[i]+m[i+1])*u[i]-(m[i]-m[i+1])/(m[i]+m[i+1])*u[i+1];
		x[i]+=v[i]*(timestep-t);
		x[i+1]+=v[i+1]*(timestep-t);
		}
	} 
	*/
	
	if(flag<0)
	{
		if(flag==-1)
		{
			if(x[0]-radius<=0)
			{
				boundingl++;
				printf("\nBounding at left. Number:%d\n", boundingl);
				
				if(fscanf(RandVelT1, "%lf ", &randvelocity)!=EOF)
				{
						v[0]=randvelocity; 				
				}
				else
				exit(0);					
				
			}						
		}
		if(flag==-2)
		{
			if(x[no_moles-1]+radius>=L)
			{
				boundingr++;
				printf("\nBounding at right. Number: %d\n", boundingr);
								
				if(fscanf(RandVelT2, "%lf ", &randvelocity)!=EOF)
				{
					v[no_moles-1]=-randvelocity; 
					
				}
				else
				exit(0);
				
			}	
		}
	}
	
	
	for(i=0;i<no_moles;i++)
	fprintf(pos, "%d %lf ", i, x[i]);
	//if(flag>0)
	//calc_timestep(x, v);
	*time=*time+dt;
	return(dt);
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

void TempProfile(double timestep, double Tdt[no_moles], double Edt[no_moles], double xdt[no_moles], double *tottime, double m[no_moles], double x[no_moles], double v[no_moles])
{
	
	int i=0, no_bins=10, bin, j, count[i];
	for(i=0;i<no_moles;i++)
	{
	Tdt[i]=0.0;
	xdt[i]=0.0;
	Edt[i]=0.0;
	for(j=-5;j<=5;j++)
		{
		if(i+j>=0&&i+j<no_moles)
		{
		if( ( x[i+j] <= L*(i+1)/no_moles ) &&( x[i+j] > L*i/no_moles ) )
		{
		Tdt[i]+=0.5*m[i+j]*v[i+j]*v[i+j]/**timestep*//(0.5*1);
		//*tottime+=timestep;
		//xdt[i]=x[i]*timestep+0.5*v[i]*timestep*timestep;
		xdt[i]+=x[i+j];
		Edt[i]+=0.5*m[i+j]*v[i+j]*v[i+j]/**timestep*/;
		count[i]++;
		}
		}
		}
		if(count[i])
		{
		Tdt[i]=Tdt[i]/(double)count[i];
		xdt[i]=xdt[i]/(double)count[i];
		Edt[i]=Edt[i]/(double)count[i];
		}
		
	}
	FILE* tempprof=fopen("tempprofileamass.txt", "w");
	for(i=0;i<no_moles;i++)
	fprintf(tempprof, "%d %lf ", i, Tdt[i]); 
	fclose(tempprof);
	/*	
		
	int count[no_bins];
	double unit=L/no_bins, T[no_bins], E[no_bins];
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
		E[bin]=0;
		
			for(i=0; i<no_moles;i++)			
			{
				if(x[i]>=((double)bin)*unit && x[i]<((double)(bin+1)*unit))
				{
					T[bin]+=m[i]*v[i]*v[i]/k;
					E[bin]+=0.5*m[i]*v[i]*v[i];
					count[bin]++;
				 }
				 
			}		
		
		T[bin]=T[bin]/count[bin];
		E[bin]=E[bin]/count[bin];
		printf("\n%f", T[bin]);

		fprintf(tempprof, "%d %f, ", bin, T[bin]); 
		fprintf(energyprof, "%d %f, ", bin, E[bin]);
	}
	*/
	return;
}
	

void main()
{

double m[no_moles], v[no_moles], x[no_moles], K, Temp, tottime=0.0, Tdt[no_moles], Edt[no_moles], xdt[no_moles];
int i, j;
for(i=0;i<no_moles;i++)
{
	Tdt[i]=0.0;
	xdt[i]=0.0;
	Edt[i]=0.0;
}


/*file naming*/
FILE *pos=fopen("position.txt", "w");
FILE *vel=fopen("velocity.txt", "w");
FILE *ene=fopen("energy1.txt", "w");
FILE *temp=fopen("temperature.txt", "w");
//FILE* RandVelT1=fopen("RandVel1_0.txt", "r");
//FILE* RandVelT2=fopen("RandVel3_0.txt", "r");
FILE* RandVelT1=fopen("RaylRnd.txt", "r");
FILE* RandVelT2=fopen("RaylRnd3.txt", "r");
/*file naming end*/


initialize(pos, vel, m, x, v, &K); 

//calc_timestep(m, x, v);

for (; time<=10.0;time)
{	
	timestep=position(RandVelT1, RandVelT2, pos, m, x, v, &time);	
	velocity(vel, m, v, &K);
	
	Temp=K/(0.5*(double)no_moles*k);
   	printf("Time: %lf K: %lf ",time, K);
   	printf("Temp=%lf K\n", K/(0.5*k*(double)no_moles));
   	fprintf(temp,"\n %lf, %lf", time, K/(0.5*k*(double)no_moles));
	fprintf(ene, "%lf %lf %lf",  time, Temp, K);
	//if(time>=80.0)
	//{
		//TempProfile(timestep, Tdt, Edt, xdt, &tottime, m, x, v);
	//}
	
}
TempProfile(timestep, Tdt, Edt, xdt, &tottime, m, x, v);
printf("\nboundingl=%d, boundingr=%d collision=%d\n", boundingl, boundingr, collision);

fclose(pos);
fclose(vel);
fclose(ene);
fclose(temp);
fclose(RandVelT1);
fclose(RandVelT2);

/*FILE* tempprof=fopen("tempprofileamass.txt", "w");
for(i=0;i<no_moles;i++)
fprintf(tempprof, "%lf %lf %lf ", xdt[i]/tottime, Tdt[i]/tottime, Edt[i]/tottime); 
fclose(tempprof);
*/


printf("\ndone\n flag=%d", flag);
printf("\ntimestep=%lf", timestep);

}
