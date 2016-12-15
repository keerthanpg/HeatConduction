/*1. initialise positions to forced coordinates
  2. initialise velocities to uniform random with mean 0, and std_dev sqrt(k*T/m)
  3. calculate accerelations and potentials at time 0, using the positions at time 0
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
double T=1;//at room temperature
double sigma=1.0;//for hydrogen
double epsilon=1.0;//for hydrogen
#define dis_cut 2.5*sigma
#define dis_far_nbour 2.0*dis_cut
double pi=3.14;
int no_moles=1000; 
double lattice=20.000;//Number density = 0.728
double timestep=0.001, time=0.0;
 
double randn (double mu, double sigma)
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;
 
  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }
 
  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);
 
  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;
 
  call = !call;
 
  return (mu + sigma * (double) X1);
}
	
void initialize_pos(FILE* pos, double x[][3])
{
int p=1, q=1, r=1;
int i=0;
double x1, x2, x3;
double v_mean=0;
double v_stddev=sqrt(k*T/m);

fprintf(pos, "time=0\n");

for(p=1;p<=10;p++)
{
	for(q=1;q<=10;q++)
	{
		for(r=1;r<=10;r++)
		{
			x1=1.1*p; x2=1.1*q; x3=1.1*r;
			fprintf(pos, "\nMole%d %lf %lf %lf ", i, x1, x2, x3);
			x[i][0]=x1; x[i][1]=x2; x[i][2]=x3;
			i++;
			if(i>=no_moles)
			{
				fprintf(pos, "\n");
				return;
			}
		}
	}
}
}

void initialize_pos_random(FILE* pos, double x[no_moles][3])
{
	int i, k;
	for(i=0;i<no_moles; i++)
	{
		fprintf(pos, "\nMole%d ", i+1);
		for(k=0;k<3;k++)
		{
			x[i][k]=(double)rand()/(double)RAND_MAX*lattice;
			fprintf(pos, "%lf ", x[i][k]);
		}
		fprintf(pos, "\n");
	}

return;	
}

void initialize(FILE* pos, FILE* vel, double x[][3], double v[][3]) 
{
int i, k;

fprintf(vel, "\ntime=0\n");
double v1, v2, v3, K=0;

initialize_pos(pos, x);


for(i=0;i<no_moles;i++)
{
	v1=randn(0,1); v2=randn(0,1); v3=randn(0,1);
	fprintf(vel, "\nMole%d %lf %lf %lf ", i, v1, v2, v3);
	v[i][0]=v1; v[i][1]=v2; v[i][2]=v3;
	for(k=0;k<3;k++)
	K+=0.5*m*v[i][k]*v[i][k];
	
}

fprintf(vel, "\n");

return;
}

void position(FILE* pos, double x[no_moles][3], double v[no_moles][3], double a[no_moles][3][2], double time) 
{
	int i, k;
	fprintf(pos, "\nTime=%lf\n", time);
	for(i=0;i<no_moles;i++)
    {         
		for(k=0;k<3;k++)
		{			
			x[i][k]+=v[i][k]*timestep+0.5*a[i][k][0]*timestep*timestep;
			if(x[i][k]>=lattice)
				x[i][k]-=lattice;
			else if (x[i][k]<0.0)
				x[i][k]+=lattice;//periodic boundary condition
		}
		fprintf(pos, "\nmole%d %lf %lf %lf ", i, x[i][0], x[i][1], x[i][2]);
	}
	return;
}
	
	
void force(double x[no_moles][3], FILE* acc, double a[no_moles][3][2], double *U, double time)
{
	int i, k, j;
	double x_diff[3], dis;
	
	fprintf(acc, "\ntime=%lf\n ", time); 
	for(i=0;i<no_moles;i++)
	{
		for(k=0;k<3;k++)
		{
			a[i][k][0]=a[i][k][1];
			a[i][k][1]=0;
		}
	}
	*U = 0.0;
	

	for(i=0;i<no_moles;i++)
    {         
        for(j=i+1;j<no_moles;j++)
        {
			for(k=0;k<3;k++)
			{
				x_diff[k]=x[i][k]-x[j][k];
				if (x_diff[k]>lattice/2)
				x_diff[k]-=lattice;
				else if(x_diff[k]<-lattice/2)
				x_diff[k]+=lattice;//closest neighbour condition
				
			}	
			dis=sqrt(pow(x_diff[0], 2)+pow(x_diff[1], 2)+pow(x_diff[2], 2));
			if(dis<dis_far_nbour)              
			{
				*U += 4.0*epsilon*(pow((sigma/dis),12)-pow((sigma/dis),6));
				for(k=0;k<3;k++)
				{					
						a[i][k][1] -= 4.0*epsilon/m*(-12*pow(sigma/dis, 14)+6*pow(sigma/dis, 8))*x_diff[k];
						a[j][k][1] += 4.0*epsilon/m*(-12*pow(sigma/dis, 14)+6*pow(sigma/dis, 8))*x_diff[k];				
				}
			}			
			//potential energy due to interaction of two molecules is
			//assumed to be equally divided between them
        }       
		fprintf(acc, "\nmole%d %lf %lf %lf ", i, a[i][0][1], a[i][1][1], a[i][2][1]);
	}
return;
}

void velocity(FILE* vel, double v[no_moles][3], double a[no_moles][3][2], double *K)
{
	int i,k;
	*K = 0.0;
	fprintf(vel, "\ntime=%lf\n ", time);
	for(i=0;i<no_moles;i++)
    {         	       
        for(k=0;k<3;k++)
	    {
			v[i][k]+=0.5*timestep*(a[i][k][0]+a[i][k][1]);	
			*K +=0.5*m*v[i][k]*v[i][k];		
			a[i][k][0]=a[i][k][1];
		}
		
		fprintf(vel, "\nmole%d %lf %lf %lf ", i, v[i][0], v[i][1], v[i][2]);
	}
	return;
}
	

void main()
{

double v[no_moles][3], x[no_moles][3], U, a[no_moles][3][2], K, T;
int i, j, k;

/*file naming*/
FILE *pos=fopen("position.txt", "w+");
FILE *vel=fopen("velocity.txt", "w+");
FILE *acc=fopen("acceleration.txt", "w+");
FILE *ene=fopen("energy.txt", "w+");
/*file naming end*/

initialize(pos, vel, x, v);    
force(x, acc, a, &U, time);

for (time=timestep; time<=1000*timestep; time=time+timestep)
{	
	position(pos, x, v, a, time);
	force(x, acc, a, &U, time);
	velocity(vel, v, a, &K);
	T=K+U;
   	printf("Time: %lf Potential: %lf Kinetic: %lf Tot: %lf\n",time,U, K, T);
	fprintf(ene, "\nTime: %lf Potential: %lf Kinetic: %lf Tot: %lf", time, U, K, T);
}

fclose(pos);
fclose(vel);
fclose(acc);
fclose(ene);
printf("\ndone\n");

}
