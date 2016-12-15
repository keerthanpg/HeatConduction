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
double timestep=0.01, time=0.0;
int boundingl=0, boundingr=0;
 
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
 
   U1 = (double) rand () / (double)RAND_MAX;
  U2 = (double) rand () / (double)RAND_MAX;
    
  X1 = sqrt(-2*log(U1))*cos(2*3.14*U2);
  X2= sqrt(-2*log(U1))*sin(2*3.14*U2);
 
  call = !call;
 
  return (mu + sigma * (double) X1);
}
	
void initialize_pos(FILE* pos, double x[])
{

int i=0;
double x1;

fprintf(pos, "time=0\n");

for(i=0;i<=no_moles;i++)
{
	x[i]=L/(2.0*(double)no_moles)+ L/((double)no_moles+1)*(double)i;
	fprintf(pos, "\n%d %lf ", i, x[i]);
				
}

fprintf(pos, "\n");
return;

}

void initialize_pos_random(FILE* pos, double x[no_moles])
{
	int i;
	for(i=0;i<no_moles; i++)
	{
		fprintf(pos, "\n%d ", i+1);
		x[i]=(double)rand()/(double)RAND_MAX*L;
		fprintf(pos, "%lf ", x[i]);
		fprintf(pos, "\n");
	}

return;	
}

void initialize(FILE* pos, FILE* vel, double x[], double v[], double *K) 
{
int i=0, j;
double v_mean=0.0;
double v_stddev=sqrt(k*Tl/m);
double velocity;

fprintf(vel, "\ntime=0\n");
*K=0.0;
printf("v_stddev=%lf", v_stddev);
initialize_pos(pos, x);
FILE* Inivel=fopen("RandVel1_5.txt", "r");
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
double LeftRecursiveX(FILE* RandVel, double* x, double* v, double a, double t)
{
	double bct, t1, t2;//before collision time
	double randvelocity;
	boundingl++;
	*x-=(*v)*(t)+0.5*a*(t)*(t);
	
	printf("\n Bounding at left!!! No. %d\n", boundingl);
	t1=(-(*v)+sqrt((*v)*(*v)-2.0*a*(*x)))/a;
	t2=(-(*v)-sqrt((*v)*(*v)-2.0*a*(*x)))/a;
	if(t1>=0)
	bct=t1;
	else bct=t2;
	if(fscanf(RandVel, "%lf ", &randvelocity)!=EOF)
		{
			*v=randvelocity; 
			if (*v<0)
			*v=fabs(*v);				
		}
	else
		exit(0);
	*x=0.0;
	*x+=(*v)*(t-bct)+0.5*a*(t-bct)*(t-bct);
	if(boundingl==4)
	{
	printf("\n %lf, %lf, %lf\n ", *x, *v, a);
	//exit(0);
	
	}
	
	return(t-bct);
}

double RightRecursiveX(FILE* RandVel, double *x, double *v, double a, double t)
{
	double bct, t1, t2;//before collision time
	double randvelocity;
	boundingr++;
	*x-=(*v)*(t)+0.5*a*(t)*(t);
	printf("\n Bounding at right!!! No. %d\n", boundingr);
	t1=(-(*v)+sqrt((*v)*(*v)+2.0*a*(*x)))/a;
	t2=(-(*v)-sqrt((*v)*(*v)+2.0*a*(*x)))/a;
	printf("%lf \n", a);
	if(t1>=0)
	bct=t1;
	else bct=t2;
	if(fscanf(RandVel, "%lf ", &randvelocity)!=EOF)
		{
			*v=randvelocity; 
			if (*v>0)
			*v=-*v;				
		}
	else
		exit(0);
	*x=L;
	*x+=(*v)*(t-bct)+0.5*a*(t-bct)*(t-bct);
	printf(" %lf", *x);
	//exit(0);
	return(t-bct);
}

void position(FILE*RandVel, FILE* pos, double x[no_moles], double v[no_moles], double a[no_moles][2], double time) 
{
	int i, j=0;
	double vl_mean=0.0;
	double vl_stddev=sqrt(k*Tl/m);
	double vr_mean=0.0;
	double vr_stddev=sqrt(k*Tr/m);
	double randvelocity;
	fprintf(pos, "\nTime=%lf\n", time);
	double t;
	
	for(i=0;i<no_moles;i++)
    {         
		x[i]+=v[i]*timestep+0.5*a[i][1]*timestep*timestep;
		if(x[i]<=0)
		{
			t=timestep;
			do
			{
			t=LeftRecursiveX(RandVel,&x[i], &v[i], a[i][1], t);
			}while(x[i]<=0);
			
		}
		else if(x[i]>=L)
		{
			t=timestep;
			//do
			//{
			t=RightRecursiveX(RandVel, &x[i], &v[i], a[i][1], t);
			//}while(x[i]>=L);		
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
	fprintf(pos, "\n%d %lf ", i, x[i]);
	return;
}

void velocity(FILE* vel, double v[no_moles], double a[no_moles][2], double *K)
{
	int i,k=0;	
	*K = 0.0;
	
	fprintf(vel, "\ntime=%lf\n ", time);
	for(i=0;i<no_moles;i++)
    {   
		v[i]+=0.5*timestep*(a[i][0]+a[i][1]);	
		*K +=0.5*m*v[i]*v[i];		
		a[i][0]=a[i][1];
		fprintf(vel, "\n%d %lf ", i, v[i]);
	}
	
	return;
}
	

void main()
{

double v[no_moles], x[no_moles], U, a[no_moles][2], K, TE, Temp;
int i, j;
for(i=0;i<no_moles;i++)
{
a[i][0]=0;
a[i][1]=0;
}

/*file naming*/
FILE *pos=fopen("position.txt", "w+");
FILE *vel=fopen("velocity.txt", "w+");
FILE *acc=fopen("acceleration.txt", "w+");
FILE *ene=fopen("energy.txt", "w+");
FILE *temp=fopen("temperature.txt", "w+");
FILE* RandVel=fopen("RandVel.txt", "r");
/*file naming end*/


initialize(pos, vel, x, v, &K);    


for (time=timestep; time<=10000000*timestep; time=time+timestep)
{	
	position(RandVel, pos, x, v, a, time);
	
	velocity(vel, v, a, &K);
	TE=K+U; 
	Temp=K/(0.5*(double)no_moles*k);
   	printf("Time: %lf U: %lf K: %lf T: %lf",time,U, K, TE);
   	printf("Temp=%lf K\n", K/(0.5*k*(double)no_moles));
   	fprintf(temp,"\n %lf, %lf", time, K/(0.5*k*(double)no_moles));
	fprintf(ene, "\nTime: %lf Temp: %lf K U: %lf K: %lf T: %lf",  time, Temp, U, K, TE);
	
}

printf("\nboundingl=%d, boundingr=%d", boundingl, boundingr);
fclose(pos);
fclose(vel);
fclose(acc);
fclose(ene);
fclose(temp);
fclose(RandVel);

printf("\ndone\n");

}
