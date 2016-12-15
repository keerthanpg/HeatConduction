/*Classical 4th order Runge Kutta Method*/
#include<iostream.h>

void main()
{
double x[1000000], p[1000000], n[1000000], k[4];
x[0]=0.0;p[0]=0.0;n[0]=0.0;
double h=0.01;
int i;
for(i=1;i<1000000;i++)
{
k[0]=p[i-1];
k[1]=p[i-1]+h/2*k[0];
k[2]=p[i-1]+h/2*k[1];
k[3]=p[i-1]+h*k[2];
x[i]=x[i-1]+h/6*(k[0]+2*k[1]+2*k[2]+k[3]);

k[0]=p[i-1]*p[i-1]-1;
k[1]=(p[i-1]+h/2*k[0])*(p[i-1]+h/2*k[0])-1;
k[2]=(p[i-1]+h/2*k[1])*(p[i-1]+h/2*k[1])-1;
k[3]=(p[i-1]+h*k[2])*(p[i-1]+h*k[2])-1;
n[i]=n[i-1]+h/6*(k[0]+2*k[1]+2*k[2]+k[3]);

k[0]=-x[n-1]-n[i-1]*p[i-1];
k[1]=-(x[n-1]+h/2*k[0])-(n[i-1]+h/2*k[0])*(p[i-1]+h/2*k[0]);
k[2]=-(x[n-1]+h/2*k[1])-(n[i-1]+h/2*k[1])*(p[i-1]+h/2*k[1]);
k[3]=-(x[n-1]+h*k[2])-(n[i-1]+h*k[2])*(p[i-1]+h*k[2]);
p[i]=p[i-1]+h/6*(k[0]+2*k[1]+2*k[2]+k[3]);
}

print("x=%lf, p=%lf, n=%lf\n", x[i], p[i], n[i]);
}

