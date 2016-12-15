/*Adaptive Runge Kutta
Runge Kutta Fehlberg method*/
#include<iostream.h>
long double f(int a, long double x, long double p, long double n)
{
if(a==1)
return(p);
else if(a==2)
return(-x-np);
else if (a==3)
return(p*p-1);
}

void main()
{
double x[2], p[2], n[2],k[6];
x[0]=0.0;p[0]=0.0;n[0]=0.0;
double h=0.01;
int i;
for(i=1;i<1000000;i++)
{
for(j=1;j<=3;j++)
{
k[0]=h*f(j, x[1], p[1], n[1]);
k[1]=h*f(j, x[1]+h*
}
}
printf("x=%lf, p=%lf, n=%lf", x[i], p[i], n[i]);

}
