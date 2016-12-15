#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define k 1.0
#define m 1.0
#define timestep 0.01
#define sigma 1.0
#define epsilon 1.0
#define T 1.0
#define lattice 20.0
int no_moles=1000;

double randn(double mu, double sig)
{
	
	double v1,v2,s;
	
	
	 
	do 
	{
	    v1 = 2.0 * ((double) rand()/RAND_MAX) - 1;
	    v2 = 2.0 * ((double) rand()/RAND_MAX) - 1;
	    s = v1*v1 + v2*v2;
	} while ( s >= 1.000000 );
	if (s == 0.0)
		return( mu+ 0.0);
	else
	return( mu+ sig*(v1*sqrt(-2.0 * log(s) / s)));
}	
void initialize_pos(double x[][3])
{
int p=1, q=1, r=1;
int i=0;
double x1, x2, x3;
double v_mean=0;
double v_stddev=sqrt(k*T/m);


for(p=1;p<=10;p++)
{
	for(q=1;q<=10;q++)
	{
		for(r=1;r<=10;r++)
		{
			x1=1.8*p;    
			x2=1.8*q;
			x3=1.8*r;           

			x[i][0]=x1; x[i][1]=x2; x[i][2]=x3;
			i++;
			if(i>=no_moles)
			{
				return;
			}
		}
	}
}
}

void initialise_velocity(double v[][3],double a[][3][2])
{    
 
      int i; 
      double sig;
      sig=sqrt((k*T)/m); 
     	
        
        
      for (i=0;i<no_moles;i++)
      { 
          v[i][0]=randn(0,1);
          v[i][1]=randn(0,1);
          v[i][2]=randn(0,1);
           a[i][0][1]=0;
           a[i][1][1]=0;
           a[i][2][1]=0;
           a[i][0][0]=0;
           a[i][1][0]=0;
           a[i][2][0]=0;
     }
 } 
  
void position(double x[no_moles][3], double v[no_moles][3],double a[no_moles][3][2])
{     int i,j;
  
           for( i=0;i<no_moles;i++)
           { 
               x[i][0] += v[i][0]*timestep+0.5*(a[i][0][1])*timestep*timestep;
               x[i][1] += v[i][1]*timestep+0.5*(a[i][1][1])*timestep*timestep;
               x[i][2] += v[i][2]*timestep+0.5*(a[i][2][1])*timestep*timestep;
               
               for(j=0; j <3; j++)
               {
				   if (x[i][j] > lattice)
					 x[i][j] -= lattice;
				   else if(x[i][j] < 0.0)
					 x[i][j] += lattice;
			   }
           } 
 }   
 
 void velocity(double v[no_moles][3], double a[no_moles][3][2] )
{
	 int i;
    

      for(i=0;i<no_moles;i++)
      {     v[i][0] += 0.5 * ( a[i][0][0] + a[i][0][1] ) * timestep;
            v[i][1] += 0.5 * ( a[i][1][0] + a[i][1][1] ) * timestep;
            v[i][2] += 0.5 * ( a[i][2][0] + a[i][2][1] ) * timestep;
       }
}
 
double potential(double x[no_moles][3])
   { double U=0; 
    
    double r;int a;
    int i,j,n;
     double x_diff[3];
    
    
    for(i=0; i<no_moles; i++)
      { for (j=i+1; j<no_moles; j++)
        {    	
			for(n=0;n<3;n++)
		    {
				x_diff[n]=x[i][n]-x[j][n];
				if (x_diff[n]>lattice/2.0)
					x_diff[n]-=lattice;
				else if (x_diff[n]<-lattice/2.0)
					x_diff[n]+=lattice;
				         			
			}
			r=sqrt((x_diff[0]*x_diff[0])+(x_diff[1]*x_diff[1])+(x_diff[2]*x_diff[2]));  
            if(r<2.5)
			{
				U +=  4.0*epsilon*( pow((sigma/r),12.0)- pow((sigma/r),6.0)); 
			}
        }
     }    
     return U;
     
 }  
 
void acceleration( double x[no_moles][3] ,double a[no_moles][3][2] )
   {  double r;
      double x_diff[3];
      int i,j,n,temp;
    
       for( i=0;i<no_moles;i++)
        {
			a[i][0][0]=a[i][0][1];
			a[i][1][0]=a[i][1][1];
			a[i][2][0]=a[i][2][1];
			a[i][0][1]=0.0;
			a[i][1][1]=0.0;
			a[i][2][1]=0.0;
       
        } 
           
       for(i=0; i<no_moles; i++)
         {
         for(j=0; j<no_moles;j++)
                {
                 	if(i!=j)
                     { 
						 for(n=0;n<3;n++)
		    	         {
							x_diff[n]=x[i][n]-x[j][n];
							if (x_diff[n]>lattice/2.0)
								x_diff[n]-=lattice;
							else if (x_diff[n]<-lattice/2.0)
								x_diff[n]+=lattice;
				         			
					      }
						r=sqrt((x_diff[0]*x_diff[0])+(x_diff[1]*x_diff[1])+(x_diff[2]*x_diff[2]));  
						
						if(r < 0.1) 
						{
							printf("%d %d",i,j);
							printf("Enter a number to continue");
							scanf("%d",&temp);
						}
                        if  (r<2.5)
 			            {
						    a[i][0][1] -= 4.0*epsilon*(x_diff[0]/r)*( 6.0*pow((sigma/r),7)-12.0* pow((sigma/r),13));
						    a[i][1][1] -= 4.0*epsilon*(x_diff[1]/r)*( 6.0*pow((sigma/r),7)-12.0* pow((sigma/r),13));
                            a[i][2][1] -= 4.0*epsilon*(x_diff[2]/r)*( 6.0*pow((sigma/r),7)-12.0* pow((sigma/r),13));
						}
                    }
             
           }       
      }

        
}   


  double kinetic(double v[no_moles][3])
  {    double V=0; 
       int i;
      
       for(i=0;i<no_moles;i++)
         {
			V+=0.5*m*(v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2]); 
         }
         
       
        return V;
       
 }  
 
  int main()
  {    
       double v[no_moles][3], TE; double time;
     double x[no_moles][3],a[no_moles][3][2],U,K;
   
       initialize_pos(x);
	
       initialise_velocity(v,a);
       acceleration(x,a);
	   
	   for ( time=0; time<=10000*timestep; time=time+timestep)
		{ 
		   position(x,v,a);
		   acceleration(x,a);
		   velocity(v,a);
		  
       
			U=potential(x);
			K=kinetic(v);
			TE=U+K;
			printf("potential: %lf\tkinetic: %lf\ttotal:%lf\n ",U,K,TE);
       }
      
      return 0;
  }      
    
         
  
 


        
      
                 
      
           
           
      
       
                
