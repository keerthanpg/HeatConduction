#include<stdio.h>
#include<math.h>
#include<stdlib.h>
int     no_moles=100;                        // Total number of particles
int 	collisionL=0,collisionR=0;
//int 	collision1=0;
int counterL=0,counterR=0;
double  Length=100;                // Length of the 1-D box 
double  Tl=1.5;                   // Temperature at the left end
double  Ti=1;						//Temperature of the rod
double  Tr=0.5;                  // Temperature at the right end   
double  m=1;                       // Mass of each atom  
double  kb=1;                       // Boltzman constant
double  dt=0.0001;                   // Value of the time step
double r_cut=0.01;
int     ind=0;
double  ran_vel[10000000];
#define PI 3.14159265359
#define COR 1 

//double Ns=100;                     // Number of time steps   

double randn(double mu,double sigma )
{
	double v1,v2,s;

	do 
	{
	    v1 = 2.0 * ((double) rand()/RAND_MAX) - 1;
	    v2 = 2.0 * ((double) rand()/RAND_MAX) - 1;

	    s = v1*v1 + v2*v2;
	} while ( s >= 1.000000 );

	if (s == 0.0)
		return mu+ 0.0;
	else
		return mu+ sigma*(v1*sqrt(-2.0 * log(s) / s));
}

void initialize_pos(FILE* pos, double x[no_moles])
{
    int i;
    //fprintf(pos,"time:0\n");
    for(i=0;i<no_moles;i++)
    {
    	x[i]= Length/(no_moles*2) + (Length/no_moles)*i;
       // fprintf(pos,"\nmol:%d %lf",i,x[i]);
    }
    //fprintf(pos,"\n");
    return;
}

void initialize(FILE* pos, FILE* vel,double x[no_moles], double v[no_moles],double* K,double* J) 
{
	 int i;
  	*K=0.0;
  	*J=0.0;

  	//fprintf(vel, "\ntime=0\n");

  	initialize_pos(pos, x);
  	for(i=0;i<no_moles;i++)
  	{
  		v[i]=randn(0,sqrt(kb*Ti/m));//ran_vel[ind]*sqrt((kb*Ti)/m); ind=(ind+1)%1000000;
  	//	fprintf(vel, "\nMole%d %lf", i, v[i]);
  		
  		*K +=0.5*m*v[i]*v[i];
  		*J +=0.5*m*v[i]*v[i]*v[i];
	 }

  	//fprintf(vel, "\n");

  	return;
}

void correction(double x[no_moles],double v[no_moles],int low_i)
{
    double ti,vdiff,xdiff,foo;
    vdiff=v[low_i+1]-v[low_i];
    xdiff=x[low_i+1]-x[low_i];
    if(vdiff>0)
        ti=(xdiff-vdiff*dt-2*r_cut)/vdiff;
    else ti=-(xdiff-vdiff*dt-2*r_cut)/vdiff;

    x[low_i]+=(v[low_i+1]-v[low_i])*(dt-ti);
    x[low_i+1]-=(v[low_i+1]-v[low_i])*(dt-ti);
    foo=v[low_i+1];
    v[low_i+1]=v[low_i];
    v[low_i]=foo;
    return;        
}

void position(FILE* pos,FILE* vel,double x[no_moles], double v[no_moles],double* K,double* J,int count) 
{
	int i;
	double ti,foo;
	*K=0;
	/*if(count%10000==0)
    {
        fprintf(pos, "\nTime=%lf\n",count*dt);
        fprintf(vel, "\nTime=%lf\n",count*dt);
    }*/
	
    for(i=0;i<no_moles;i++)
       x[i]+=v[i]*dt ;
      
  	   	
    // Temperature Controlled Boundary
  	if(x[no_moles-1]>=Length)
	 {
        ti=(Length-x[no_moles-1])/v[no_moles-1]+dt;
        foo=v[no_moles-1];
        //do{
            v[no_moles-1]=/*randn(0,sqrt(kb*Tr/m));*/-ran_vel[ind]*sqrt((kb*Tr)/m); ind=(ind+1)%10000;
        //}while(v[no_moles-1]>=0);
        if(v[no_moles-1]*v[no_moles-1]<foo*foo)
			counterR++;
        //fprintf(vel,"%lf\n",v[no_moles-1]);
        //  v[no_moles-1]=-1;
     	x[no_moles-1]=Length + v[no_moles-1]*(dt - ti);
      	collisionR++ ;
  	}

  	if(x[0]<=0.0)
  	{
    	ti=-x[0]/v[0]+dt;
    	foo=v[0];
    	//do{
           v[0]=/*randn(0,sqrt(kb*Tl/m));*/ran_vel[ind]*sqrt((kb*Tl)/m); ind=(ind+1)%10000;
       	//}while(v[0]<=0);
          if(v[0]*v[0]<foo*foo)
			counterL++;
           //v[0]=1;
        //fprintf(pos,"%lf\n",v[0]);
      	x[0]= v[0]*(dt - ti);
    	collisionL++;
  	}
  	
  	for(i=0;i<no_moles-1;i++)
  	{
		if((x[i+1]-x[i]<2*r_cut)&&(x[i]<x[i+1]))
		{
			foo=v[i+1];
			v[i+1]=v[i];
			v[i]=foo;
			//collision1++;
		}
        if(x[i+1]<x[i])
        {
            correction(x,v,i);
        }
	}

    for(i=0;i<no_moles-1;i++)
    {
        if(x[i]>x[i+1])
        {
            printf("%d crossed %d",i,i+1);
            scanf("%d",&i);
        }
    }
	
    *J *=count;
	
    for(i=0;i<no_moles;i++)
	{
		*J +=0.5*m*v[i]*v[i]*v[i];
		*K +=0.5*m*v[i]*v[i];
	}
	*J/=(count+1);	
    /*if(count%10000==0)
    {
        for(i=0;i<no_moles;i++)
        {
            fprintf(pos, "\nMole:%d %lf ", i, x[i]);
            fprintf(vel, "\nMole:%d %lf ", i, v[i]);
        }
    }*/
 	
 	printf(" \ncollisions wallL=(%d,%d) wallR=(%d,%d) ",collisionL,counterL,collisionR,counterR);
	return;
}


void Temperature(FILE* temp,double x[no_moles],double v[no_moles],double T[no_moles],double count[no_moles],int cnt)
{
    int i,k;
    if(cnt%1000000==0)
      fprintf(temp,"\nTime=%lf\n",cnt*dt);
    if(cnt==0)
    {
        for(i=0;i<no_moles;i++)
        T[i]=0;
        count[i]=0;
    }
    else
    {
        for(i=0;i<no_moles;i++)
            T[i]*=count[i];
    }    

    for(i=0;i<no_moles;i++)    
    {
        for(k=-4;k<4;k++)
        {
            if(i+k>=0 && i+k<no_moles)
            {
                if( ( x[i+k] <= Length*(i+1)/no_moles ) &&( x[i+k] > Length*i/no_moles ) )
                { 
                    T[i]+=m*v[i+k]*v[i+k]/kb;
                    count[i]=count[i]+1;
                }
            }
        }
    }
    for(i=0;i<no_moles;i++)
    {
        if(count[i]==0)
            T[i]=0;
        else T[i]/=(count[i]);
    }
    if(cnt%1000000==0)
     {
        for(i=0;i<no_moles;i++)
             fprintf(temp,"Mole%d %lf\n",i,T[i]);
     }  
     return;
}


int main()
{
  	double x[no_moles],v[no_moles],K,J,T[no_moles],Count[2*no_moles];
  	int i,dum;

    /*file naming*/
  	FILE *pos=fopen("position.txt", "w");
  	FILE *vel=fopen("velocity.txt", "w");
  	FILE *ene=fopen("energy.txt", "w");
  	FILE* ran=fopen("rayleigh.dat","r");
    FILE *temp=fopen("Temperature.txt","w"); 
  	//FILE *tra=fopen("trajectory.lammpstrj","w+");
  	/*file naming end*/
  	for(i=0;i<10000;i++)
      fscanf(ran,"%lf",&ran_vel[i]);
        
	   i=0;

  	initialize(pos,vel,x,v,&K,&J);    
  	Temperature(temp,x,v,T,Count,0);
  	
	printf("Time: %lf  Kinetic:%lf J:%1.9lf\n",i*dt,K,J);

  	for (i=1; i<=30000000;i++)
  	{	
  		position(pos,vel,x,v,&K,&J,i);
        Temperature(temp,x,v,T,Count,i);
 	/* 	
 		fprintf(tra,"ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS p p p\n0.000000 %lf\n0.000000 %lf\n0.000000 %lf\n",j,totalParticles,unitLength,unitLength,unitLength);
  		fprintf(tra,"ITEM: ATOMS id type xu yu zu vx vy vz\n");
  		for(k=0;k<no_moles;k++)
  		fprintf(tra,"%d 0 %lf %lf %lf %lf %lf %lf\n",k,x[k][0],x[k][1],x[k][2],v[k][0],v[k][1],v[k][2]);
  	*/
    	
  		printf("Time: %lf Kinetic:%lf J:%1.8lf\n",i*dt,K,J);
  		if(i%1000==0)
        {
          fprintf(ene, "Time: %lf Kinetic: %lf J: %lf\n",i*dt, K,J);
        }
       
  	}

  	//fclose(pos);
  	fclose(vel);
  	fclose(temp);
  	fclose(ene);
  	printf("\ndone\n");
  	return 0;
}
























