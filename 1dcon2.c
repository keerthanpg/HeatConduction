//**********************************HEAT CONDUCTION OF A 1 DIMENTIONAL ROD COMPOSED OF PARTICLES WITH SAME MASSES******************************************//
//*********************************************************************************************************************************************************//
//  1)Initialize the positions of molecule equispaced along the 1D ROD                                                                                     //     
//  2)Initialize the velocities of the molecules as normal ranom Numbers                                                                                   //
//  3)Evolve the system using the formula  x(t+dt)=x(t)+v(t)*dt                                                                                            // 
//  4)The molecules have a Hard Potential with a radius "r_cut"                                                                                            //     
//  5)When The molecules have a distance less than 2*r_cut the velocities gets exchanged                                                                   //
//  6)When the molecules cross the wall they are reflected back in with a velocity from the rayleigh distribution                                          // 
//*********************************************************************************************************************************************************//
//*********************************************************************************************************************************************************// 


#include<stdio.h>
#include<math.h>
#include<stdlib.h>
int     no_moles=1281;              // Total number of particles
double  Length=1281;                // Length of the 1-D box 
double  Tl=8;                      // Temperature at the left end
double  Ti=1;					             // Temperature of the rod
double  Tr=2;                      // Temperature at the right end   
double  m1=1,m2=1.22;            // Mass of each atom  
double  kb=1;                      // Boltzman constant
double  time=0.000000;             //
double  r_cut=0.00001;             // Hard Shell radius of an molecule 
int     ind=0;                     // Index for the random velocity array                 
double  ran_vel[10000];            // Random Velocity array 
#define PI 3.14159265359           // Constant Pi 

//Generates Normal Random Numbers with Mean "mu" and Standard Deviation "sigma"
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

//Initializes the position of each molecule; equally spaced
void initialize_pos(FILE* pos, double x[no_moles])
{
    int i;
    fprintf(pos,"time:0\n");
    for(i=0;i<no_moles;i++)
    {
    	x[i]= Length/(no_moles*2) + (Length/no_moles)*i;
            fprintf(pos,"\nmol:%d %lf",i,x[i]);
    }
    fprintf(pos,"\n");
    return;
}

//Initializes the system of molecules
void initialize(FILE* pos, FILE* vel,double x[no_moles], double v[no_moles],double* K) 
{
	 int i;
     //*K stores the value of the Kinetic Energy of the system of molecules in a time instance
  	*K=0.0;

  	fprintf(vel, "\ntime=0\n");

  	initialize_pos(pos, x);
  	for(i=0;i<no_moles;i++)
  	{
        //Intializing the velocities as a normal random number 
        //Initializing the system at Temperature Ti
        if(i%2==0)
          v[i]=randn(0,sqrt(kb*Ti/m1));
        else v[i]=randn(0,sqrt(kb*Ti/m2));                        
  		
        fprintf(vel, "\nMole%d %lf", i, v[i]);
       if(i%2==0) 
  		  *K +=0.5*m1*v[i]*v[i];
      else *K+=0.5*m2*v[i]*v[i];
	 }

  	fprintf(vel, "\n");

  	return;
}


//Evaluating and storing the density function n(x,t)=<sum(delta(x-xl))>
void density(FILE* dense,double x[no_moles],double v[no_moles],double den[no_moles],int cnt,double dt_min)
{
   int i,i_temp,k=0;
   double dt,dt_sum=0;;
  
  if(cnt==0)
    {
      for(i=0;i<no_moles;i++)
        den[i]=0;
      return;
    }

  //************WRITING DENSITY*****************//
    if(cnt==20000000)
      fprintf(dense,"\nTime=%lf\n",time);
  //********************************************//

     
   for(i=0;i<no_moles;i++)
     den[i]*=time;

    for(i=0;i<no_moles;i++)    
    { 
      i_temp=(int)(x[i]*Length/no_moles);
      while((dt_sum<dt_min)&&(i_temp+k<no_moles)&&(i_temp+k>=0))
      {
        if(k==0)
        {
          if(v[i]>0)
            dt=(Length*(i_temp+1)/no_moles-x[i])/v[i];
          else dt=(Length*(i_temp-1)/no_moles-x[i])/v[i];
          dt_sum+=dt;
        }
        else 
        {
          if(v[i]>0)
            dt=Length/(no_moles*v[i]);
          else dt=-Length/(no_moles*v[i]);
          dt_sum+=dt;
        }

        if(dt_sum<dt_min)
          den[i_temp+k]+=dt;
        else den[i_temp+k]+=(dt_min + dt - dt_sum);
        if(v[i]>0)
          k++;
        else k--;
      }
      dt_sum=0;
      k=0;
    }

    for(i=0;i<no_moles;i++)
      {
        if(cnt==0)
          den[i]=0;
        else den[i]/=(time+dt_min);
      }
  //***********PRINTING ENERGY DENSITY***************//
    if(cnt==20000000)
    {
      for(i=0;i<no_moles;i++)
        fprintf(dense,"Mole%d %lf\n",i,den[i]);
    }
  //*************************************************//
  return;
}

//Evaluating and storing the energy density function e(x,t)=<sum(0.5*ml*ul^2delta(x-xl))>
void energy_density(FILE* ene,double x[no_moles],double v[no_moles],double e_den[no_moles],int cnt,double dt_min)
{
  int i,i_temp,k=0;
  double dt,dt_sum=0;
  
   if(cnt==0)
    {
      for(i=0;i<no_moles;i++)
        e_den[i]=0;
      return;
    }

  //************WRITING ENERGY DENSITY************//
      if(cnt==20000000)
        fprintf(ene,"\nTime=%lf\n",time);
  //********************************************//

  
  for(i=0;i<no_moles;i++)
    e_den[i]*=time;

  for(i=0;i<no_moles;i++)
  {
    i_temp=(int)(x[i]*Length/no_moles);
    while((dt_sum<dt_min)&&(i_temp+k<no_moles)&&(i_temp+k>=0))
    {
        if(k==0)
        {
          if(v[i]>0)
            dt=(Length*(i_temp+1)/no_moles-x[i])/v[i];
          else dt=(Length*(i_temp-1)/no_moles-x[i])/v[i];
          dt_sum+=dt;
        }
        else 
        {
          if(v[i]>0)
            dt=Length/(no_moles*v[i]);
          else dt=-Length/(no_moles*v[i]);
          dt_sum+=dt;
        }

        if(dt_sum<dt_min)
        {
         if(i%2==0)
            e_den[i_temp+k]+=0.5*m1*v[i]*v[i]*dt;
         if(i%2==1)
            e_den[i_temp+k]+=0.5*m2*v[i]*v[i]*dt; 
        }
        else
        {
          if(i%2==0)
            e_den[i_temp+k]+=0.5*m1*v[i]*v[i]*(dt_min+ dt -dt_sum );
          if(i%2==1)
            e_den[i_temp+k]+=0.5*m2*v[i]*v[i]*(dt_min+ dt -dt_sum );      
        }
        if(v[i]>0)
          k++;
        else k--;
    }
      dt_sum=0;
      k=0;
  }
 

  for(i=0;i<no_moles;i++)
    e_den[i]/=(time+dt_min);
  //***********PRINTING ENERGY DENSITY***************//
    if(cnt==20000000)
    {
      for(i=0;i<no_moles;i++)
        fprintf(ene,"Mole%d %lf\n",i,e_den[i]);
    }
  //*************************************************//
    
  return;
}

//Evaluating and storing the heat current function j(x,t)=<sum(0.5*ml*ul^3delta(x-xl))>
void heat_current(FILE* heat,double x[no_moles],double v[no_moles],double h_cur[no_moles],int cnt,double dt_min)
{
  int i,i_temp,k=0;
  double dt,dt_sum=0;
  
  if(cnt==0)
    {
      for(i=0;i<no_moles;i++)
        h_cur[i]=0;
      return;
    }

  //************WRITING HEAT CURRENT************//
    if(cnt==20000000)
      fprintf(heat,"\nTime=%lf\n",time);
  //********************************************//
    
  for(i=0;i<no_moles;i++)
    h_cur[i]*=time;

  for(i=0;i<no_moles;i++)
  {
    i_temp=(int)(x[i]*Length/no_moles);
    while((dt_sum<dt_min)&&(i_temp+k<no_moles)&&(i_temp+k>=0))
    {
        if(k==0)
        {
          if(v[i]>0)
            dt=(Length*(i_temp+1)/no_moles-x[i])/v[i];
          else dt=(Length*(i_temp-1)/no_moles-x[i])/v[i];
          dt_sum+=dt;
        }
        else 
        {
          if(v[i]>0)
            dt=Length/(no_moles*v[i]);
          else dt=-Length/(no_moles*v[i]);
          dt_sum+=dt;
        }

        if(dt_sum<dt_min)
        {
         if(i%2==0)
            h_cur[i_temp+k]+=0.5*m1*v[i]*v[i]*v[i]*dt;
         if(i%2==1)
            h_cur[i_temp+k]+=0.5*m2*v[i]*v[i]*v[i]*dt; 
        }
        else
        {
          if(i%2==0)
            h_cur[i_temp+k]+=0.5*m1*v[i]*v[i]*v[i]*(dt_min+ dt -dt_sum );
          if(i%2==1)
            h_cur[i_temp+k]+=0.5*m2*v[i]*v[i]*v[i]*(dt_min+ dt -dt_sum );      
        }
        if(v[i]>0)
          k++;
        else k--;
    }
      dt_sum=0;
      k=0;
  }

  for(i=0;i<no_moles;i++)
    h_cur[i]/=(time+dt_min);
  //***********PRINTING HEAT CURRENT*****************//
    if(cnt==20000000)
    {
      for(i=0;i<no_moles;i++)
      fprintf(heat,"Mole%d %lf\n",i,h_cur[i]);
    }
  //*************************************************//
    
  return;
}

//Evaluates and stores the Temperature function T(x,t)=2*e(x,t)/n(x,t) 
void Temperature(FILE* temp,double e_den[no_moles],double den[no_moles],double T[no_moles],int cnt)
{
    int i;
    
    if(cnt==0)
    {
      for(i=0;i<no_moles;i++)
        T[i]=0;
      return;
    }

    //************WRITING TEMPERATURE************//
      if(cnt==20000000)
      fprintf(temp,"\nTime=%lf\n",time);
    //********************************************//
    for(i=0;i<no_moles;i++)
      T[i]=2*e_den[i]/den[i];
    

    //***********PRINTING TEMPERATURE*****************//
        if(cnt==20000000)
        {
          for(i=0;i<no_moles;i++)
             fprintf(temp,"Mole%d %lf\n",i,T[i]);
        }
    //*************************************************//
     
     return;
}


//Evolves the position and velocity of the molecule
void position(FILE* pos,FILE* vel,FILE* dense,FILE* ene,FILE* heat,FILE* temp,double x[no_moles], double v[no_moles],double den[no_moles],double e_den[no_moles],double h_cur[no_moles],double T[no_moles] ,double* K,int cnt) 
{
  double dt_min=10,del_x[no_moles+1],rel_v[no_moles+1],foo;
  int i,i_min=-1;
  //dt_min is the least time in which a pair of molecules has a distance of 2*r_cut between them
  //i_min is the lowest index of the particle which comes close to 2*r_cut with another
  //del_x required displacement
  //rel_v(i) is the relative velocity of the ith particle w.r.t i+1th particle ; 0th particle being the left wall and no_mole+1th particle being the right wall
  //foo is a temporary variable
  *K=0;
  

//**************Calculating dt_min begins**************//
  del_x[0]=x[0]-r_cut;
  if(del_x[0]<0)
  {//Precautionary check
    printf("error 0");
    scanf("%lf",&foo);
  }
  rel_v[0]=-v[0];
  if(del_x[0]/rel_v[0]>0)
    dt_min=del_x[0]/rel_v[0];

  for(i=0;i<no_moles-1;i++)
  {
    del_x[i+1]=x[i+1]-x[i]-2*r_cut;
    rel_v[i+1]=-v[i+1]+v[i];
    if( (del_x[i+1]/rel_v[i+1]>0) && (del_x[i+1]/rel_v[i+1]<dt_min))
      {
        dt_min=del_x[i+1]/rel_v[i+1];
        i_min=i;
      }
  }

  del_x[no_moles]=Length-x[no_moles-1]-r_cut;
  if(del_x[no_moles]<0)
  {//Precautionary check
    printf("error Length");
    scanf("%lf",&foo);
  }
  rel_v[no_moles]=v[no_moles-1];
  if( (del_x[no_moles]/rel_v[no_moles]>0) && (del_x[no_moles]/rel_v[no_moles]<dt_min) )
    {
      dt_min=del_x[no_moles]/rel_v[no_moles];
      i_min=no_moles-1;
    }
  //***********Calculating dt_min ends*********//  

    density(dense,x,v,den,cnt,dt_min);
    energy_density(ene,x,v,e_den,cnt,dt_min);
    heat_current(heat,x,v,h_cur,cnt,dt_min);
    Temperature(temp,e_den,den,T,cnt);

    //Precautionary check
    if(dt_min<0)
        printf("dt_min<0");

//Calculating new position
    for(i=0;i<no_moles;i++)
        x[i]+=v[i]*(dt_min);
  
//*****************Calculating new velocity begins  ****************//
    if(i_min==-1)
    {
        v[0]=ran_vel[ind]*sqrt((kb*Tl)/m1); 
        ind=(ind+1)%10000;
        x[0]=r_cut*1.000001;
    }
    else if(i_min==no_moles-1)
    {
        if((no_moles-1)%2==0)      
          v[no_moles-1]=-ran_vel[ind]*sqrt((kb*Tr)/m1); 
        else v[no_moles-1]=-ran_vel[ind]*sqrt((kb*Tr)/m2);
        ind=(ind+1)%10000;
        x[no_moles-1]=Length-r_cut*1.000001;
    }
    else
    {
        foo=v[i_min+1];
        if(i_min%2==0)
        {//Elastic Collision

          v[i_min+1]=(2*m1*v[i_min]-(m1-m2)*foo)/(m1+m2);
          v[i_min]=((m1-m2)*v[i_min]+2*m2*foo)/(m1+m2);
        } 
        else
        {
          v[i_min+1]=(2*m2*v[i_min]-(m2-m1)*foo)/(m1+m2);
          v[i_min]=((m2-m1)*v[i_min]+2*m1*foo)/(m1+m2);
        } 
          x[i_min]-=0.000001*r_cut;
          x[i_min+1]+=0.00001*r_cut;
        
    }
    for(i=0;i<no_moles-1;i++)
    {
        if(x[i+1]-x[i]<2*r_cut-0.000001)
         {//Precautionary check
            printf("dist=%lf",x[i+1]-x[i]);
            scanf("%lf",&foo);
        }
    }

//************Calculating new velocity ends****************//

    //Calculating Kinetic Energy
    for(i=0;i<no_moles;i++)
    {  
      if(i%2==0)
        *K+=0.5*m1*v[i]*v[i];
      else *K+=0.5*m2*v[i]*v[i];
    }
    //Updating time
    time+=dt_min;

    for(i=0;i<no_moles-1;i++)
    {//Precautionary Check
        if(x[i]>x[i+1])
        {
            printf("%d:%lf crossed %d:%lf",i,x[i],i+1,x[i+1]);
            scanf("%d",&i);
        }
    }
  return;
}





int main()
{
  	double x[no_moles],v[no_moles],K,T[no_moles],den[no_moles],e_den[no_moles],h_cur[no_moles];
  	int i,dum;

    /*file naming*/
  	FILE *pos=fopen("position.txt", "w");                //writes the position of all the molecules in different time instances 
  	FILE *vel=fopen("velocity.txt", "w");                //writes the velocity of all the molecules in different time instances
  	FILE *dense=fopen("density.txt","w");
    FILE *ene=fopen("energy_density.txt", "w");          //writes the kinetic energy of all the molecules in different time instances
  	FILE *heat=fopen("heat_current.txt","w");
    FILE* ran=fopen("rayleigh.dat","r");                 //reads random rayleigh generated numbers 
    FILE *temp=fopen("Temperature.txt","w");             //writes steady state temperature(continuum) after different time instances    
  	//FILE *tra=fopen("trajectory.lammpstrj","w+");      //writes the position and velocities of all particles at different time instances to a lammpstrj file
  	/*file naming end*/

    // Reading Random Rayleigh Numbers
  	for(i=0;i<10000;i++)
      fscanf(ran,"%lf",&ran_vel[i]);
        
	i=0;
    //**********Initialization**********//
  	initialize(pos,vel,x,v,&K);    
  	density(dense,x,v,den,0,0);
    energy_density(ene,x,v,e_den,0,0);
    heat_current(heat,x,v,h_cur,0,0);
    Temperature(temp,e_den,den,T,0);
	  printf("Time: %lf  Kinetic:%lf\n",time,K);
    //***********************************//
  	
    for (i=1; i<=20000000;i++)
  	{	
  		position(pos,vel,dense,ene,heat,temp,x,v,den,e_den,h_cur,T,&K,i);
      
 	
    //********************************************************LAMMPSTRJ*******FILE******************************************************************************// 	
 	  /*fprintf(tra,"ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS p p p\n0.000000 %lf\n0.000000 %lf\n0.000000 %lf\n",j,totalParticles,unitLength,unitLength,unitLength);
  		fprintf(tra,"ITEM: ATOMS id type xu yu zu vx vy vz\n");
  		for(k=0;k<no_moles;k++)
  		fprintf(tra,"%d 0 %lf %lf %lf %lf %lf %lf\n",k,x[k][0],x[k][1],x[k][2],v[k][0],v[k][1],v[k][2]);                                                           */                
  	//*************************************************************************************************************************************************************//
    	
  		printf("%d:Time: %lf Kinetic:%lf\n",i,time,K);
  	}

  	//fclose(pos);
  	fclose(vel);
  	fclose(temp);
  	fclose(dense);
    fclose(ene);
    fclose(heat);
  	printf("\ndone\n");
  	return 0;
}
























