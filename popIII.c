#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"
#ifdef OMP_NUM_THREADS
#include <pthread.h>
#endif

#if defined(POPIII) && defined(GALSF) && defined(METALS)

//ergs per event
#define ENERGY_TYPE2 1.0e51
#define ENERGY_PISN  1.0e52
#define ENERGY_BH    0.0
#define ENERGY_OTHER 0.0

// Bubble properties 
#define HII_DENSITY_ITERATION 1
#define HII_DENSITY_POPIII 0.01
#define I_FRONT_RADIUS 2.0     //kpc/h PHYSICAL
#define TEMP_BUBBLE 2.0e4  //Kelvin
#define E_FRAC_BUBBLE 0.10  //n_e/n_H

// yield per type of SN event
#define YIELD_TYPE2 0.10
#define YIELD_PISN  0.50
#define YIELD_BH    0.0
#define YIELD_OTHER 0.0

#define LOWER 1.0     //lower limit to IMF
#define UPPER 150.0   //upper limit to IMF
#define DM 0.1        //dM for IMF table
//#define N_VALS 1490   //(UPPER-LOWER)/DM
#define N_VALS (int)((UPPER-LOWER)/DM)

static double prob[N_VALS];
static double masses[N_VALS];
double trap(double (*fm)(double, double, double),int, double, double, double, double);  //trapezoidal integration for m-larson imf
double fm(double, double, double); //m-larson function (mass)
double fn(double, double, double); //m-larson function (number)

void InitLegacy(){
  int i;
  double IMF_norm, N,sum=0.0, M1,M2;

  srand48(ThisTask);
  
  double Msun_to_Internal  = All.HubbleParam * SOLAR_MASS / All.UnitMass_in_g;
  //double Upper = UPPER * Msun_to_Internal;
  //double Lower = LOWER * Msun_to_Internal;
  //double dm    = DM * Msun_to_Internal;
  
  //All.PopIII_Mass = POPIII_MASS * Msun_to_Internal;
  All.PopIII_Mass = 0.0;
  if(RestartFlag == 0){
    All.POPIII_SFR = 0.0;
    All.POPIII_SFRD = 0.0;
    All.POPIII_C_SFRD = 0.0;
    All.POPIII_TIME = 0.0;
    All.POPIII_SUM = 0.0;
    All.POPIII_flag = 0;
    //All.POPII_flag = 0;
}
  // if power-law IMF mcut=0.0 ==> exp(term)=1.0 in fn() and fm()
  #ifndef POPIII_MLARSON
  double POPIII_MCUT=0.0;
  #endif

  //normalize the IMF 
  IMF_norm = POPIII_MASS / trap(fm,1000,LOWER,UPPER,POPIII_ALPHA,POPIII_MCUT);
  N = IMF_norm * trap(fn,1000,LOWER,UPPER,POPIII_ALPHA,POPIII_MCUT);
  

  //create the CDF from IMF
  for(i=0; i<N_VALS-1; i++){
    M1 = LOWER+DM*i;
    M2 = LOWER+DM*(i+1);

    masses[i] = M1 * Msun_to_Internal;
    sum += (IMF_norm * trap(fn,1000,M1,M2,POPIII_ALPHA,POPIII_MCUT)) / N;
    prob[i] = sum;
  }
  masses[N_VALS-1] = UPPER * Msun_to_Internal;
  prob[N_VALS-1] = 1.0;

  if(ThisTask==0){
    printf("\n");
    //printf("POPIII gasM %e\n",All.OrigGasMass);
    printf("POPIII Mass %e\n",POPIII_MASS);
    printf("POPIII Norm %e\n",IMF_norm);
    printf("POPIII N    %e\n",N);
    printf("\n");   
    fflush(stdout);
  }
  
}

#ifdef POPII_PROXY
int ionized(int i)
{
    srand(P[i].ID);
    double Qion,ion_rand;
    int ionized = 1;
    int not_ionized = 0;

    if (1./All.Time-1.<= 6) Qion=1.0;
    else Qion = exp(-1.*pow((1./All.Time-1.-6),2)/5.0);  //from Greif:06 5.0(late), 25.0 (mid), 57(early)

    //ion_rand = drand48();
    ion_rand = rand()/(double)RAND_MAX;
    //printf("************Ionization Check: %0.5f\t%0.5f ===> ",Qion,ion_rand );
    if(ion_rand <= Qion) {
      //printf("IONIZED\n");
      return ionized;   // region is ionized therefore can NOT form Pop III star
    }
    else {
      //printf("NOT IONIZED\n");
      return not_ionized;                   //  Pop III star formation CAN happen
    }
}
#endif

//linearly interpolate
double GetMass(int i){
  int n,k=0;
  double rand, di, mu;

  rand = drand48();
  
  for(n=0; n<N_VALS-1; n++){
    if(rand >= prob[n] && rand < prob[n+1]){
      k = n;
      break;
    }
  }

  di = prob[k+1] - prob[k];
  mu = (rand - prob[k])/di;
  
  return (masses[k] * (1.0-mu) + masses[k+1]*mu);
}

//trapezoidal integration
double trap(double (*f)(double, double, double),int N, double a, double b, double alpha, double mcut)
{
   double sum=0.0;
   int i=0;
   double start=a;
   double step=(b-a)/N;
   for(i=0; i<N; i++)
   {
      sum+=(f(start+i*step,alpha,mcut)+f(start+(i+1)*step,alpha,mcut))*step/2.;
   }

   return sum;
}

// imf function (mass)
double fm (double x,double alpha,double mcut)
{ 
   return x*pow(x,-1.*alpha)*exp(-1.*mcut/(x*x));             
} 

// imf function (number)
double fn (double x,double alpha,double mcut)
{ 
   return pow(x,-1.*alpha)*exp(-1.*mcut/(x*x));               
}

double Legacy_PopIII(int i){

  double radius, radius_Ifront, log_radius;
  int n, j, numngb, numngb_Ifront, startnode, dummy;
  double sumwk=0.0, norm_j;
  //double dt=0.0, sm_popIII=0.0, tot_sm_popIII = 0.0, SFR_PopIII=0.0, SFR_PopIII_tot;
  double sum_popIII=0.0;
  double Msun_to_Internal  = All.HubbleParam * SOLAR_MASS / All.UnitMass_in_g;
  double massj_fact, meanweight, energy_increase;
  double ufac = All.UnitPressure_in_cgs/All.UnitDensity_in_cgs;
  double densfac = All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * All.cf_a3inv;

  //printf("######### LP3 called at P[index]: %d with denisty of : %e\n",i,(SphP[i].Density * densfac * 0.76) / PROTONMASS);

#ifdef POPII_PROXY
  if(ionized(i)) return 0.0;
#endif
  MyFloat *Ngbweight;


  /* VALUES TO DISTRIBUTE */
  double sum=0.0, ran_mass=0.0;
  int N_Type2 = 0,N_PISN = 0,N_BH = 0,N_Other = 0;
  double Mass_Type2=0.0,Mass_PISN=0.0,Mass_BH=0.0,Mass_Other=0.0;
  do{
    ran_mass = GetMass(i);

    if(ran_mass < 8.0*Msun_to_Internal){
      N_Other += 1;
      Mass_Other += ran_mass;
    }
    else if(8.0*Msun_to_Internal <= ran_mass && ran_mass < 40.0*Msun_to_Internal){
      N_Type2 += 1;
      Mass_Type2 += ran_mass;
    }
    else if(40.0*Msun_to_Internal <= ran_mass && ran_mass < 140.0*Msun_to_Internal){
      N_BH += 1;
      Mass_BH += ran_mass;
    }
    else{
      N_PISN += 1;
      Mass_PISN += ran_mass;
    }    
    sum += ran_mass;
  }while(sum < POPIII_MASS * Msun_to_Internal);

  
  //sum up energy
  double total_energy = 0.0;
  total_energy += N_Type2 * ENERGY_TYPE2;
  total_energy += N_PISN  * ENERGY_PISN;
  total_energy += N_BH    * ENERGY_BH;
  total_energy += N_Other * ENERGY_OTHER;
  total_energy /= All.UnitEnergy_in_cgs;

  //sum up mass
  double total_metal_mass = 0.0;
  total_metal_mass += Mass_Type2 * YIELD_TYPE2;
  total_metal_mass += Mass_PISN  * YIELD_PISN;
  total_metal_mass += Mass_BH    * YIELD_BH;
  total_metal_mass += Mass_Other * YIELD_OTHER;

  //determine radius
  //radius  = 1.0;       //kpc/h PHYSICAL
  log_radius = -17.02 + 0.383 * log10(total_energy*All.UnitEnergy_in_cgs); //fitting formula from Jaacks+2016 in [pc/h] physical 
  radius = pow(10,log_radius) / 1000.;
  radius /= All.Time;  //to    COMOVING

  /* NEIGHBOR STUFF */
  Ngblist = (int *) mymalloc("popIII list",NumPart * sizeof(int));  
  Ngbweight = (MyFloat *) mymalloc("popIII weights",NumPart * sizeof(MyFloat));
  
  startnode = All.MaxPart;
  numngb = ngb_treefind_variable_threads(P[i].Pos, radius, -1, &startnode, 0, &dummy, &dummy, &dummy, Ngblist);


  //weight by a kernel, otherwise distribute evenly
#ifdef POPIII_KERNEL_DIST
  double wk=0.0, dwk=0.0;
  double dx,dy,dz,ngb_r2,r2,u,rinv4,rinv3;
  double radius_squared = radius*radius;

  //loop through neighbors and calculate kernel contribution
  for( n=0; n<numngb; n++ ){
    j = Ngblist[n];  //index of neighbor n

    dx = P[i].Pos[0] - P[j].Pos[0];
    dy = P[i].Pos[1] - P[j].Pos[1];
    dz = P[i].Pos[2] - P[j].Pos[2];
#ifdef PERIODIC   /*  distance to the closest image in the given box size  */
    NEAREST_XYZ(dx,dy,dz,1);
#endif
                                                
    ngb_r2 = dx*dx + dy*dy + dz*dz;
    r2 = radius_squared;
    
    u = sqrt(ngb_r2/r2);
    rinv4 = 1/(r2*r2);
    rinv3 = rinv4*radius;
    kernel_main(u, rinv3, rinv4, &wk, &dwk, -1);
    sumwk += wk;
    Ngbweight[n] = wk;
  }
#else
  for( n=0; n<numngb; n++ )
    Ngbweight[n] = 1.0;
  sumwk = (double)numngb;
#endif

  
  /*    if total_energy > 0: loop through neighbors and give METAL feedback  else: skip metals and only do thermal    */
  if(total_energy>0.0){
    for( n=0; n<numngb; n++){
        j = Ngblist[n];
        norm_j = 1./(double)numngb; // Ngbweight[n]/sumwk //

        /* Metals only!  Temp now has larger radius handled below */
        massj_fact = 1./(P[j].Mass+total_metal_mass * norm_j);
        //no weighting by kernel
        
        //add to PopIII metal tracer
        P[j].Metallicity[2] = (P[j].Mass * P[j].Metallicity[2] + total_metal_mass*norm_j) * massj_fact;
        //add to Total metal tracer
        P[j].Metallicity[0] = (P[j].Mass * P[j].Metallicity[0] + total_metal_mass*norm_j) * massj_fact;
        
      }
  }  
  myfree(Ngbweight);
  myfree(Ngblist);

  double mean_density, I_front_radius, temp_density;
  double HII_dens_ref = HII_DENSITY_POPIII  * PROTONMASS / (HYDROGEN_MASSFRAC * All.UnitDensity_in_cgs * All.HubbleParam*All.HubbleParam);

  //Determine neighbors for fixed thermal radius defined above (I-front).
  radius_Ifront = I_FRONT_RADIUS / All.Time;  //to COMOVING
  startnode=All.MaxPart;
  Ngblist_Ifront = (int *) mymalloc("popIII I-front",NumPart * sizeof(int));
  numngb_Ifront = ngb_treefind_variable_threads(P[i].Pos, radius_Ifront, -1, &startnode, 0, &dummy, &dummy, &dummy, Ngblist_Ifront);

  I_front_radius = I_FRONT_RADIUS;
  for(int k=0; k<HII_DENSITY_ITERATION; k++){
    temp_density = 0.0;
    mean_density = 0.0;
    for(n=0; n<numngb_Ifront; n++){
      j = Ngblist_Ifront[n];
      temp_density += 1.0/(SphP[j].Density * All.cf_a3inv);
    }
    mean_density = numngb_Ifront/temp_density;
    I_front_radius = I_front_radius * pow(mean_density/HII_dens_ref, -2.0/3.0);
    myfree(Ngblist_Ifront);

    radius_Ifront = I_front_radius / All.Time;  //to COMOVING
    startnode = All.MaxPart;
    Ngblist_Ifront = (int *) mymalloc("PopII radius",NumPart * sizeof(int));
    numngb_Ifront = ngb_treefind_variable_threads(P[i].Pos, radius_Ifront, -1, &startnode, 0, &dummy, &dummy, &dummy, Ngblist_Ifront);
  }

  /*                    loop through neighbors and give THERMAL feedback              */
  for( n=0; n<numngb_Ifront; n++){
      j = Ngblist_Ifront[n];

      // Temperature - converting temperature to energy to add to each particle
      //meanweight = 4. / (1 + 3 * HYDROGEN_MASSFRAC);
      //energy_increase = 1. / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * TEMP_BUBBLE * P[j].Mass * All.UnitMass_in_g;
      //SphP[j].InternalEnergy += energy_increase / All.UnitEnergy_in_cgs / P[j].Mass;
      energy_increase = (TEMP_BUBBLE*BOLTZMANN)/(GAMMA_MINUS1*PROTONMASS *ufac);
      //SphP[j].InternalEnergy += energy_increase / All.UnitEnergy_in_cgs / P[j].Mass;
      SphP[j].InternalEnergy += energy_increase; //test to set all to 5e4K
      SphP[j].InternalEnergyPred = SphP[j].InternalEnergy;
      SphP[j].Pressure = get_pressure(j);

      SphP[j].P3_heating_time = POPIII_PE_HEATING; //added with P2L model

      //Tnew = convert_u_to_temp(SphP[j].InternalEnergy/ufacinv, Particle_density_for_energy_i(j)*densfac,&SphP[j].Ne,j);
      // Don't believe this is needed as the Ne is updated during the next cooling timestep based on the new temp
      /* not sure if this is needed but seems to be done each time feedback 
      is implemented (i.e. metals.c) to update the electron abundance */
      //convert_u_to_temp(SphP[j].InternalEnergy*ufac, Particle_density_for_energy_i(j)*densfac,&SphP[j].Ne,j);

      
      // Ionization fraction
      //SphP[j].Ne = E_FRAC_BUBBLE;
  }
  
  printf("######### POPIII PID[%d]: %d %d %d %d %d %d %e \n",
        P[i].ID,numngb,numngb_Ifront, N_Other, N_Type2, N_BH, N_PISN, sum / Msun_to_Internal);  //mass has little h included here
  fflush(stdout);
  int star_id_type;
  star_id_type = N_Other*1000 + N_Type2*100 + N_BH*10 + N_PISN; // this encodes the Star_Type with info reagarding stellar population
  P[i].Star_Type += star_id_type;

  //write to output, each processor has own output file
  FILE *FdPopIII;
  char buf[200];
  sprintf(buf, "%sPopIII", All.OutputDir);
  mkdir(buf, 02755);
  sprintf(buf, "%sPopIII/popIII.%d.txt", All.OutputDir, ThisTask);
  if(!(FdPopIII = fopen(buf, "a"))){
    printf("error in opening file '%s'\n",buf);
    endrun(1);
  }
  
  /* P[index] #ngb #ngbI  #OTHER #Type2 #BH #PISN  Etot  Msun Mmass  radius  z */
  fprintf(FdPopIII,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%e\t%e\t%e\t%e\t%e\t%e\n",
          i,numngb,numngb_Ifront,N_Other,N_Type2,N_BH,N_PISN,
          total_energy * All.UnitEnergy_in_cgs,  //ergs
          sum / Msun_to_Internal, total_metal_mass / Msun_to_Internal,     //Msun little h included
          radius * All.Time, radius_Ifront * All.Time,                     //PHYSICAL kpc/h
          1./All.Time - 1.);
  fflush(FdPopIII);

  myfree(Ngblist_Ifront);


  return (Mass_Other + Mass_Type2 + Mass_PISN + Mass_BH) / Msun_to_Internal;     // total mass in Msun no little h
   
}

#endif
