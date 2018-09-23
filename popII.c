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


#if defined(POPII_P2L) && defined(GALSF) && defined(METALS)

//ergs per event
#define ENERGY_TYPE2 1.0e51
#define ENERGY_PISN  1.0e52
#define ENERGY_BH    0.0
#define ENERGY_OTHER 0.0

// Bubble properties
#define HII_DENSITY_ITERATION 5
#define HII_DENSITY_POPII 1.0
#define I_FRONT_RADIUS 0.240     //kpc/h PHYSICAL
#define TEMP_BUBBLE 2.0e4  //Kelvin
#define P2_METAL_RADIUS 1.0     //kpc/h PHYSICAL
#define P2_METAL_MASS 16.0     //Msun

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




void P2L(int i, int mode)
{
  double mass_of_star;
  int n, j, numngb_metal, numngb_Ifront, startnode, dummy;
  double energy_radius, metal_radius, energy_increase, total_metal_mass;
  double Msun_to_Internal  = All.HubbleParam * SOLAR_MASS / All.UnitMass_in_g;
  double massj_fact, norm_j;
  double ufac = All.UnitPressure_in_cgs/All.UnitDensity_in_cgs;
  double densfac = All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * All.cf_a3inv;
  double mean_density, I_front_radius, temp_density;
  double HII_dens_ref = HII_DENSITY_POPII  * PROTONMASS / (HYDROGEN_MASSFRAC * All.UnitDensity_in_cgs * All.HubbleParam*All.HubbleParam);
  //printf("############################### POPII P2L called at PID: %d, heating_time = %e, P2L_flag = %0.01f ################################## \n",i,SphP[i].PE_heating_time,SphP[i].P2L_flag);
  //fflush(stdout);

  // Only identify particle within radius of initial star burst and set PE_heating_time.  NO FEEDBACK
  if(mode==0){
  energy_radius = I_FRONT_RADIUS / All.Time;  //to COMOVING
  startnode = All.MaxPart;
  Ngblist_Ifront = (int *) mymalloc("PopII radius",NumPart * sizeof(int));
  numngb_Ifront = ngb_treefind_variable_threads(P[i].Pos, energy_radius, -1, &startnode, 0, &dummy, &dummy, &dummy, Ngblist_Ifront);

  /*                    loop through neighbors and set PE_heating_time             */
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

    energy_radius = I_front_radius / All.Time;  //to COMOVING
    startnode = All.MaxPart;
    Ngblist_Ifront = (int *) mymalloc("PopII radius",NumPart * sizeof(int));
    numngb_Ifront = ngb_treefind_variable_threads(P[i].Pos, energy_radius, -1, &startnode, 0, &dummy, &dummy, &dummy, Ngblist_Ifront);
  }

  for( n=0; n<numngb_Ifront; n++){
      j = Ngblist_Ifront[n];
      SphP[j].P2_heating_time = POPII_PE_HEATING;
  }
  myfree(Ngblist_Ifront);
}//end mode=0


// if mode=1 then do the thermal and metal feedback at star death
if(mode==1){
  /*                              P2L THERMAL                            */
  energy_radius = I_FRONT_RADIUS / All.Time;  //to COMOVING
  startnode = All.MaxPart;
  Ngblist_Ifront = (int *) mymalloc("PopII radius",NumPart * sizeof(int));
  numngb_Ifront = ngb_treefind_variable_threads(P[i].Pos, energy_radius, -1, &startnode, 0, &dummy, &dummy, &dummy, Ngblist_Ifront);

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

    energy_radius = I_front_radius / All.Time;  //to COMOVING
    startnode = All.MaxPart;
    Ngblist_Ifront = (int *) mymalloc("PopII radius",NumPart * sizeof(int));
    numngb_Ifront = ngb_treefind_variable_threads(P[i].Pos, energy_radius, -1, &startnode, 0, &dummy, &dummy, &dummy, Ngblist_Ifront);
  }

  /*                    loop through neighbors and give THERMAL feedback              */
  for( n=0; n<numngb_Ifront; n++){
      j = Ngblist_Ifront[n];
      // no feedback on original star formation particle
      if(j!=i){
      //SphP[j].InternalEnergy += energy_increase / All.UnitEnergy_in_cgs / P[j].Mass;
      energy_increase = (TEMP_BUBBLE*BOLTZMANN)/(GAMMA_MINUS1*PROTONMASS *ufac);
      SphP[j].InternalEnergy += energy_increase; 
      SphP[j].InternalEnergyPred = SphP[j].InternalEnergy;
      SphP[j].Pressure = get_pressure(j);
    }
  }
  myfree(Ngblist_Ifront);


  /*                              P2L METAL                            */
  metal_radius = P2_METAL_RADIUS / All.Time;
  Ngblist = (int *) mymalloc("popII metal",NumPart * sizeof(int));  
  startnode = All.MaxPart;
  numngb_metal = ngb_treefind_variable_threads(P[i].Pos, metal_radius, -1, &startnode, 0, &dummy, &dummy, &dummy, Ngblist);
  total_metal_mass = P2_METAL_MASS * Msun_to_Internal;

  for( n=0; n<numngb_metal; n++){
      j = Ngblist[n];
      // do not add metals to the original star formation particle!!
      if(j!=i){
      norm_j = 1./(double)numngb_metal;
      massj_fact = 1./(P[j].Mass+total_metal_mass * norm_j);

      //add to PopII metal tracer
      P[j].Metallicity[1] = (P[j].Mass * P[j].Metallicity[1] + total_metal_mass*norm_j) * massj_fact;
      //add to Total metal tracer
      P[j].Metallicity[0] = (P[j].Mass * P[j].Metallicity[0] + total_metal_mass*norm_j) * massj_fact;
      }
    }  
  myfree(Ngblist);
}//end mode=1

}//end void

#endif
