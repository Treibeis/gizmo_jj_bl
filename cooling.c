#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

#include "./cooling.h"

/*
 * This file contains the routines for optically-thin cooling (generally aimed towards simulations of the ISM, 
 *   galaxy formation, and cosmology). A wide range of heating/cooling processes are included, including 
 *   free-free, metal-line, Compton, collisional, photo-ionization and recombination, and more. Some of these 
 *   are controlled by individual modules that need to be enabled or disabled explicitly.
 *
 * This file was originally part of the GADGET3 code developed by
 *   Volker Springel (volker.springel@h-its.org). The code has been modified heavily by 
 *   Phil Hopkins (phopkins@caltech.edu) for GIZMO; everything except the original metal-free free-free and 
 *   photo-ionization heating physics has been added (or re-written), and the iteration routine to converge to 
 *   temperatures has been significantly modified.
 */


#ifdef COOLING



#define NCOOLTAB  2000

#define JAMPL 1.0   /* amplitude factor relative to input table */
#define TABLESIZE 250   /* Max # of lines in TREECOOL */

static float inlogz[TABLESIZE];
static float gH0[TABLESIZE], gHe[TABLESIZE], gHep[TABLESIZE];
static float eH0[TABLESIZE], eHe[TABLESIZE], eHep[TABLESIZE];
static int nheattab;    /* length of table */

#define SMALLNUM 1.0e-60
#define COOLLIM  0.2
#define HEATLIM	 10.0

#define eV_to_K   11606.0
#define eV_to_erg 1.60184e-12

/* CAFG: H number density above which we assume no ionizing bkg (proper cm^-3) */
#define NH_SS 0.0123

static double yhelium = (1. - HYDROGEN_MASSFRAC) / (4. * HYDROGEN_MASSFRAC);
static double mhboltz = PROTONMASS / BOLTZMANN;         /* hydrogen mass over Boltzmann constant */
static double XH = HYDROGEN_MASSFRAC;

static double ethmin;		/* minimum internal energy for neutral gas */

static double Tmin = 0.4313637642;	/* in log10 */
static double Tmax = 6.0;
static double CMB_Temp=2.73;
static double deltaT;
static double Color_Temp = 20000;  //SSP color temp

static double *BetaH0, *BetaHep, *Betaff;
static double *AlphaHp, *AlphaHep, *Alphad, *AlphaHepp;
static double *GammaeH0, *GammaeHe0, *GammaeHep;
#ifdef COOL_METAL_LINES_BY_SPECIES
/* if this is enabled, the cooling table files should be in a folder named 'spcool_tables' in the run directory.
 cooling tables can be downloaded at: https://dl.dropbox.com/u/16659252/spcool_tables.tgz */
static float *SpCoolTable0;
static float *SpCoolTable1;
#endif

static double J_UV = 0, gJH0 = 0, gJHep = 0, gJHe0 = 0, epsH0 = 0, epsHep = 0, epsHe0 = 0;

static double ne, necgs, nHcgs;
static double bH0, bHep, bff, aHp, aHep, aHepp, ad, geH0, geHe0, geHep;
static double gJH0ne, gJHe0ne, gJHepne;
static double nH0, nHp, nHep, nHe0, nHepp;
static double k[30],ny_out[12];

static double DoCool_u_old_input, DoCool_rho_input, DoCool_dt_input, DoCool_ne_guess_input;




/* this is just a simple loop if all we're doing is cooling (no star formation) */
void cooling_only(void)
{
    int i;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].Type == 0 && P[i].Mass > 0)
        {
            do_the_cooling_for_particle(i);
        } // if(P[i].Type == 0 && P[i].Mass > 0)
    } // for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
} // void cooling_only(void)





/* subroutine which actually sends the particle data to the cooling routine and updates the entropies */
void do_the_cooling_for_particle(int i)
{
    double unew;
    double dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval;
    double dtime = dt / All.cf_hubble_a; /*  the actual time-step */
   //printf("%e\t%e\n",dt,dtime);

    if((P[i].TimeBin)&&(dt>0)&&(P[i].Mass>0)&&(P[i].Type==0))  // upon start-up, need to protect against dt==0 //
    {
        
        double ne = SphP[i].Ne;	/* electron abundance (gives ionization state and mean molecular weight) */
// #ifdef PRIMORDIAL_COOLING
//         ne = SphP[i].Primordial_Chem[5];
// #endif
        double uold = DMAX(All.MinEgySpec, SphP[i].InternalEnergy);

#ifdef GALSF_FB_HII_HEATING
        double u_to_temp_fac = PROTONMASS / BOLTZMANN * GAMMA_MINUS1 * All.UnitEnergy_in_cgs / All.UnitMass_in_g;
        double uion = HIIRegion_Temp / u_to_temp_fac;
        if(SphP[i].DelayTimeHII > 0) if(uold<uion) uold=uion; /* u_old should be >= ionized temp if used here */
#endif // GALSF_FB_HII_HEATING
        
#ifndef COOLING_OPERATOR_SPLIT
        /* do some prep operations on the hydro-step determined heating/cooling rates before passing to the cooling subroutine */
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        /* calculate the contribution to the energy change from the mass fluxes in the gravitation field */
        double grav_acc; int k;
        for(k = 0; k < 3; k++)
        {
            grav_acc = All.cf_a2inv * P[i].GravAccel[k];
#ifdef PMGRID
            grav_acc += All.cf_a2inv * P[i].GravPM[k];
#endif
            SphP[i].DtInternalEnergy -= SphP[i].GravWorkTerm[k] * All.cf_atime * grav_acc;
        }
#endif
        /* limit the magnitude of the hydro dtinternalenergy */
       double du = SphP[i].DtInternalEnergy * dtime;
	     if(du < (COOLLIM-1)*SphP[i].InternalEnergy) {SphP[i].DtInternalEnergy = (COOLLIM-1)*SphP[i].InternalEnergy / dtime;}
	     if(du >  (HEATLIM-1)*SphP[i].InternalEnergy) {SphP[i].DtInternalEnergy =  (HEATLIM-1)*SphP[i].InternalEnergy / dtime;}


#ifndef GRACKLE
        /* and convert to cgs before use in the cooling sub-routine */
        SphP[i].DtInternalEnergy *= All.HubbleParam * All.UnitEnergy_in_cgs / (All.UnitMass_in_g * All.UnitTime_in_s) * (PROTONMASS/XH);
#endif
#endif // COOLING_OPERATOR_SPLIT

#ifndef RT_COOLING_PHOTOHEATING
        unew = DoCooling(uold, Particle_density_for_energy_i(i) * All.cf_a3inv, dtime, &ne, i); /* All.cf_a3inv * Density_code = Density_physical */
#else
        double fac_entr_to_u = pow(SphP[i].Density * All.cf_a3inv, GAMMA_MINUS1) / GAMMA_MINUS1;
        unew = uold + dt * fac_entr_to_u * (rt_DoHeating(i, dt) + rt_DoCooling(i, dt));
#endif // RT_COOLING_PHOTOHEATING
        
        
#ifdef GALSF_FB_HII_HEATING
        /* set internal energy to minimum level if marked as ionized by stars */
        if(SphP[i].DelayTimeHII > 0)
        {
            if(unew<uion)
            {
                unew=uion;
                if(SphP[i].DtInternalEnergy<0) SphP[i].DtInternalEnergy=0;
                //if(SphP[i].dInternalEnergy<0) SphP[i].dInternalEnergy=0; //manifest-indiv-timestep-debug//
            }
            SphP[i].Ne = 1.0 + 2.0*yhelium;
        }
#endif // GALSF_FB_HII_HEATING
        
        
#if defined(BH_THERMALFEEDBACK)
        if(SphP[i].Injected_BH_Energy)
		{
            unew += SphP[i].Injected_BH_Energy / P[i].Mass;
            SphP[i].Injected_BH_Energy = 0;
		}
#endif
        

#if defined(COSMIC_RAYS) && !defined(COSMIC_RAYS_DISABLE_COOLING)
        /* cosmic ray interactions affecting the -thermal- temperature of the gas are included in the actual cooling/heating functions; 
            they are solved implicitly above. however we need to account for energy losses of the actual cosmic ray fluid, here. The 
            timescale for this is reasonably long, so we can treat it semi-explicitly, as we do here.
            -- We use the estimate for combined hadronic + Coulomb losses from Volk 1996, Ensslin 1997, as updated in Guo & Oh 2008: */
        double ne_cgs = ((0.78 + 0.22*ne*XH) / PROTONMASS) * (SphP[i].Density * All.cf_a3inv * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam);
        double CR_coolingrate_perunitenergy = -7.51e-16 * ne_cgs * (All.UnitTime_in_s / All.HubbleParam); // converts cgs to code units //
        double CR_Egy_new = SphP[i].CosmicRayEnergyPred * exp(CR_coolingrate_perunitenergy * dtime);
        SphP[i].CosmicRayEnergyPred = SphP[i].CosmicRayEnergy = CR_Egy_new;
#endif
        
        /* InternalEnergy, InternalEnergyPred, Pressure, ne are now immediately updated; however, if COOLING_OPERATOR_SPLIT
         is set, then DtInternalEnergy carries information from the hydro loop which is only half-stepped here, so is -not- updated. 
         if the flag is not set (default), then the full hydro-heating is accounted for in the cooling loop, so it should be re-zeroed here */
        
     

        SphP[i].InternalEnergy = unew;
        SphP[i].Ne = ne;
        SphP[i].InternalEnergyPred = SphP[i].InternalEnergy;
        SphP[i].Pressure = get_pressure(i);


#ifndef COOLING_OPERATOR_SPLIT
        SphP[i].DtInternalEnergy = 0;
#endif
        
        
#ifdef GALSF_FB_HII_HEATING
        /* count off time which has passed since ionization 'clock' */
        if(SphP[i].DelayTimeHII > 0) SphP[i].DelayTimeHII -= dtime;
        if(SphP[i].DelayTimeHII < 0) SphP[i].DelayTimeHII = 0;
#endif // GALSF_FB_HII_HEATING
        
    } // closes if((P[i].TimeBin)&&(dt>0)&&(P[i].Mass>0)&&(P[i].Type==0)) check
}


/* checks if the particle temp is below the Jeans Temperature
 * if so update the internal energy and abundances
 * rho should be in code units (physical)
 * u   should be in code units
 */
int JeansMinT(double *u, double rho, double *ne_guess, int target)
{
#ifdef GALSF_JEANS_MIN_T
  double logTmincool, logTmaxcool, logT;
  
  if( rho > All.MinRhoThresh && All.cf_atime > 0.02 ) {
    logTmincool = 4 + 0.3333333 * log10(rho / All.MinRhoThresh);
    logTmaxcool = logTmincool + 0.5;
    logT = LogTemp(*u * (All.UnitPressure_in_cgs / All.UnitDensity_in_cgs),*ne_guess,target);
    if( logT < logTmaxcool ){	// only adjust energy for particles below this T
      // set internal energy to ensure Jeans mass is resolved
      *u *= pow(10., logTmincool - logT);
      // set ionization so that the mass-weighted average is the JMT temperature
      logT = convert_u_to_temp(*u * (All.UnitPressure_in_cgs / All.UnitDensity_in_cgs), 
			       rho * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam, 
			       ne_guess, target);
      return 1;
    }
  }
#endif
  return 0;
}

/* returns new internal energy per unit mass. 
 * Arguments are passed in code units, density is proper density.
 */

double DoCooling(double u_old, double rho, double dt, double *ne_guess, int target)
{
  double u, du;
  double u_lower, u_upper, fix_up;
  double ratefact;
  double LambdaNet;
  int iter=0, iter_upper=0, iter_lower=0, nyct=0;

// #ifdef GRACKLE
// #ifndef COOLING_OPERATOR_SPLIT
//   // ## OUR ROUTINE ##
//   u = u_old+SphP[target].DtInternalEnergy*dt;
//   double energy=u;
//   if( ! JeansMinT(&energy, rho, ne_guess, target) )
//     energy = CallGrackle(u, rho, dt, ne_guess, target, 0);
//   if( energy > HEATLIM * u ) energy = HEATLIM * u;
//   if( energy < COOLLIM * u ) energy = COOLLIM * u;
//   /*double lognH = log10(rho * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * HYDROGEN_MASSFRAC / PROTONMASS);
//   double logTold = LogTemp(u_old * (All.UnitPressure_in_cgs / All.UnitDensity_in_cgs),*ne_guess);
//   double logTnew = LogTemp(energy * (All.UnitPressure_in_cgs / All.UnitDensity_in_cgs),*ne_guess);
//   if( (logTnew<3.8 || logTold<3.8) && lognH > -1 ) printf("DOCOOLING: P%d %d %d %g %g %g %g %g %g %d %g\n",ThisTask,target,P[target].ID,1./All.Time-1,lognH,logTold,logTnew,SphP[target].DtInternalEnergy*dt/u_old,rho/All.MinRhoThresh,SphP[target].N_windLaunches,log10(All.MinRhoThresh* All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam * HYDROGEN_MASSFRAC / PROTONMASS));
//   if( lognH!=lognH || logTold!=logTold ) endrun(5801);*/
//   return DMAX(energy,All.MinEgySpec);

//   // ## PHIL'S ROUTINE ##  (CURRENTLY BYPASSED by above return statement)
//      because grackle uses a pre-defined set of libraries, we can't properly incorporate the hydro heating
//      into the cooling subroutine. instead, we will use the approximate treatment below
//      to split the step 
//     du = dt * SphP[target].DtInternalEnergy / (All.HubbleParam * All.UnitEnergy_in_cgs / (All.UnitMass_in_g * All.UnitTime_in_s) * (PROTONMASS/XH));
//     u_old += 0.5*du;
//     u = CallGrackle(u_old, rho, dt, ne_guess, target, 0);
//     /* now we attempt to correct for what the solution would have been if we had included the remaining half-step heating
//      term in the full implicit solution. The term "r" below represents the exact solution if the cooling function has
//      the form d(u-u0)/dt ~ -a*(u-u0)  around some u0 which is close to the "ufinal" returned by the cooling routine,
//      to which we then add the heating term from hydro and compute the solution over a full timestep */
//     double r=u/u_old; if(r>1) {r=1/r;} if(fabs(r-1)>1.e-4) {r=(1-r)/log(r);} r=DMAX(0,DMIN(r,1));
//     du *= 0.5*r; if(du<-0.5*u) {du=-0.5*u;} u+=du;
//     /* with full operator splitting we just call grackle normally. note this is usually fine,
//      but can lead to artificial noise at high densities and low temperatures, especially if something
//      like artificial pressure (but not temperature) floors are used such that the temperature gets
//      'contaminated' by the pressure terms */
//     u = CallGrackle(u_old, rho, dt, ne_guess, target, 0);
// #endif
//     return DMAX(u,All.MinEgySpec);
// #endif
  
  // if(JeansMinT(&u_old, rho, ne_guess, target))
  //   return u_old;


  DoCool_u_old_input = u_old;
  DoCool_rho_input = rho;
  DoCool_dt_input = dt;
  DoCool_ne_guess_input = *ne_guess;


  rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;	/* convert to physical cgs units */
  u_old *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;
  dt *= All.UnitTime_in_s / All.HubbleParam;

  nHcgs = XH * rho / PROTONMASS;	/* hydrogen number dens in cgs units */
  ratefact = nHcgs * nHcgs / rho;

  u = u_old;
  u_lower = u;
  u_upper = u;


  LambdaNet = CoolingRateFromU(u, rho, ne_guess, target, dt);

  /* bracketing */

  if(u - u_old - ratefact * LambdaNet * dt < 0)	/* heating */
    {
      u_upper *= sqrt(1.1);
      u_lower /= sqrt(1.1);
      while((iter_upper<MAXITER)&&(u_upper - u_old - ratefact * CoolingRateFromU(u_upper, rho, ne_guess, target, dt) * dt < 0))
	{
	  u_upper *= 1.1;
	  u_lower *= 1.1;
        iter_upper++;
	}

    }

  if(u - u_old - ratefact * LambdaNet * dt > 0)
    {
      u_lower /= sqrt(1.1);
      u_upper *= sqrt(1.1);
      while((iter_lower<MAXITER)&&(u_lower - u_old - ratefact * CoolingRateFromU(u_lower, rho, ne_guess, target, dt) * dt > 0))
	{
	  u_upper /= 1.1;
	  u_lower /= 1.1;
        iter_lower++;
	}
    }

  do
    {
      u = 0.5 * (u_lower + u_upper);

      LambdaNet = CoolingRateFromU(u, rho, ne_guess, target, dt);

      if(u - u_old - ratefact * LambdaNet * dt > 0)
	{
	  u_upper = u;
	}
      else
	{
	  u_lower = u;
	}

      du = u_upper - u_lower;

      iter++;

      if(iter >= (MAXITER - 10))
	printf("u= %g\n", u);
    }
    while(((fabs(du/u) > 1.0e-3)||((fabs(du/u) > 1.0e-6)&&(iter < 10))) && (iter < MAXITER));
    //while(fabs(du / u) > 1.0e-6 && iter < MAXITER);

  if(iter >= MAXITER)
    {
      printf("failed to converge in DoCooling()\n");
      printf("DoCool_u_old_input=%g\nDoCool_rho_input= %g\nDoCool_dt_input= %g\nDoCool_ne_guess_input= %g\n",
	     DoCool_u_old_input, DoCool_rho_input, DoCool_dt_input, DoCool_ne_guess_input);
      endrun(10);
    }

    
#ifdef PRIMORDIAL_COOLING
  double mu,xn,meanweight,CMB_u, TP;
  //mu = (1 + 4 * yhelium) / (1 + yhelium + *ne_guess);
  mu = (1 + 4 * yhelium) / (1 + yhelium + SphP[target].Primordial_Chem[5]);
  meanweight = mu * PROTONMASS;
  xn =  rho / meanweight;  /* total number dens in cgs units (physical)*/
  TP = convert_u_to_temp(u, rho, ne_guess, target);//GAMMA_MINUS1 / BOLTZMANN * u * meanweight;

  /* check that minimum internal energy is >= cmb floor if not set = cmb floor 
                              this is a temporary fix for the particle "drip" issue*/
  CMB_u = (1.0 / GAMMA_MINUS1) * (BOLTZMANN / meanweight) * (CMB_Temp/All.Time-1.0);
  if(u < CMB_u && xn > 1.0 && (1.0/All.Time-1. < 30.0)) u=CMB_u;
  //fix_up = (1.0 / GAMMA_MINUS1) * (BOLTZMANN / meanweight) * (CMB_Temp/All.Time+100.0);
  //if(u > fix_up && xn > 100.0*All.HubbleParam && (1.0/All.Time-1. < 30.0)) u=CMB_u;

  //if(ThisTask == 0) printf("@@@@@@@@@@@@@@@@@@@@@@%e\t%e\n",xn,TP);

  /* get new abundance ratios using final Temp at end of cooling; update Primordial_Chem array */
  P_Chem(target, dt , xn , rho, TP, 0.0, 0.0, ne_guess);  
  for(nyct=0;nyct<12;nyct++) SphP[target].Primordial_Chem[nyct]=ny_out[nyct];
  
  // if(xn < 1.e3 && TP < 1.e4 && (1.0/All.Time-1. < 15.)){
  //   P_Chem(target, dt , xn , rho, TP, 0.0, 0.0, ne_guess); 
  //   for(nyct=0;nyct<12;nyct++) SphP[target].Primordial_Chem[nyct]=ny_out[nyct];
  // } // this if statement is here to keep the chemistry routine from running on high temp/density particles to avoid slow down.
  
  // if(xn > 1.e3 && TP > 1.e4){
  //   u *= All.UnitDensity_in_cgs / All.UnitPressure_in_cgs;  /* back to internal units */
  //   return u;
  // }
  // else{
  //   P_Chem(target, dt , xn , rho, TP, 0.0, 0.0, ne_guess); 
  //   for(nyct=0;nyct<12;nyct++) SphP[target].Primordial_Chem[nyct]=ny_out[nyct];
  // }


#endif

  u *= All.UnitDensity_in_cgs / All.UnitPressure_in_cgs;	/* back to internal units */
  return u;
}



/* returns cooling time. 
 * NOTE: If we actually have heating, a cooling time of 0 is returned.
 */
 // double GetCoolingTime(double u_old, double rho, double *ne_guess, int target)
 // {
 //   double u;
 //   double ratefact;
 //   double LambdaNet, coolingtime;

 //   /*  currently not functional for EFFECTIVE_EQS
 // #ifdef GRACKLE
 //   coolingtime = CallGrackle(u_old, rho, 0.0, ne_guess, target, 1);
 //   if(coolingtime >= 0) coolingtime = 0.0;
 //   coolingtime *= All.HubbleParam / All.UnitTime_in_s;
 //   return coolingtime;
 // #endif
 //   */

 //   DoCool_u_old_input = u_old;
 //   DoCool_rho_input = rho;
 //   DoCool_ne_guess_input = *ne_guess;

 //   rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;	/* convert to physical cgs units */
 //   u_old *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;

 //   nHcgs = XH * rho / PROTONMASS;	/* hydrogen number dens in cgs units */
 //   ratefact = nHcgs * nHcgs / rho;
 //   u = u_old;
 //   LambdaNet = CoolingRateFromU(u, rho, ne_guess, target, dt);

 //   /* bracketing */

 //   if(LambdaNet >= 0)		/* ups, we have actually heating due to UV background */
 //     return 0;

 //   coolingtime = u_old / (-ratefact * LambdaNet);

 //   coolingtime *= All.HubbleParam / All.UnitTime_in_s;

 //   return coolingtime;
 // }


/* returns new internal energy per unit mass. 
 * Arguments are passed in code units, density is proper density.
 */
// double DoInstabilityCooling(double m_old, double u, double rho, double dt, double fac, double *ne_guess, int target)
// {
//   double m, dm;
//   double m_lower, m_upper;
//   double ratefact;
//   double LambdaNet;
//   int iter = 0;

//   DoCool_u_old_input = u;
//   DoCool_rho_input = rho;
//   DoCool_dt_input = dt;
//   DoCool_ne_guess_input = *ne_guess;

//   if(fac <= 0)			/* the hot phase is actually colder than the cold reservoir! */
//     {
//       return 0.01 * m_old;
//     }

//   rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;	/* convert to physical cgs units */
//   u *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;
//   dt *= All.UnitTime_in_s / All.HubbleParam;
//   fac *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

//   nHcgs = XH * rho / PROTONMASS;	/* hydrogen number dens in cgs units */
//   ratefact = nHcgs * nHcgs / rho * fac;

//   m = m_old;
//   m_lower = m;
//   m_upper = m;

//   LambdaNet = CoolingRateFromU(u, rho, ne_guess, target,dt);

//   /* bracketing */

//   if(m - m_old - m * m / m_old * ratefact * LambdaNet * dt < 0)	/* heating */
//     {
//       m_upper *= sqrt(1.1);
//       m_lower /= sqrt(1.1);
//       while(m_upper - m_old -
//           m_upper * m_upper / m_old * ratefact * CoolingRateFromU(u, rho * m_upper / m_old,
//                                                                   ne_guess, target, dt) * dt < 0)
//       {
// 	m_upper *= 1.1;
// 	m_lower *= 1.1;
//       }
//     }

//   if(m - m_old - m_old * ratefact * LambdaNet * dt > 0)
//     {
//       m_lower /= sqrt(1.1);
//       m_upper *= sqrt(1.1);
//       while(m_lower - m_old -
//           m_lower * m_lower / m_old * ratefact * CoolingRateFromU(u, rho * m_lower / m_old,
//                                                                   ne_guess, target, dt) * dt > 0)
//       {
// 	m_upper /= 1.1;
// 	m_lower /= 1.1;
//       }
//     }

//   do
//     {
//       m = 0.5 * (m_lower + m_upper);

//         LambdaNet = CoolingRateFromU(u, rho * m / m_old, ne_guess, target, dt);

//       if(m - m_old - m * m / m_old * ratefact * LambdaNet * dt > 0)
// 	{
// 	  m_upper = m;
// 	}
//       else
// 	{
// 	  m_lower = m;
// 	}

//       dm = m_upper - m_lower;

//       iter++;

//       if(iter >= (MAXITER - 10))
// 	printf("m= %g\n", m);
//     }
//   while(fabs(dm / m) > 1.0e-6 && iter < MAXITER);

//   if(iter >= MAXITER)
//     {
//       printf("failed to converge in DoCooling()\n");
//       printf("DoCool_u_old_input=%g\nDoCool_rho_input= %g\nDoCool_dt_input= %g\nDoCool_ne_guess_input= %g\n",
// 	     DoCool_u_old_input, DoCool_rho_input, DoCool_dt_input, DoCool_ne_guess_input);
//       printf("m_old= %g\n", m_old);
//       endrun(11);
//     }

//   return m;
// }






/* this function determines the electron fraction, and hence the mean 
 * molecular weight. With it arrives at a self-consistent temperature.
 * Element abundances and the rates for the emission are also computed
 */
double convert_u_to_temp(double u, double rho, double *ne_guess, int target)
{
  double temp, temp_old, temp_new, max = 0, ne_old;
  double mu, meanweight;
  int iter = 0;

  double u_input, rho_input, ne_input;

#ifdef PRIMORDIAL_COOLING
  
  mu = (1 + 4 * yhelium) / (1 + yhelium + SphP[target].Primordial_Chem[5]);
  meanweight = mu * PROTONMASS;

  temp = GAMMA_MINUS1 / BOLTZMANN * u * meanweight;
  
  if(temp<=0) temp=pow(10.0,Tmin);
  if(log10(temp)<Tmin) temp=pow(10.0,Tmin);
  if(log10(temp)>Tmax) temp=pow(10.0,Tmax);
  //find_abundances_and_rates(log10(temp), rho, ne_guess, target);
  //*ne_guess = SphP[target].Primordial_Chem[5];
  
  return temp;
#endif

  u_input = u;
  rho_input = rho;
  ne_input = *ne_guess; 

  mu = (1 + 4 * yhelium) / (1 + yhelium + *ne_guess);
  temp = GAMMA_MINUS1 / BOLTZMANN * u * PROTONMASS * mu;

  do
    {
      ne_old = *ne_guess;

      find_abundances_and_rates(log10(temp), rho, ne_guess, target);
      
      temp_old = temp;

      mu = (1 + 4 * yhelium) / (1 + yhelium + *ne_guess);

      temp_new = GAMMA_MINUS1 / BOLTZMANN * u * PROTONMASS * mu;

      max =
	DMAX(max,
	     temp_new / (1 + yhelium + *ne_guess) * fabs((*ne_guess - ne_old) / (temp_new - temp_old + 1.0)));

      temp = temp_old + (temp_new - temp_old) / (1 + max);
      iter++;

      if(iter > (MAXITER - 10))
	  printf("-> temp= %g ne=%g\n", temp, *ne_guess);
    }
    while(
          ((fabs(temp - temp_old) > 0.1 * temp) ||
           ((fabs(temp - temp_old) > 1.0e-3 * temp) && (temp > 200.))) && iter < MAXITER);

  if(iter >= MAXITER)
      {
	printf("failed to converge in convert_u_to_temp(%d)\n",target);
	printf("u_input= %g\nrho_input=%g\n ne_input=%g\n", u_input, rho_input, ne_input);
	printf
	  ("DoCool_u_old_input=%g\nDoCool_rho_input= %g\nDoCool_dt_input= %g\nDoCool_ne_guess_input= %g\n",
	   DoCool_u_old_input, DoCool_rho_input, DoCool_dt_input, DoCool_ne_guess_input);

	endrun(12);
      }

    if(temp<=0) temp=pow(10.0,Tmin);
    if(log10(temp)<Tmin) temp=pow(10.0,Tmin);

  return temp;
}



/* this function computes the actual abundance ratios 
 */
void find_abundances_and_rates(double logT, double rho, double *ne_guess, int target)
{
  double neold, nenew;
  int j, niter;
  double flow, fhi, t;
  double logT_input, rho_input, ne_input;
  double NH_SS_z=NH_SS, shieldfac;

  logT_input = logT;
  rho_input = rho;
  ne_input = *ne_guess;

  if(isnan(logT)) logT=Tmin;    /* nan trap (just in case) */
    
  if(logT <= Tmin)		/* everything neutral */
    {
      nH0 = 1.0;
      nHe0 = yhelium;
      nHp = 0;
      nHep = 0;
      nHepp = 0;
      ne = 0;
      *ne_guess = 0;
      return;
    }

  if(logT >= Tmax)		/* everything is ionized */
    {
      nH0 = 0;
      nHe0 = 0;
      nHp = 1.0;
      nHep = 0;
      nHepp = yhelium;
      ne = nHp + 2.0 * nHepp;
      *ne_guess = ne;		/* note: in units of the hydrogen number density */
      return;
    }

  t = (logT - Tmin) / deltaT;
  j = (int) t;
    if(j<0){j=0;}
    if(j>NCOOLTAB){
        printf("warning: j>NCOOLTAB : j=%d t %g Tlow %g Thi %g logT %g Tmin %g deltaT %g \n",j,t,Tmin+deltaT*j,Tmin+deltaT*(j+1),logT,Tmin,deltaT);fflush(stdout);
        j=NCOOLTAB;
    }
  fhi = t - j;
  flow = 1 - fhi;

  if(*ne_guess == 0)
    *ne_guess = 1.0;

    double local_gammamultiplier=1;
#ifdef GALSF_FB_LOCAL_UV_HEATING
    if ((target >= 0) && (gJH0 > 0))
    {
    local_gammamultiplier = SphP[target].RadFluxEUV * 2.29e-10; // converts to GammaHI for typical SED (rad_uv normalized to Habing)
    local_gammamultiplier = 1 + local_gammamultiplier/gJH0;
    }
#endif
    
    /* CAFG: this is the density that we should use for UV background threshold */
    nHcgs = XH * rho / PROTONMASS;	/* hydrogen number dens in cgs units */
    if(gJH0>0)
        NH_SS_z = NH_SS*pow(local_gammamultiplier*gJH0/1.0e-12,0.66)*pow(10.,0.173*(logT-4.));
    else
        NH_SS_z = NH_SS*pow(10.,0.173*(logT-4.));
    if(nHcgs<100.*NH_SS_z) shieldfac=exp(-nHcgs/NH_SS_z); else shieldfac=0;
#ifdef COOL_LOW_TEMPERATURES
    if(logT < Tmin+1) shieldfac *= (logT-Tmin); // make cutoff towards Tmin more continuous //
    //shieldfac *= 1 - exp(Tmin-logT);
#endif
#ifdef GALSF_EFFECTIVE_EQS
    shieldfac = 1; // self-shielding is implicit in the sub-grid model already //
#endif

  ne = *ne_guess;
  neold = ne;
  niter = 0;
  necgs = ne * nHcgs;

  /* evaluate number densities iteratively (cf KWH eqns 33-38) in units of nH */
  do
    {
      niter++;

      aHp = flow * AlphaHp[j] + fhi * AlphaHp[j + 1];
      aHep = flow * AlphaHep[j] + fhi * AlphaHep[j + 1];
      aHepp = flow * AlphaHepp[j] + fhi * AlphaHepp[j + 1];
      ad = flow * Alphad[j] + fhi * Alphad[j + 1];
      geH0 = flow * GammaeH0[j] + fhi * GammaeH0[j + 1];
      geHe0 = flow * GammaeHe0[j] + fhi * GammaeHe0[j + 1];
      geHep = flow * GammaeHep[j] + fhi * GammaeHep[j + 1];
#ifdef COOL_LOW_TEMPERATURES
        // make cutoff towards Tmin more continuous //
        if(logT < Tmin+1) {
            geH0 *= (logT-Tmin);
            geHe0 *= (logT-Tmin);
            geHep *= (logT-Tmin);
        }
#endif

      if(necgs <= 1.e-25 || J_UV == 0)
	{
	  gJH0ne = gJHe0ne = gJHepne = 0;
	}
      else
	{
        /* CAFG: if density exceeds NH_SS, ignore ionizing background. */
        gJH0ne = gJH0 * local_gammamultiplier / necgs * shieldfac;
        gJHe0ne = gJHe0 * local_gammamultiplier / necgs * shieldfac;
        gJHepne = gJHep * local_gammamultiplier / necgs * shieldfac;
	}

      nH0 = aHp / (aHp + geH0 + gJH0ne);	/* eqn (33) */
      nHp = 1.0 - nH0;		/* eqn (34) */

      if((gJHe0ne + geHe0) <= SMALLNUM)	/* no ionization at all */
	{
	  nHep = 0.0;
	  nHepp = 0.0;
	  nHe0 = yhelium;
	}
      else
	{
	  nHep = yhelium / (1.0 + (aHep + ad) / (geHe0 + gJHe0ne) + (geHep + gJHepne) / aHepp);	/* eqn (35) */
	  nHe0 = nHep * (aHep + ad) / (geHe0 + gJHe0ne);	/* eqn (36) */
	  nHepp = nHep * (geHep + gJHepne) / aHepp;	/* eqn (37) */
	}

      neold = ne;

      ne = nHp + nHep + 2 * nHepp;	/* eqn (38) */
      necgs = ne * nHcgs;

      if(J_UV == 0)
	break;

      nenew = 0.5 * (ne + neold);
      ne = nenew;
      necgs = ne * nHcgs;

      if(fabs(ne - neold) < 1.0e-4)
	break;

      if(niter > (MAXITER - 10))
	printf("ne= %g  niter=%d\n", ne, niter);
    }
  while(niter < MAXITER);

  if(niter >= MAXITER)
    {
      printf("no convergence reached in find_abundances_and_rates()\n");
      printf("logT_input= %g  rho_input= %g  ne_input= %g\n", logT_input, rho_input, ne_input);
      printf("DoCool_u_old_input=%g\nDoCool_rho_input= %g\nDoCool_dt_input= %g\nDoCool_ne_guess_input= %g\n",
	     DoCool_u_old_input, DoCool_rho_input, DoCool_dt_input, DoCool_ne_guess_input);
      endrun(13);
    }

#if (GRACKLE_CHEMISTRY >  0) // non-tabular; this is used to set abundances with JMT
  SphP[target].grHI = nH0;
  SphP[target].grHII = nHp;
  SphP[target].grHeI = nHe0;
  SphP[target].grHeII = nHep;
  SphP[target].grHeIII = nHepp;
  SphP[target].grHM  = 1.0e-20; //cannot be zero
#endif

  bH0 = flow * BetaH0[j] + fhi * BetaH0[j + 1];
  bHep = flow * BetaHep[j] + fhi * BetaHep[j + 1];
  bff = flow * Betaff[j] + fhi * Betaff[j + 1];

  *ne_guess = ne;
}




/*  this function first computes the self-consistent temperature
 *  and abundance ratios, and then it calculates 
 *  (heating rate-cooling rate)/n_h^2 in cgs units 
 */
double CoolingRateFromU(double u, double rho, double *ne_guess, int target, double dt)
{
  /* rho is passed in physical cgs units, dt is system timestep in seconds */
  double Temp, mu, ncgs,xn, meanweight,nHcgs2,dt_years;
  dt_years = dt / 3.154e7;  //system timestep in years

  Temp = convert_u_to_temp(u,rho,ne_guess,target);
  //mu = (1 + 4 * yhelium) / (1 + yhelium + *ne_guess);
  mu = (1 + 4 * yhelium) / (1 + yhelium + SphP[target].Primordial_Chem[5]);
  //meanweight = 4.0/(3 * XH + 1 + 4 * XH * *ne_guess) * PROTONMASS;
  //meanweight = mu * PROTONMASS;   
  meanweight = mu * PROTONMASS;
  xn =  rho / meanweight;  /* total number dens in cgs units (physical)*/
  nHcgs = XH * rho / PROTONMASS;
  nHcgs2 = nHcgs*nHcgs;

#ifdef PRIMORDIAL_COOLING
  double Lambda, lam_mol, lam_line,lam_metal,Heating,PE2_Heat,PE3_Heat, UVB_Heat;
  //  first update chemical abundances based on particle internal energy (temp) and dt 
  // P_Chem(target, dt , xn, rho, Temp, 0.0, 0.0, ne_guess); 


  /* molecular cooling rate */
  lam_mol = hmcool(target, 2.7e0/All.Time, Temp, xn);

  /* H, He line cooling + CI + IC + Bremsstrahlung rates */
  lam_line = hlcool(target, 1./All.Time-1.0, Temp, xn);

#ifdef JJ_METAL_COOLING
  if(P[target].Metallicity[0] == 0.0) lam_metal=0.0;
  else lam_metal = JJ_metal_cooling(target,Temp,xn);
#else
  lam_metal=0.0;
#endif

//we do this P3 since metal enrichment is instantaneous...negates metal impact until after PE heating is done.
#if defined(JJ_METAL_COOLING) && defined(POPIII_PE_HEATING)
if(SphP[target].P3_heating_time > 0.0) lam_metal=0.0;  
#endif 

  /* total cooling rate */
  Lambda = lam_line/nHcgs2 + lam_mol/nHcgs2 + lam_metal/nHcgs2;

#ifdef POPII_P2L
  if(SphP[target].P2_heating_time > 0.0){
    PE2_Heat = PE2_heating(xn);
  }
  else PE2_Heat = 0.0;

  if(SphP[target].P3_heating_time > 0.0){
    PE3_Heat = PE3_heating(xn);
  }
  else PE3_Heat = 0.0;
#else
  PE3_Heat = 0.0;
  PE2_Heat = 0.0;
#endif

#ifdef POPII_UVB
  if( 1. / All.Time - 1 > 25.) UVB_Heat=0.0;
  else UVB_Heat = UVB_heating(xn,target);
#else
  UVB_Heat = 0.0;
#endif

//  REALLY HATE units here!!  dtIE is [erg s^-1] and PE & UVB are [erg s^-1 cm^-3] but need [erg/s cm^3] thus /nHgs and /nHgs2
  Heating = SphP[target].DtInternalEnergy/nHcgs + PE2_Heat/nHcgs2 + PE3_Heat/nHcgs2 + UVB_Heat/nHcgs2;    

  /* return Heating - Cooling */
  //if(target==10000) printf("Cooling Check: %g\t%g\t%g\n",lam_mol/nHcgs2,lam_line/nHcgs2,Heating);

  //return -1 * Lambda + (SphP[target].DtInternalEnergy*ratefact_inv);
  return Heating - Lambda;

#endif

#ifndef PRIMORDIAL_COOLING
    /* Normal Gizmo cooling rate */ 
    return CoolingRate(log10(Temp), rho, ne_guess, target);
#endif
}

//photo-electric heating
double PE2_heating(double xnh){
  //double alphaB,psi,Tc,kb;
  //alphaB = 2.59e-13;  //case B recombination
  //psi = 1.38;
  double Tc = 20000;  //ssp temp

  // [cm^3/s] * [cm^-6] * [erg/K] * [K]  ===>  [erg s^-1 cm^-3]
  return  2.59e-13 * (xnh * xnh) * 1.38 * BOLTZMANN * Tc;
}

double PE3_heating(double xnh){
  //double alphaB,psi,Tc,kb;
  //alphaB = 2.59e-13;  //case B recombination
  //psi = 1.38;
  double Tc = 30000;  //ssp temp
  
  // [cm^3/s] * [cm^-6] * [erg/K] * [K]  ===>  [erg s^-1 cm^-3]
  return  2.59e-13 * (xnh * xnh) * 1.38 * BOLTZMANN * Tc;
}

double UVB_heating(double xnh, int target){
  int i, ilow;
  double logz, dzlow, dzhi;
  double redshift;
  double zeta, tau;
  double xnHI,io_frac;
  
  io_frac = 1.0 - SphP[target].Primordial_Chem[5];
  if(io_frac <= 0.0) return 0.0; // if completely ionized no UV heating
  
  xnHI = io_frac * xnh;  //neutral H fraction


  if(All.ComovingIntegrationOn)
    redshift = 1 / All.Time - 1;
  else
    {
    /* in non-cosmological mode, still use, but adopt z=0 background */
    redshift = 0;
    /*
         gJHe0 = gJHep = gJH0 = 0;
         epsHe0 = epsHep = epsH0 = 0;
         J_UV = 0;
         return;
    */
    }

  logz = log10(redshift + 1.0);
  ilow = 0;
  for(i = 0; i < nheattab; i++)
    {
      if(inlogz[i] < logz)
  ilow = i;
      else
  break;
    }

  dzlow = logz - inlogz[ilow];
  dzhi = inlogz[ilow + 1] - logz;


  if(logz > inlogz[nheattab - 1] || gH0[ilow] == 0 || gH0[ilow + 1] == 0 || nheattab == 0)
    {
      return 0.0;
    }

  //gJH0 = JAMPL * pow(10., (dzhi * log10(gH0[ilow]) + dzlow * log10(gH0[ilow + 1])) / (dzlow + dzhi));
  zeta = pow(10., (dzhi * log10(gH0[ilow]) + dzlow * log10(gH0[ilow + 1])) / (dzlow + dzhi)); //same as done below

  tau = 6.30e-18 * xnh * 3.086e21;  
  if(tau>10.) return 0.0;

  // zeta [s^-1] * nh [cm^-3] * Kb [erg/K] * Tc [K]
  return   exp(-tau) * zeta * xnHI  * 1.38 * BOLTZMANN * Color_Temp;  //erg s^-1 cm^-3
}


double JJ_metal_cooling(int i, double temp, double xn)
/*=======================================================================
** Metal fine structure cooling from CII, OI, SiII and FeII
** added by JdJ 03/2017 based on CSS code and equations in Maio:07

  xn - total number density in physical cgs
=======================================================================*/
{
  
  double CII_gamma_H_21, CII_gamma_H_12, CII_gamma_e_21, CII_gamma_e_12;
  double SiII_gamma_H_21, SiII_gamma_H_12, SiII_gamma_e_21, SiII_gamma_e_12;
  double OI_gamma_H_21, OI_gamma_H_12, OI_gamma_e_21, OI_gamma_e_12;
  double FeII_gamma_H_21, FeII_gamma_H_12, FeII_gamma_e_21, FeII_gamma_e_12;
  double CII_A21, CII_delta_E21, CII_g1, CII_g2;
  double SiII_A21, SiII_delta_E21, SiII_g1, SiII_g2;
  double OI_A21, OI_delta_E21, OI_g1, OI_g2;
  double FeII_A21, FeII_delta_E21, FeII_g1, FeII_g2;
  double lamda_metal,lambda_CII,lambda_SiII,lambda_OI,lambda_FeII;
  double ne,nh,xe,beta;
  double x_CII,x_SiII,x_OI,x_FeII;  //n_x/n_h
  double n_CII,n_SiII,n_OI,n_FeII; 

  nh = xn * 0.93;
  xe = SphP[i].Primordial_Chem[5];   // ne/nh
  //Asplund:09 solar abundances

  x_CII = 3.26e-3 * P[i].Metallicity[0];    //3.26e-3
  x_OI = 8.65e-3 * P[i].Metallicity[0];    //8.65e-3
  x_SiII = 1.08e-3 * P[i].Metallicity[0];    //1.08e-3
  x_FeII = 1.73e-3 * P[i].Metallicity[0];    //1.73e-3
  
  //total number density for each species
  n_CII = x_CII * nh;
  n_OI = x_OI * nh;
  n_SiII = x_SiII * nh;
  n_FeII = x_FeII * nh;

  //if Z=0 or T>20000 K no need to calculate lambda metals
  if (P[i].Metallicity[0]==0.0 || temp>20000.0) return 0.0;

  // CII
  CII_gamma_H_21 = 8e-10 * pow((temp/100.0),0.07);
  CII_gamma_e_21 = 2.8e-7 * pow((temp/100.0),-0.5);
  CII_A21 = 2.4e-6;
  CII_delta_E21 = 1.259e-14;
  CII_g2=4.;
  CII_g1=2.;

  beta = pow((BOLTZMANN*temp),-1.0);
  CII_gamma_e_12 = CII_g2/CII_g1 * CII_gamma_e_21 * exp(-1*beta*CII_delta_E21);
  CII_gamma_H_12 = CII_g2/CII_g1 * CII_gamma_H_21 * exp(-1*beta*CII_delta_E21);

  lambda_CII = lam_2level(CII_gamma_H_21,CII_gamma_H_12,CII_gamma_e_21,CII_gamma_e_12,CII_A21,CII_delta_E21,xe,nh,n_CII);
  /////////////////////

  // SiII
  SiII_gamma_H_21 = 8e-10 * pow((temp/100.0),-0.07);
  SiII_gamma_e_21 = 1.7e-6 * pow((temp/100.0),-0.5);
  SiII_A21 = 2.1e-4;
  SiII_delta_E21 = 5.71e-14;
  SiII_g2=4.;
  SiII_g1=2.;

  
  SiII_gamma_e_12 = SiII_g2/SiII_g1 * SiII_gamma_e_21 * exp(-1*beta*SiII_delta_E21);
  SiII_gamma_H_12 = SiII_g2/SiII_g1 * SiII_gamma_H_21 * exp(-1*beta*SiII_delta_E21);

  lambda_SiII = lam_2level(SiII_gamma_H_21,SiII_gamma_H_12,SiII_gamma_e_21,SiII_gamma_e_12,SiII_A21,SiII_delta_E21,xe,nh,n_SiII);
  ////////////////////

  // OI 1-->2  NOTE: just treating OI and FeII as 2-level system...need to update to 3-level but should be minimal impact
  OI_gamma_H_21 = 9.2e-11 * pow((temp/100.0),0.67);
  OI_gamma_e_21 = 1.4e-8; 
  OI_A21 = 8.9e-5;
  OI_delta_E21 = 3.144e-14;
  OI_g2=5.;
  OI_g1=3.;

  OI_gamma_e_12 = OI_g2/OI_g1 * OI_gamma_e_21 * exp(-1*beta*OI_delta_E21);
  OI_gamma_H_12 = OI_g2/OI_g1 * OI_gamma_H_21 * exp(-1*beta*OI_delta_E21);

  lambda_OI = lam_2level(OI_gamma_H_21,OI_gamma_H_12,OI_gamma_e_21,OI_gamma_e_12,OI_A21,OI_delta_E21,xe,nh,n_OI);
  ////////////////////

  // FeI 1-->2
  FeII_gamma_H_21 = 9.5e-10; 
  FeII_gamma_e_21 = 1.8e-6 * pow((temp/100.),-0.5); 
  FeII_A21 = 2.15e-3;
  FeII_delta_E21 = 7.64e-14;
  FeII_g2=8.;
  FeII_g1=10.;

  FeII_gamma_e_12 = FeII_g2/FeII_g1 * FeII_gamma_e_21 * exp(-1*beta*FeII_delta_E21);
  FeII_gamma_H_12 = FeII_g2/FeII_g1 * FeII_gamma_H_21 * exp(-1*beta*FeII_delta_E21);

  lambda_FeII = lam_2level(FeII_gamma_H_21,FeII_gamma_H_12,FeII_gamma_e_21,FeII_gamma_e_12,FeII_A21,FeII_delta_E21,xe,nh,n_FeII);
  ////////////////////
  

  lamda_metal = lambda_CII + lambda_SiII + lambda_OI + lambda_FeII;

  //return lamda_metal/(nh*nh);
  return lamda_metal;
}

double lam_2level(double gamma_H_21,double gamma_H_12,double gamma_e_21,double gamma_e_12,double A21,double delta_E21,double xe,double nh,double ntot)
{
  double top, bottom;
  
  top = gamma_H_12 + gamma_e_12*xe;
  bottom = gamma_H_12 + gamma_H_21 + (gamma_e_12 + gamma_e_21)*xe + A21/nh;

  return top/bottom * (ntot * A21 * delta_E21);
}


double  hmcool(int i, double Tcmb, double temp, double xn)
/*=======================================================================
** Evaluates H2 cooling function!
**    include H2 cooling
**    (Galli & Palla 1998)
=======================================================================*/
{

     
  double LH2,LLTE,XNCRIT,LLOW,T3, lamH2, lamHD;
  double LOGLOW,LRLTE,LVLTE;
  double T3C,LOGLOWC,templogc,templog;
  double LLOWC,LRLTEC,LVLTEC,LLTEC,XNCRITC;
  double w_hd,dw_hddT,XHDm,XHDg,w_hdg,xdummy;
  double xnh2,xnh,xnd,xnhd;
  double w_hdC,hdLTE,hdLTEC,g10,g10C,g21,g21C;
  double kB,E10,E21,hdLOW,hdLOWC,lamHDC,lamhm;

  kB = 1.3806e-16;   
  //temp=P[i].Temp;
  //xnh=P[i].nh;
  //xnh2=P[i].H2I*P[i].nh;
  //xnh =  SphP[i].Primordial_Chem[0] * rho / PROTONMASS ;  /* hydrogen number dens in cgs units */
  xnh = xn * 0.93;
  xnd = 4.3e-5 * xnh;
  xnh2  = SphP[i].Primordial_Chem[3] * xnh;  /* H2 number dens in cgs units */
  xnhd  = SphP[i].Primordial_Chem[11] * xnd;  /* HD number dens in cgs units */
/*=======================================================================
**    include H2 cooling
**    (Galli & Palla 1998)
=======================================================================*/
  T3=temp/1000.e0;
  T3C=Tcmb/1000.e0;
  templog=log10(temp);
  templogc=log10(Tcmb);
  LOGLOW=-103.0 + 97.59*templog - 48.05*pow(templog,2) +
        10.80*pow(templog,3) - 0.9032*pow(templog,4); 
  LOGLOWC=-103.0 + 97.59*templogc - 48.05*pow(templogc,2) +
        10.80*pow(templogc,3) - 0.9032*pow(templogc,4);
  LLOW=pow(10.e0,LOGLOW);
  LLOWC=pow(10.e0,LOGLOWC);
  LRLTE=(9.5e-22*pow(T3,3.76)*exp(-pow(0.13/T3,3))/
        (1.e0+0.12*pow(T3,2.1))+3.e-24*exp(-0.51/T3))/xnh;
  LRLTEC=(9.e-22*pow(T3C,3.76)*exp(-pow(0.13/T3C,3))/
        (1.e0+0.12*pow(T3C,2.1))+3.e-24*exp(-0.51/T3C))/xnh;
  LVLTE=(6.7e-19*exp(-5.86/T3) + 1.6e-18*exp(-11.7/T3))/xnh;
  LVLTEC=(6.7e-19*exp(-5.86e0/T3C) + 
         1.6e-18*exp(-11.7e0/T3C))/xnh;
  LLTE=LRLTE + LVLTE;
  LLTEC=LRLTEC + LVLTEC;
  XNCRIT=xnh*(LLTE/LLOW);
  XNCRITC=xnh*(LLTEC/LLOWC);
  LH2=LLTE/(1.e0+XNCRIT/xnh) - LLTEC/(1.e0+XNCRITC/xnh);

  lamH2= xnh2 * xnh * LH2;  //units: erg cm^-3 s^-1
  //lamH2= LH2;

/*==================================================================
**     include HD cooling 
**     (Galli & Palla 1998, Coppola et al. 2011)
==================================================================*/
  
  w_hd = -55.5725+56.649*templog-37.9102*pow(templog,2)+12.698*pow(templog,3)-2.02424*pow(templog,4)+0.122393*pow(templog,5);
  hdLTE = pow(10.e0,w_hd);
  w_hdC = -55.5725+56.649*templogc-37.9102*pow(templogc,2)+12.698*pow(templogc,3)-2.02424*pow(templogc,4)+0.122393*pow(templogc,5);
  hdLTEC = pow(10.e0,w_hdC);
  g10 = (4.4e-12+3.6e-13*pow(temp,0.77))/1.27;
  g21 = (4.1e-12+2.1e-13*pow(temp,0.92))/1.27;
  g10C = (4.4e-12+3.6e-13*pow(Tcmb,0.77))/1.27;
  g21C = (4.1e-12+2.1e-13*pow(Tcmb,0.92))/1.27;
  E10 = 128.e0;
  E21 = 255.e0;
  hdLOW = (2.e0*g10*E10*kB*exp(-E10/temp)+5.0/3.0*g21*E21*kB*exp(-E21/temp))*xnh;
  hdLOWC = (2.e0*g10C*E10*kB*exp(-E10/Tcmb)+5.0/3.0*g21C*E21*kB*exp(-E21/Tcmb))*xnh;
  lamHD = hdLOW*hdLTE/(hdLOW+hdLTE);
  lamHDC = hdLOWC*hdLTEC/(hdLOWC+hdLTEC);
  lamhm = xnhd * (lamHD - lamHDC) + lamH2;

  return lamhm; //lamH2;
}

double hlcool(int i, double zred, double temp, double xn)
/**======================================================================
*** Evaluates atomic line cooling due to H and He!
*** (see Cen 1992, ApJS, 78, 341)
***====================================================================*/
{
  double xnhe,mu,LAM[11], lamHl, ny[9];
  double gff,XH1,XH2,XH3,XH4,XH5,T3,T5,T6;
  double Tcmb,xnh, meanweight;

  //temp=P[i].Temp;
  //xnh=P[i].nh;
  /* rho is passed in physical units */
  //mu = (1 + 4 * yhelium) / (1 + yhelium + *ne_guess);
  //mu = 1.27;
  //xnh = rho /  (PROTONMASS*mu);
  //xnh = XH * rho /  PROTONMASS;  /* hydrogen number dens in physical cgs units */ 
  //xnh = SphP[i].Primordial_Chem[0];
  //meanweight = 1.22 * PROTONMASS;
  xnh = xn * 0.93;
  xnhe=0.07*xnh/0.93;
  //xnhe = yhelium * xnh;
  //xnhe   = SphP[i].Primordial_Chem[6]*nHcgs;

  Tcmb= 2.7e0 * (1.e0 + zred);

  ny[0]=SphP[i].Primordial_Chem[0] * xnh;     // HI
  ny[1]=SphP[i].Primordial_Chem[1] * xnh;     // HII
  ny[2]=0.e0;                                 // H-
  ny[3]=0.e0;                                 // H2
  ny[4]=0.e0;                                 // H2+
  ny[5]=SphP[i].Primordial_Chem[5] * xnh;     // e-

  //ny[6]=SphP[i].Primordial_Chem[6] * xnhe;    // HeI
  ny[7]=SphP[i].Primordial_Chem[7] * xnhe;    // HeII
  ny[8]=SphP[i].Primordial_Chem[8] * xnhe;    // HeIII

/* Bremsstrahlung cooling  */
  if (temp > 5.e3) {
    gff=1.1e0+0.34e0*exp(-pow(5.5e0-log10(temp),2.e0)/3.e0);
    LAM[0]=1.42e-27*gff*sqrt(temp)*(ny[1]+ny[7]
                  +4.e0*ny[8])*ny[5];
  } else {
    LAM[0]=1.e-52;
  }

/* Collisional ionization cooling (HI) */
  T5=temp/1.e5;
  XH1=sqrt(temp)/(1.e0+sqrt(T5));
  if (temp > 5.e3) {
    LAM[1]=1.27e-21*XH1*exp(-157809.1e0/temp)*ny[0]*ny[5];
  } else {
    LAM[1]=1.e-52;
  }

/* Collisional ionization cooling (HeI) */
  if (temp > 8.e3) {
    LAM[2]=9.38e-22*XH1*exp(-285335.4e0/temp)*ny[6]*ny[5];
  } else {
    LAM[2]=1.e-52;
  }

/* Collisional ionization cooling (HeII) */
  if (temp > 1.e4) {
    LAM[3]=4.95e-22*XH1*exp(-631515.0e0/temp)*ny[7]*ny[5];
  } else {
    LAM[3]=1.e-52;
  }

/* Recombination cooling (HII) */
  T3=temp/1.e3;
  T6=temp/1.e6;
  XH4=1.e0/(1.e0+pow(T6,0.7e0));
  if (temp > 5.e3) {
    LAM[4]=8.70e-27*XH4*sqrt(temp)*pow(T3,-0.2e0)*ny[1]*ny[5];
  } else {
    LAM[4]=1.e-52;
  }

/* Recombination cooling (HeII) */
  if (temp > 5.e3) {
    LAM[5]=1.55e-26*pow(temp,0.3647e0)*ny[7]*ny[5];
  } else {
    LAM[5]=1.e-52;
  }

/* Recombination cooling (HeIII) */
  if (temp > 5.e3) {
    LAM[6]=3.48e-26*XH4*sqrt(temp)*pow(T3,-0.2e0)*ny[8]*ny[5];
  } else {
    LAM[6]=1.e-52;
  }

/* Dielectronic recombination cooling  */
  XH5=exp(-470000.e0/temp)*(1.e0+0.3e0*exp(-94000.e0/temp));
  if (temp > 5.e3) {
    LAM[7]=1.24e-13*XH5*pow(temp,-1.5e0)*ny[7]*ny[5];
  } else {
    LAM[7]=1.e-52;
  }

/* Collisional excitation cooling (HI) */
  XH2=1.e0/(1.e0+sqrt(T5));
  if (temp > 5.e3) {
    LAM[8]=7.50e-19*XH2*exp(-118348.e0/temp)*ny[0]*ny[5];
  } else {
    LAM[8]=1.e-52;
  }

/* Collisional excitation cooling (HeII) */
  XH3=XH2*pow(temp,-0.397e0);
  if (temp > 8.e3) {
    LAM[9]=5.54e-17*XH3*exp(-473638.e0/temp)*ny[7]*ny[5];
  } else {
    LAM[9]=1.e-52;
  }

/* Compton cooling  */
  LAM[10]=5.4e-36*pow(1.e0+zred,4.e0)*ny[5]*(temp-Tcmb);

  lamHl=LAM[0]+LAM[1]+LAM[2]+LAM[3]+LAM[4]+LAM[5]+LAM[6]+
       LAM[7]+LAM[8]+LAM[9]+LAM[10];
  
  return lamHl; //ergs cm^-3 s^-1

}

void P_Chem(int target,double dt_tot,double xn, double rho,double tempr,double J_21,double del_r,double *ne_guess)  
/****************************************************************
** Solves primordial chemistry network for the 9 species : 
** H,H+,H-,H2,H2+,e-,He,He+,He++
** with approximate BDF method! 
** See methods papers: Bromm, Coppi & Larson 2002, ApJ, 564, 23
**                     Johnson & Bromm 2006, MNRAS, 366, 247
** Deuterium network (species: D,D+,HD) is added.
**
** dt_tot: timestep over which the abundances are to be advanced; to be given
** in physical seconds!
** xn: total number density (physical, NOT comoving), in units of cm^-3
** tempr: Gas temperature (in units of K)
**
**
** Species are indexed as follows:
** ny[0]: HI
** ny[1]: HII
** ny[2]: H-
** ny[3]: H2
** ny[4]: H2+
** ny[5]: e-
** ny[6]: HeI
** ny[7]: HeII
** ny[8]: HeIII
** nyd[0]: DI
** nyd[1]: DII
** nyd[2]: HD
****************************************************************/
{
  double dt,dt_cum;
  double delt_x,XNUM1,XDENOM1,xnd;
  double XNUM2,XDENOM2;
  double kd1,kd3,kd4,kd8,kd10,xnh,xnhe;
  double CHD,DHD,CDp,DDp;
  int kl;
  double nyd[3],ny[9],Cr[9],Ds[9];
  double EPS=1.0e-4;
  double zred = 1.0 / All.Time - 1;
  double Tcmb = 2.7e0 * (1.e0 + zred);
  double eta_LW,alpha_kde,beta_kdi;
  double JLW_popIII,JLW_popII,popII_sfrd,meanweight;



/* xn is TOTAL num density in physical units */
  xnh = 0.93 * xn;
  xnhe=0.07*xnh/0.93;
  xnd = 4.3e-5 * xnh;

  ny[0] = SphP[target].Primordial_Chem[0] * xnh;  // HI
  ny[1] = SphP[target].Primordial_Chem[1] * xnh;  // HII
  ny[2] = SphP[target].Primordial_Chem[2] * xnh;  // H-
  ny[3] = SphP[target].Primordial_Chem[3] * xnh;  // H2
  ny[4] = SphP[target].Primordial_Chem[4] * xnh;  // H2+
  ny[5] = SphP[target].Primordial_Chem[5] * xnh;  // e-

  ny[6] = SphP[target].Primordial_Chem[6] * xnhe; // HeI
  ny[7] = SphP[target].Primordial_Chem[7] * xnhe; // HeII
  ny[8] = SphP[target].Primordial_Chem[8] * xnhe; // HeIII
  
  nyd[0] = SphP[target].Primordial_Chem[9] * xnd;  // DI
  nyd[1] = SphP[target].Primordial_Chem[10] * xnd; // DII
  nyd[2] = SphP[target].Primordial_Chem[11] * xnd; // HD


       
  P_Rates(tempr,xnh,J_21,del_r);

  #ifdef LYMAN_WERNER
  /* first get the shielding factor using Draine & Bertoldi 1996 */
    SphP[target].Primordial_Chem[12] = get_fshield(target,tempr, xn, rho);
    
  /* then calculate J_LW,21 from Johnson 2011 & CSS 2012 */ 
    eta_LW = 2.e4;    // Greif & Bromm 2006
    alpha_kde = 1.71; // from Agarwal 2014
    beta_kdi = 0.97;  // from Agarwal 2014
    //Pop III 
  #ifdef POPIII
    JLW_popIII = 1. * (eta_LW/1.e4) * (All.POPIII_SFRD/1.e-2) * pow((1.+zred)/10.,3);
  #endif
    //Pop II estimate for sfrd, analytical fit in JDJ:16 based on SF:16 obs
    if(zred > 25) popII_sfrd = 0.0;
    else popII_sfrd = 0.15*exp(-0.024*pow(zred+1.,2));   
#ifdef POPII_P2L
    popII_sfrd = All.POPII_SFRD;
#endif
    eta_LW /= 5.0;  //adjustment for Pop II stars having lower photons/baryon
    JLW_popII  = 1. * (eta_LW/1.e4) * (popII_sfrd/1.e-2) * pow((1.+zred)/10.,3); 
    SphP[target].Primordial_Chem[13] = JLW_popIII + JLW_popII;  
    
  /* finally update the k_di and k_de rates */
    k[24] = 1.1e-10 * SphP[target].Primordial_Chem[13] * alpha_kde;
    //k[25] = 0.0;  set below
    //k[26] = 0.0;  set below
    k[27] = 1.38e-12 * beta_kdi * SphP[target].Primordial_Chem[13] * SphP[target].Primordial_Chem[12]; //e.g. Abel 1997
  #endif


  kd1=k[3];
  kd3=3.7e-10*pow(tempr,0.28e0)*exp(-43.e0/tempr);
  kd4=3.7e-10*pow(tempr,0.28e0);
  kd8=2.1e-9;
  kd10=1.e-9*exp(-464.e0/tempr);

  dt_cum=0.e0;
  while (dt_cum < dt_tot) {

    dt=dt_tot;

    delt_x=0.e0;
    while (dt >= delt_x) {
      Cr[5]=k[21]*ny[0]+k[22]*ny[6]+k[23]*ny[7]+
            (k[0]*ny[0]+k[1]*ny[6]+k[2]*ny[7])*ny[5];
      Ds[5]=k[3]*ny[1]+k[4]*ny[7]+
            k[5]*ny[8];
      delt_x=EPS*ny[5]/fabs(Cr[5]-Ds[5]*ny[5]);
      dt =dt/2.e0;
    }

    if ((dt_cum+dt) > dt_tot) {
      dt=dt_tot-dt_cum;
      dt_cum=dt_tot;
    } else {
      dt_cum=dt_cum+dt;
    }

/**** calculate equilibrium abundance for e- *********************/
    Cr[5]=k[21]*ny[0]+k[22]*ny[6]+k[23]*ny[7]+
          (k[0]*ny[0]+k[1]*ny[6]+k[2]*ny[7])*ny[5];
    Ds[5]=k[3]*ny[1]+k[4]*ny[7]+
          k[5]*ny[8];
    ny[5]=(ny[5]+Cr[5]*dt)/(1.e0+Ds[5]*dt);

/**** calculate equilibrium abundance for HI *********************/
    Cr[0]=k[3]*ny[1]*ny[5];
    Ds[0]=k[0]*ny[5]+k[21];
    ny[0]=(ny[0]+Cr[0]*dt)/(1.e0+Ds[0]*dt);

/**** calculate equilibrium abundance for HII *********************/
    Cr[1]=k[0]*ny[5]*ny[0]+k[21]*ny[0];
    Ds[1]=k[3]*ny[5];
    ny[1]=(ny[1]+Cr[1]*dt)/(1.e0+Ds[1]*dt);

/**** calculate equilibrium abundance for HeI *********************/
    Cr[6]=k[4]*ny[7]*ny[5];
    Ds[6]=k[1]*ny[5]+k[22];
    ny[6]=(ny[6]+Cr[6]*dt)/(1.e0+Ds[6]*dt);

/**** calculate equilibrium abundance for HeII *********************/
    Cr[7]=(k[1]*ny[5]+k[22])*ny[6]+k[5]*ny[5]*ny[8];
    Ds[7]=(k[2]+k[4])*ny[5]+k[23];
    ny[7]=(ny[7]+Cr[7]*dt)/(1.e0+Ds[7]*dt);

/**** calculate equilibrium abundance for HIII *********************/
    Cr[8]=k[2]*ny[7]*ny[5]+k[23]*ny[7];
    Ds[8]=k[5]*ny[5];
    ny[8]=(ny[8]+Cr[8]*dt)/(1.e0+Ds[8]*dt);


/**** calculate equilibrium abundance for H- *********************/
    XNUM1=(k[8]*ny[0] + k[13]*ny[3])*ny[5];
    //XNUM1=(k[8]*ny[0]*ny[5] + k[13]*ny[3]);
    XDENOM1=(k[9]+k[19])*ny[0]+(k[12]+k[20])*ny[1]+k[18]*ny[5]
         + k[24];
    ny[2]=XNUM1/XDENOM1;

/**** calculate equilibrium abundance for H2+ ********************/
    XNUM2=(k[6]*ny[0] + k[16]*ny[3]+k[20]*ny[2])*ny[1]
        + k[26]*ny[3];
    XDENOM2=k[7]*ny[0]+k[10]*ny[5]+k[11]*ny[2]+k[25];
    ny[4]=XNUM2/XDENOM2;

/**** calculate equilibrium abundance for H2 *********************/
    Ds[3]=k[13]*ny[5]+k[14]*ny[0]+k[15]*ny[3]
         +k[16]*ny[1]+k[17]*ny[5]+k[26]+k[27];
    Cr[3]=k[7]*ny[4]*ny[0]+k[9]*ny[2]*ny[0]+k[11]*ny[4]*ny[2];
    ny[3]=(ny[3]+Cr[3]*dt)/(1.e0+Ds[3]*dt);
    //tform=ny[3]/Cr[3];
    //tdest=1.e0/Ds[3];

    ny[0]=xnh-2.e0*ny[3]-ny[1]-ny[2]-2.e0*ny[4];
    if (ny[0] < 0.e0) ny[0]=0.e0;
    ny[6]=xnhe-ny[7]-ny[8];
    if (ny[6] < 0.e0) ny[6]=0.e0;
    ny[1]=ny[5]+ny[2]-ny[4]-ny[7]-2.e0*ny[8];
    if (ny[1] < 0.e0) ny[1]=0.e0;
    ny[0]=xnh-2.e0*ny[3]-ny[1]-ny[2]-2.e0*ny[4];
    if (ny[0] < 0.e0) ny[0]=0.e0;

    XNUM1=(k[8]*ny[0] + k[13]*ny[3])*ny[5];
    //XNUM1=(k[8]*ny[0]*ny[5] + k[13]*ny[3]);
    XDENOM1=(k[9]+k[19])*ny[0]+(k[12]+k[20])*ny[1]+k[18]*ny[5]
        + k[24];
    ny[2]=XNUM1/XDENOM1;
    for (kl=0; kl<9; kl++) {
      if (ny[kl] < 1.e-30) ny[kl]=0.e0;
    }

/**** Solve Deuterium network for D, D+, and HD *********************
**   Galli & Palla (1998) -> 'minimal model for D chemistry!'       */

    CDp=kd3*nyd[0]*ny[1]+kd10*nyd[2]*ny[1];
    DDp=kd1*ny[5]+kd4*ny[0]+kd8*ny[3];
    nyd[1]=(nyd[1]+CDp*dt)/(1.e0+DDp*dt);
    CHD=kd8*nyd[1]*ny[3];
    DHD=kd10*ny[1];
    nyd[2]=(nyd[2]+CHD*dt)/(1.e0+DHD*dt);
    nyd[0]=xnd-nyd[2]-nyd[1];
    if (nyd[0] < 1.e-50) {
        nyd[0]=0.e0;
        nyd[1]=xnd-nyd[0]-nyd[2];
        nyd[2]=xnd-nyd[0]-nyd[1];
    }
    if (nyd[1] < 1.e-50) {
        nyd[1]=0.e0;
        nyd[0]=xnd-nyd[1]-nyd[2];
        nyd[2]=xnd-nyd[0]-nyd[1];
    }
    if (nyd[2] < 1.e-50) {
        nyd[2]=0.e0;
    }

    for (kl=0; kl<9; kl++) {
      if (ny[kl] < 0.e0) {
        printf("%6d  %15.6e \n",kl,ny[kl]);
        exit(0); 
      }
    }
  }

  // for (kl=0; kl<9; kl++) {
  //   SphP[target].Primordial_Chem[kl]=ny[kl];
  // }

  // for (kl=9; kl<12; kl++) {
  //   SphP[target].Primordial_Chem[kl]=nyd[kl]; 
  // }

  ny_out[0] = ny[0] / xnh;
  ny_out[1] = ny[1] / xnh;
  ny_out[2] = ny[2] / xnh;
  ny_out[3] = ny[3] / xnh;
  ny_out[4] = ny[4] / xnh;
  ny_out[5] = ny[5] / xnh;
  ny_out[6] = ny[6] / xnhe;
  ny_out[7] = ny[7] / xnhe;
  ny_out[8] = ny[8] / xnhe;
  
  ny_out[9] = nyd[0] / xnd;
  ny_out[10] = nyd[1] / xnd;
  ny_out[11] = nyd[2] / xnd;

  //*ne_guess = ny[5];
}


void P_Rates(double tempr,double xnh,double J_21,double del_r)
/**========================================================================
*** Rate Coefficients k[0]-k[29]:
***======================================================================*/
{
  double k5a,k5b,ncr_1,ncr_2,kh_1,kh_2,kl_1,kl_2;
  double yab[1],tt,t5fac,x,Tgas, del_rlim;
  yab[0]=xnh;
  tt=sqrt(tempr);
  t5fac=1.e0/(1.e0+sqrt(tempr/1.e5));


/*** recombination rate coefficients: Hep from Black (1981, MNRAS 197, 553)
----                              Hp and Hepp from Cen (1992, ApJS 78, 341)*/
  k5a=1.5e-10*pow(tempr,-0.6353e0);
  if (tempr > pow(10.e0,3.6e0)) {
    k5b=1.9e-3*pow(tempr,-1.5e0)*exp(-4.7e5/tempr)*
      (1.e0+0.3e0*exp(-9.4e4/tempr));
  } else {
    k5b = 0.e0;
  }
  k[4] = k5a + k5b;
  k[3]=8.40e-11/tt*pow(tempr/1.e3,-0.2e0)/
     (1.e0+pow(tempr/1.e6,0.7e0));
  k[5]=3.36e-10/tt*pow(tempr/1.e3,-0.2e0)/
     (1.e0+pow(tempr/1.e6,0.7e0));

/*-- collisional ionization rates (from Black, with Cen fix at high-T):*/

  if (tempr > pow(10.e0,2.8e0)) {
    k[0] = 5.85e-11*tt*exp(-157809.1e0/tempr)*t5fac;
  } else {
    k[0]=0.e0;
  }

  if (tempr > pow(10.e0,3.e0)) {
    k[1] = 2.38e-11*tt*exp(-285335.4e0/tempr)*t5fac;
  } else {
    k[1] = 0.e0;
  }

  if (tempr > pow(10.e0,3.4e0)) {
    k[2] = 5.68e-12*tt*exp(-631515.e0/tempr)*t5fac;
  } else {
    k[2] = 0.e0;
  }

/*--- additional reactions from Z. Haiman: */

  if (tempr <= 4000.e0) {
    k[6] = 1.38e-23*pow(tempr,1.845e0);
  } else {
    if (tempr < pow(10.e0,5.e0)) {
      k[6] = -6.157e-17 + 3.255e-20*tempr - 4.152e-25*tempr*tempr;
    } else {
      k[6] = 0.e0;
    }
  }
  if (k[6] <= 0.e0) k[6]=0.e0;
  k[7] = 6.4e-10;

/**** Try the Abel et al. 1997 rate! */     
  Tgas=tempr;
  k[8]=1.429e-18*pow(Tgas,0.762)*pow(Tgas,0.1523*log10(Tgas))*
          pow(Tgas,-3.274e-2*pow(log10(Tgas),2.e0));

  k[9] = 1.3e-9;
  k[10] = 1.68e-8*pow(tempr/300.e0,-0.29e0);
  k[11] = 5.0e-6/tt;
  k[12] = 7.0e-7/tt;
  if (tempr > pow(10.e0,2.2e0)) {
    k[13] = (2.7e-8/(tempr*tt))*exp(-43000.e0/tempr)*t5fac;
  } else {
    k[13] = 0.e0;
  }
  if(tempr <= pow(10.e0,2.2e0)) {
    k[14]=0.e0;
    k[15]=0.e0;
  } else {
    x=log10(tempr/1.e4);
         
    ncr_1 = pow(10.e0,4.00e0 - 0.416e0*x - 0.327e0*x*x);
    ncr_2 = pow(10.e0,4.13e0 - 0.968e0*x + 0.119e0*x*x);
         
    kh_1=3.52e-9*exp(-43900.e0/tempr);
    kh_2=5.48e-9*exp(-53000.e0/tempr);

    if(tempr >= 7390.e0) {
      kl_1=6.11e-14*exp(-29300.e0/tempr);
    } else { 
      if (tempr > pow(10.e0,2.6e0)) {
        kl_1=2.67e-15*exp(-pow(6750.e0/tempr,2.e0));
      } else {
        kl_1 = 0.e0;
      }
    }
         
    if(tempr >= 7291.e0) {
      kl_2=5.22e-14*exp(-32200.e0/tempr);
    } else { 
      if (tempr > pow(10.e0,2.8e0)) {
        kl_2=3.17e-15*exp(-(4060.e0/tempr)-pow(7500.e0/tempr,2.e0));
      } else {
        kl_2 = 0.e0;
      }
    }
         
    k[14] = kh_1*pow(kl_1/kh_1,1.e0/(1.e0+yab[0]/ncr_1));
    k[15] = kh_2*pow(kl_2/kh_2,1.e0/(1.e0+yab[0]/ncr_2));
         
  }
      
/*--- More reactions from Shapiro and Kang:  */

  if (tempr > pow(10.e0,2.e0)) {
    k[16]=2.40e-9*exp(-21200.e0/tempr);
  } else {
    k[16] = 0.e0;
  }

  if (tempr > pow(10.e0,2.6e0)) {
    k[17]=4.38e-10*exp(-102000.e0/tempr)*pow(tempr,0.35e0);
  } else {
    k[17] = 0.e0;
  }

  if (tempr > pow(10.e0,1.6e0)) {
    k[18]=4.e-12*tempr*exp(-8750.e0/tempr);
  } else {
    k[18] = 0.e0;
  }

  if (tempr > pow(10.e0,1.6e0)) {
    k[19]=5.3e-20*pow(tempr,2.17e0)*exp(-8750.e0/tempr);
  } else {
    k[19] = 0.e0;
  }

  if (tempr > 1.e4) {
    k[20]=1.e-8*pow(tempr,-0.4e0);
  } else {
    if (tempr > pow(10.e0,1.8e0)) {
      k[20]=4.e-4*pow(tempr,-1.4e0)*exp(-15100.e0/tempr);
    } else {
      k[20] = 0.e0;
    }
  }

/*** Radiative processes
**      CALL photo
** Naoki: I still have to translate the photo subroutine into C */

/* JJ ionz195 */
/* JJ - Below, I've set photoionization on for whole box. */

/*  if (del_r < 5000.e0 && del_r >= 1. ){  
    k[21]=6.4e-6/pow(del_r,2.e0);
    k[22]=6.4e-6/pow(del_r,2.e0);
    k[23]=2.1e-7/pow(del_r,2.e0);
  } else if (del_r < 1.) {
    k[21]=3.6e-6;
    k[22]=6.4e-6;
    k[23]=2.1e-7;
  } else {
    k[21]=0.e0;
    k[22]=0.e0;
    k[23]=0.e0;
  }*/
    /* JJ end*/
  k[21]=0.e0;
  k[22]=0.e0;
  k[23]=0.e0;

/* set in P_chem for LW background */
  k[24]=0.e0;
  k[25]=0.e0;
  k[26]=0.e0;
  k[27]=0.e0;

/*** 3-body-reactions from Palla et al.,1983,ApJ 271,632  */

  k[28]=5.5e-29/tempr;
  k[29]=0.125e0*k[28];


  return;
}

/*   This function returns the H2 self-sheilding fraction
 *   f_shield = 1 -- no shielding
 *   f_shield = 0 -- complete shielding
 *   Need physical density in cgs 
 *  
 *   implemeted by JdJ 2016 see CSS 2012 for details
 */
double get_fshield(int target, double tempr, double xn, double rho)
{
  double f_shield,x,b,b5,alpha,N_H2,L_J;

  /* The H2 column density is estimated by N_H2 ~ L_J * n_H2  where L_J is the local Jeans length*/
  L_J = sqrt((15.*BOLTZMANN*tempr)/(4.*3.1459*6.67e-8*PROTONMASS*rho));
  N_H2 = SphP[target].Primordial_Chem[3] * xn * SphP[target].Primordial_Chem[0] * L_J;

  /* estimate for self-shielding from Draine & Bertoldi 1996, includes model for dynamic medium */
  x = N_H2 / 5.e14;
  b = 912000. * pow(tempr/1.e4,0.5);
  b5 = b / 1.e5;
  alpha = 1.1;  // modification from Wolcott-Green 2011 (orginal alpha = 2.0)

  f_shield = (0.965/pow(1.+x/b5,alpha)) + (0.035/pow(1.+x,0.5)) * exp(-8.5e-4*pow(1+x,0.5));

  return f_shield;
}

/*  this function computes the self-consistent temperature
 *  and abundance ratios 
 */
double AbundanceRatios(double u, double rho, double *ne_guess, double *nH0_pointer, double *nHeII_pointer, int target)
{
  double temp;

  DoCool_u_old_input = u;
  DoCool_rho_input = rho;
  DoCool_ne_guess_input = *ne_guess;

  rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;	/* convert to physical cgs units */
  u *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;

  temp = convert_u_to_temp(u, rho, ne_guess, target);

  *nH0_pointer = nH0;
  *nHeII_pointer = nHep;

  return temp;
}




extern FILE *fd;





/*  Calculates (heating rate-cooling rate)/n_h^2 in cgs units 
 */
double CoolingRate(double logT, double rho, double *nelec, int target)
{
  double Lambda, Heat;
  double LambdaExc, LambdaIon, LambdaRec, LambdaFF, LambdaCmptn = 0.0;
  double LambdaExcH0, LambdaExcHep, LambdaIonH0, LambdaIonHe0, LambdaIonHep;
  double LambdaRecHp, LambdaRecHep, LambdaRecHepp, LambdaRecHepd;
  double redshift;
  double T;
  double NH_SS_z=NH_SS,shieldfac;
#ifdef COOL_LOW_TEMPERATURES
  double LambdaMol=0;
#endif
#ifdef COOL_METAL_LINES_BY_SPECIES
  double LambdaMetal=0;
  double *Z;
  if(target>=0)
  {
      Z = P[target].Metallicity;
  } else {
      /* initialize dummy values here so the function doesn't crash, if called when there isn't a target particle */
      int k;
      double Zsol[NUM_METAL_SPECIES];
      for(k=0;k<NUM_METAL_SPECIES;k++) Zsol[k]=All.SolarAbundances[k];
      Z = Zsol;
  }
#endif
  double local_gammamultiplier=1;

  if(logT <= Tmin)
    logT = Tmin + 0.5 * deltaT;	/* floor at Tmin */


  nHcgs = XH * rho / PROTONMASS;	/* hydrogen number dens in cgs units */

#ifdef GALSF_FB_LOCAL_UV_HEATING
    double LambdaPElec,photoelec=0;
    if((target >= 0) && (gJH0 > 0))
    {
    local_gammamultiplier = SphP[target].RadFluxEUV * 2.29e-10; // converts to GammaHI for typical SED (rad_uv normalized to Habing)
    local_gammamultiplier = 1 + local_gammamultiplier/gJH0;
    }
    if(target >= 0) photoelec=SphP[target].RadFluxUV;
#endif
    
    /* CAFG: if density exceeds NH_SS, ignore ionizing background. */
    if(J_UV != 0)
        NH_SS_z=NH_SS*pow(local_gammamultiplier*gJH0/1.0e-12,0.66)*pow(10.,0.173*(logT-4.));
    else
        NH_SS_z=NH_SS*pow(10.,0.173*(logT-4.));
    if(nHcgs<100.*NH_SS_z) shieldfac=exp(-nHcgs/NH_SS_z); else shieldfac=0;
#ifdef GALSF_EFFECTIVE_EQS
    shieldfac = 1; // self-shielding is implicit in the sub-grid model already //
#endif
    
#ifdef BH_COMPTON_HEATING
    double AGN_LambdaPre,AGN_T_Compton;
    AGN_T_Compton = 2.0e7; /* approximate from Sazonov et al. */
    if(target < 0) {
        AGN_LambdaPre = 0;
    } else {
        AGN_LambdaPre = SphP[target].RadFluxAGN * (3.9/2.0) * All.UnitMass_in_g/(All.UnitLength_in_cm*All.UnitLength_in_cm)*All.HubbleParam*All.cf_a2inv; /* proper units */
        /* now have incident flux, need to convert to relevant pre-factor for heating rate */
        AGN_LambdaPre *= 6.652e-25; /* sigma_T for absorption */
        AGN_LambdaPre *= (4.*1.381e-16)/(9.109e-28*2.998e10*2.998e10); /* times 4*k_B/(me*c^2) */
    }
#endif

    
  T = pow(10.0, logT);
  if(logT < Tmax)
    {
      find_abundances_and_rates(logT, rho, nelec, target);
      /* Compute cooling and heating rate (cf KWH Table 1) in units of nH**2 */
        
      LambdaExcH0 = bH0 * ne * nH0;
      LambdaExcHep = bHep * ne * nHep;
      LambdaExc = LambdaExcH0 + LambdaExcHep;	/* excitation */

      LambdaIonH0 = 2.18e-11 * geH0 * ne * nH0;
      LambdaIonHe0 = 3.94e-11 * geHe0 * ne * nHe0;
      LambdaIonHep = 8.72e-11 * geHep * ne * nHep;
      LambdaIon = LambdaIonH0 + LambdaIonHe0 + LambdaIonHep;	/* ionization */

      LambdaRecHp = 1.036e-16 * T * ne * (aHp * nHp);
      LambdaRecHep = 1.036e-16 * T * ne * (aHep * nHep);
      LambdaRecHepp = 1.036e-16 * T * ne * (aHepp * nHepp);
      LambdaRecHepd = 6.526e-11 * ad * ne * nHep;
      LambdaRec = LambdaRecHp + LambdaRecHep + LambdaRecHepp + LambdaRecHepd;

      LambdaFF = bff * (nHp + nHep + 4 * nHepp) * ne;

      Lambda = LambdaExc + LambdaIon + LambdaRec + LambdaFF;

#ifdef COOL_METAL_LINES_BY_SPECIES
        //if((logT > Tmin+0.5*deltaT)&&((logT > 3.87)||(nHcgs<NH_SS_z)))
        /* can restrict to low-densities where not self-shielded, but let shieldfac (in ne) take care of this self-consistently */
        if((J_UV != 0)&&(logT > Tmin+0.5*deltaT)&&(logT > 4.00))
        {
            /* cooling rates tabulated for each species from Wiersma, Schaye, & Smith tables (2008) */
            LambdaMetal = GetCoolingRateWSpecies(nHcgs, logT, Z); //* nHcgs*nHcgs;
            /* tables normalized so ne*ni/(nH*nH) included already, so just multiply by nH^2 */
            /* (sorry, -- dont -- multiply by nH^2 here b/c that's how everything is normalized in this function) */
            LambdaMetal *= ne;
            /* (modified now to correct out tabulated ne so that calculated ne can be inserted;
             ni not used b/c it should vary species-to-species */
            Lambda += LambdaMetal;
        }
#endif

#ifdef COOL_LOW_TEMPERATURES
        //if((nHcgs>NH_SS_z)&&(logT <= 5.2)&&(logT > Tmin+0.5*deltaT))
        if((logT <= 5.2)&&(logT > Tmin+0.5*deltaT))
        {
            /* approx to cooling function for solar metallicity and nH=1 cm^(-3) -- want to do something
             much better, definitely, but for now use this just to get some idea of system with cooling to very low-temp */
            LambdaMol = 2.8958629e-26/(pow(T/125.21547,-4.9201887)+pow(T/1349.8649,-1.7287826)+pow(T/6450.0636,-0.30749082));//*nHcgs*nHcgs;
            LambdaMol *= (1-shieldfac);
            double LambdaDust = 0;
#ifdef COOL_METAL_LINES_BY_SPECIES
            LambdaMol *= (1+Z[0]/All.SolarAbundances[0])*(0.001 + 0.1*nHcgs/(1.0+nHcgs)
                            + 0.09*nHcgs/(1.0+0.1*nHcgs)
                            + (Z[0]/All.SolarAbundances[0])*(Z[0]/All.SolarAbundances[0])/(1.0+nHcgs));
            /* add dust cooling as well */
            double Tdust = 30.;
            if(T > Tdust) {LambdaDust = 0.63e-33 * (T-Tdust) * sqrt(T) * (Z[0]/All.SolarAbundances[0]);}
#endif
            Lambda += LambdaMol + LambdaDust;
            
        }
#endif
        
        
      if(All.ComovingIntegrationOn)
	{
	  redshift = 1 / All.Time - 1;
	  LambdaCmptn = 5.65e-36 * ne * (T - 2.73 * (1. + redshift)) * pow(1. + redshift, 4.) / nHcgs;

	  Lambda += LambdaCmptn;
	}
      else
	LambdaCmptn = 0;

#ifdef BH_COMPTON_HEATING
	if(T > AGN_T_Compton)
	{
        LambdaCmptn = AGN_LambdaPre * (T - AGN_T_Compton) * ne/nHcgs;
        if(LambdaCmptn > 2.19e-21/sqrt(T/1.0e8)) LambdaCmptn=2.19e-21/sqrt(T/1.0e8);
        Lambda += LambdaCmptn;
	}
#endif
        
      Heat = 0;
        if(J_UV != 0) {
            /* CAFG: if density exceeds NH_SS, ignore ionizing background. */
            Heat += local_gammamultiplier * (nH0 * epsH0 + nHe0 * epsHe0 + nHep * epsHep) / nHcgs * shieldfac;
        }
#if defined(COSMIC_RAYS) && !defined(COSMIC_RAYS_DISABLE_COOLING)
        if(SphP[target].CosmicRayEnergyPred > 0)
        {
            /* cosmic ray heating, from Guo & Oh 2008: this scales proportional to the electron number density and 
                cosmic ray energy density, both of which we quickly evaluate here (make sure we convert to the correct per-atom units) 
                - note that only 1/6 of the hadronic cooling is thermalized, according to their calculation, while all the Coulomb losses heat */
            double Gamma_CR = 1.0e-16 * (0.98 + 1.65*ne*XH) / nHcgs *
                ((SphP[target].CosmicRayEnergyPred / P[target].Mass * SphP[target].Density * All.cf_a3inv) *
                 (All.UnitPressure_in_cgs * All.HubbleParam * All.HubbleParam));
            Heat += Gamma_CR;
        }
#endif
        
#ifdef COOL_LOW_TEMPERATURES
#if !defined(COSMIC_RAYS) || defined(COSMIC_RAYS_DISABLE_COOLING)
        /* if COSMIC_RAYS is not enabled, but low-temperature cooling is on, we account for the CRs as a heating source using
         a more approximate expression (assuming the mean background of the Milky Way clouds) */
        if(logT <= 5.2)
        {
            double Gamma_CR = 1.0e-16 * (0.98 + 1.65*ne*XH) / (1.e-2 + nHcgs) * 9.0e-12;
            // multiplied by background of ~5eV/cm^3 (Goldsmith & Langer (1978),  van Dishoeck & Black (1986) //
            Heat += Gamma_CR;
        }
#endif
#ifdef COOL_METAL_LINES_BY_SPECIES
        /* add dust heating as well */
        double Tdust = 30.;
        if(T < Tdust) {Heat += 0.63e-33 * (Tdust-T) * sqrt(Tdust) * (Z[0]/All.SolarAbundances[0]);}
#endif
#endif
        
        
        
#ifdef BH_COMPTON_HEATING
        if(T < AGN_T_Compton) Heat += AGN_LambdaPre * (AGN_T_Compton - T) / nHcgs;
        /* note this is independent of the free electron fraction */
#endif
#ifdef GALSF_FB_LOCAL_UV_HEATING
        /* photoelectric heating following Bakes & Thielens 1994 (also Wolfire 1995) */
        if(T < 1.0e6) {
            LambdaPElec = 1.0e-24*photoelec/nHcgs;
            photoelec *= sqrt(T)/nHcgs;
            LambdaPElec *= 0.049/(1+pow(photoelec/1925.,0.73)) + 0.037*pow(T/1.0e4,0.7)/(1+photoelec/5000.);
            Heat += LambdaPElec;
        }
#endif
    }
  else				/* here we're outside of tabulated rates, T>Tmax K */
    {
      /* at high T (fully ionized); only free-free and Compton cooling are present.  
         Assumes no heating. */

      Heat = 0;

      LambdaExcH0 = LambdaExcHep = LambdaIonH0 = LambdaIonHe0 = LambdaIonHep =
	LambdaRecHp = LambdaRecHep = LambdaRecHepp = LambdaRecHepd = 0;

      /* very hot: H and He both fully ionized */
      nHp = 1.0;
      nHep = 0;
      nHepp = yhelium;
      ne = nHp + 2.0 * nHepp;
      *nelec = ne;		/* note: in units of the hydrogen number density */

      LambdaFF = 1.42e-27 * sqrt(T) * (1.1 + 0.34 * exp(-(5.5 - logT) * (5.5 - logT) / 3)) * (nHp + 4 * nHepp) * ne;

      if(All.ComovingIntegrationOn)
	{
	  redshift = 1 / All.Time - 1;
	  /* add inverse Compton cooling off the microwave background */
	  LambdaCmptn = 5.65e-36 * ne * (T - 2.73 * (1. + redshift)) * pow(1. + redshift, 4.) / nHcgs;
	}
      else
	LambdaCmptn = 0;

#ifdef BH_COMPTON_HEATING
        //LambdaCmptn += AGN_LambdaPre * (T - AGN_T_Compton) * ne/nHcgs;
        /* actually at these temperatures want approximation to relativistic compton cooling */
        LambdaCmptn += AGN_LambdaPre * (T - AGN_T_Compton) * (T/1.5e9)/(1-exp(-T/1.5e9)) * ne/nHcgs;
#endif
        
      Lambda = LambdaFF + LambdaCmptn;

      /* per CAFG's calculations, we should note that at very high temperatures, the rate-limiting step may be
         the Coulomb collisions moving energy from protons to e-; which if slow will prevent efficient e- cooling */
      if(Lambda > 2.19e-21/sqrt(T/1.0e8)) Lambda=2.19e-21/sqrt(T/1.0e8);
    }

  /*      
     printf("Lambda= %g\n", Lambda);

     fprintf(fd,"%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", pow(10, logT),Lambda,
     LambdaExcH0, LambdaExcHep, 
     LambdaIonH0, LambdaIonHe0, LambdaIonHep,
     LambdaRecHp, LambdaRecHep, LambdaRecHepp, LambdaRecHepd,
     LambdaFF, LambdaCmptn, Heat,
     ne, nHp, nHep, nHepp);
   */
    
    double Q = Heat - Lambda;
#ifdef COOL_LOW_TEMPERATURES
    /* if we are in the optically thick limit, we need to modify the cooling/heating rates according to the appropriate limits; 
        this flag does so by using a simple approximation. we consider the element as if it were a slab, with a column density 
        calculated from the simulation properties and the Sobolev approximation. we then assume it develops an equilibrium internal 
        temperature structure on a radiative diffusion timescale much faster than the dynamical time, and so the surface radiation 
        from a photosphere can be simply related to the local density by the optical depth to infinity. the equations here follow 
        Rafikov, 2007 (ApJ, 662, 642): 
            denergy/dt/dArea = sigma*T^4 / fc(tau)
            fc(tau) = tau^eta + 1/tau (taking chi, phi~1; the second term describes the optically thin limit, which is calculated above 
                more accurately anyways - that was just Kirchoff's Law; so we only need to worry about the first term)
            eta = 4*(gamma-1) / [gamma*(1+alpha+beta*(gamma-1)/gamma)], where gamma=real polytropic index, and alpha/beta follow
                an opacity law kappa=kappa_0 * P^alpha * T^beta. for almost all the regimes of interest, however, eta~1, which is also 
                what is obtained for a convectively stable slab. so we will use this.
            now, this gives sigma*T^4/tau * Area_eff / nHcgs as the 'effective' cooling rate in our units of Heat or Lambda above. 
                the nHcgs just puts it in the same volumetric terms. The Area_eff must be defined as ~m_particle/surface_density
                to have the same meaning for a slab as assumed in Rafikov (and to integrate correctly over all particles in the slab, 
                if/when the slab is resolved). We estimate this in our usual fashion with the Sobolev-type column density
            tau = kappa * surface_density; we estimate kappa ~ 5 cm^2/g * (0.001+Z/Z_solar), as the frequency-integrated kappa for warm 
                dust radiation (~150K), weighted by the dust-to-gas ratio (with a floor for molecular absorption). we could make this 
                temperature-dependent, though, fairly easily - for this particular problem it won't make much difference
        This rate then acts as an upper limit to the net heating/cooling calculated above (restricts absolute value)
     */
    //if(nHcgs > 0.1) /* don't bother at very low densities, since youre not optically thick */  
    if( (nHcgs > 0.1) && (target >= 0) )  // DAA: protect from target=-1 with GALSF_EFFECTIVE_EQS
    {
        double surface_density = evaluate_NH_from_GradRho(SphP[target].Gradients.Density,PPP[target].Hsml,SphP[target].Density,PPP[target].NumNgb,1);
        surface_density *= All.cf_a2inv * All.UnitDensity_in_cgs * All.HubbleParam * All.UnitLength_in_cm; // converts to cgs
        double effective_area = 2.3 * PROTONMASS / surface_density; // since cooling rate is ultimately per-particle, need a particle-weight here
        double kappa_eff; // effective kappa, accounting for metal abundance, temperature, and density //
        if(T < 1500.)
        {
            if(T < 150.) {kappa_eff=0.0027*T*sqrt(T);} else {kappa_eff=5.;}
            kappa_eff *= P[target].Metallicity[0]/All.SolarAbundances[0];
            if(kappa_eff < 0.1) {kappa_eff=0.1;}
        } else {
            /* this is an approximate result for high-temperature opacities, but provides a pretty good fit from 1.5e3 - 1.0e9 K */
            double k_electron = 0.2 * (1. + HYDROGEN_MASSFRAC); //0.167 * ne; /* Thompson scattering (non-relativistic) */
            double k_molecular = 0.1 * P[target].Metallicity[0]; /* molecular line opacities */
            double k_Hminus = 1.1e-25 * sqrt(P[target].Metallicity[0] * rho) * pow(T,7.7); /* negative H- ion opacity */
            double k_Kramers = 4.0e25 * (1.+HYDROGEN_MASSFRAC) * (P[target].Metallicity[0]+0.001) * rho / (T*T*T*sqrt(T)); /* free-free, bound-free, bound-bound transitions */
            double k_radiative = k_molecular + 1./(1./k_Hminus + 1./(k_electron+k_Kramers)); /* approximate interpolation between the above opacities */
            double k_conductive = 2.6e-7 * ne * T*T/(rho*rho); //*(1+pow(rho/1.e6,0.67) /* e- thermal conductivity can dominate at low-T, high-rho, here it as expressed as opacity */
            kappa_eff = 1./(1./k_radiative + 1./k_conductive); /* effective opacity including both heat carriers (this is exact) */
        }
        double tau_eff = kappa_eff * surface_density;
        double Lambda_Thick_BlackBody = 5.67e-5 * (T*T*T*T) * effective_area / ((1.+tau_eff) * nHcgs);
        if(Q > 0) {if(Q > Lambda_Thick_BlackBody) {Q=Lambda_Thick_BlackBody;}} else {if(Q < -Lambda_Thick_BlackBody) {Q=-Lambda_Thick_BlackBody;}}
    }
#endif
    
#ifndef COOLING_OPERATOR_SPLIT
    /* add the hydro energy change directly: this represents an additional heating/cooling term, to be accounted for 
        in the semi-implicit solution determined here. this is more accurate when tcool << tdynamical */
    if(target >= 0) Q += SphP[target].DtInternalEnergy / nHcgs;
#endif

  return Q;
}





double LogTemp(double u, double ne, int target)	/* ne= electron density in terms of hydrogen density */
{
  double temp, mu, meanweight;

  // if(u < ethmin)
  //   u = ethmin;

  // T = log10(GAMMA_MINUS1 * u * mhboltz * (1 + 4 * yhelium) / (1 + ne + yhelium));
  mu = (1 + 4 * yhelium) / (1 + yhelium + SphP[target].Primordial_Chem[5]);
  meanweight = mu * PROTONMASS;

  //temp = GAMMA_MINUS1 / BOLTZMANN * u * PROTONMASS * mu;
  temp = GAMMA_MINUS1 / BOLTZMANN * u * meanweight;
  
  if(temp<=0) temp=pow(10.0,Tmin);
  if(log10(temp)<Tmin) temp=pow(10.0,Tmin);
  if(log10(temp)>Tmax) temp=pow(10.0,Tmax);

  return log10(temp);
}



void InitCoolMemory(void)
{
  BetaH0 = (double *) mymalloc("BetaH0", (NCOOLTAB + 1) * sizeof(double));
  BetaHep = (double *) mymalloc("BetaHep", (NCOOLTAB + 1) * sizeof(double));
  AlphaHp = (double *) mymalloc("AlphaHp", (NCOOLTAB + 1) * sizeof(double));
  AlphaHep = (double *) mymalloc("AlphaHep", (NCOOLTAB + 1) * sizeof(double));
  Alphad = (double *) mymalloc("Alphad", (NCOOLTAB + 1) * sizeof(double));
  AlphaHepp = (double *) mymalloc("AlphaHepp", (NCOOLTAB + 1) * sizeof(double));
  GammaeH0 = (double *) mymalloc("GammaeH0", (NCOOLTAB + 1) * sizeof(double));
  GammaeHe0 = (double *) mymalloc("GammaeHe0", (NCOOLTAB + 1) * sizeof(double));
  GammaeHep = (double *) mymalloc("GammaeHep", (NCOOLTAB + 1) * sizeof(double));
  Betaff = (double *) mymalloc("Betaff", (NCOOLTAB + 1) * sizeof(double));

#ifdef COOL_METAL_LINES_BY_SPECIES
  long i_nH=41; long i_T=176; long kspecies=(long)NUM_METAL_SPECIES-1;
#ifdef GALSF_FB_RPROCESS_ENRICHMENT
    //kspecies -= 1;
    kspecies -= NUM_RPROCESS_SPECIES;
#endif
  SpCoolTable0 = (float *) mymalloc("SpCoolTable0",(kspecies*i_nH*i_T)*sizeof(float));
  if(All.ComovingIntegrationOn)
    SpCoolTable1 = (float *) mymalloc("SpCoolTable1",(kspecies*i_nH*i_T)*sizeof(float));
#endif
}



void MakeCoolingTable(void)
     /* Set up interpolation tables in T for cooling rates given in KWH, ApJS, 105, 19 
        Hydrogen, Helium III recombination rates and collisional ionization cross-sections are updated */
{
    int i;
    double T,Tfact;
    
    if(All.MinGasTemp > 0.0)
        Tmin = log10(All.MinGasTemp); // Tmin = log10(0.1 * All.MinGasTemp);
    else
        Tmin = 1.0;
    
    deltaT = (Tmax - Tmin) / NCOOLTAB;
    ethmin = pow(10.0, Tmin) * (1. + yhelium) / ((1. + 4. * yhelium) * mhboltz * GAMMA_MINUS1);
    /* minimum internal energy for neutral gas */
    for(i = 0; i <= NCOOLTAB; i++)
    {
        BetaH0[i] = BetaHep[i] = Betaff[i] = AlphaHp[i] = AlphaHep[i] = AlphaHepp[i] = Alphad[i] = GammaeH0[i] = GammaeHe0[i] = GammaeHep[i] = 0;
        T = pow(10.0, Tmin + deltaT * i);
        Tfact = 1.0 / (1 + sqrt(T / 1.0e5));
        
        if(118348 / T < 70) BetaH0[i] = 7.5e-19 * exp(-118348 / T) * Tfact;
        if(473638 / T < 70) BetaHep[i] = 5.54e-17 * pow(T, -0.397) * exp(-473638 / T) * Tfact;
        
        Betaff[i] = 1.43e-27 * sqrt(T) * (1.1 + 0.34 * exp(-(5.5 - log10(T)) * (5.5 - log10(T)) / 3));
        //AlphaHp[i] = 8.4e-11 * pow(T / 1000, -0.2) / (1. + pow(T / 1.0e6, 0.7)) / sqrt(T);	/* old Cen92 fit */
        //AlphaHep[i] = 1.5e-10 * pow(T, -0.6353); /* old Cen92 fit */
        //AlphaHepp[i] = 4. * AlphaHp[i];	/* old Cen92 fit */
        AlphaHp[i] = 7.982e-11 / ( sqrt(T/3.148) * pow((1.0+sqrt(T/3.148)), 0.252) * pow((1.0+sqrt(T/7.036e5)), 1.748) ); /* Verner & Ferland (1996) [more accurate than Cen92] */
        AlphaHep[i]= 9.356e-10 / ( sqrt(T/4.266e-2) * pow((1.0+sqrt(T/4.266e-2)), 0.2108) * pow((1.0+sqrt(T/3.676e7)), 1.7892) ); /* Verner & Ferland (1996) [more accurate than Cen92] */
        AlphaHepp[i] = 2. * 7.982e-11 / ( sqrt(T/(4.*3.148)) * pow((1.0+sqrt(T/(4.*3.148))), 0.252) * pow((1.0+sqrt(T/(4.*7.036e5))), 1.748) ); /* Verner & Ferland (1996) : ~ Z*alphaHp[1,T/Z^2] */
        
        if(470000 / T < 70) Alphad[i] = 1.9e-3 * pow(T, -1.5) * exp(-470000 / T) * (1. + 0.3 * exp(-94000 / T));
        if(157809.1 / T < 70) GammaeH0[i] = 5.85e-11 * sqrt(T) * exp(-157809.1 / T) * Tfact;
        if(285335.4 / T < 70) GammaeHe0[i] = 2.38e-11 * sqrt(T) * exp(-285335.4 / T) * Tfact;
        if(631515.0 / T < 70) GammaeHep[i] = 5.68e-12 * sqrt(T) * exp(-631515.0 / T) * Tfact;
        
    }
}


#ifdef COOL_METAL_LINES_BY_SPECIES

void LoadMultiSpeciesTables(void)
{
    if(All.ComovingIntegrationOn) {
        int i;
        double z;
        if(All.Time==All.TimeBegin) {
            All.SpeciesTableInUse=48;
            ReadMultiSpeciesTables(All.SpeciesTableInUse);
        }
        z=log10(1/All.Time)*48;
        i=(int)z;
        if(i<48) {
            if(i<All.SpeciesTableInUse) {
                All.SpeciesTableInUse=i;
                ReadMultiSpeciesTables(All.SpeciesTableInUse);
            }}
    } else {
        if(All.Time==All.TimeBegin) ReadMultiSpeciesTables(0);
    }
}

void ReadMultiSpeciesTables(int iT)
{
    /* read table w n,T for each species */
    long i_nH=41; long i_Temp=176; long kspecies=(long)NUM_METAL_SPECIES-1; long i,j,k,r;
#ifdef GALSF_FB_RPROCESS_ENRICHMENT
    //kspecies -= 1;
    kspecies -= NUM_RPROCESS_SPECIES;
#endif
    /* int i_He=7;  int l; */
    FILE *fdcool; char *fname;
    
    fname=GetMultiSpeciesFilename(iT,0);
    if(ThisTask == 0) printf("Opening Cooling Table %s \n",fname);
    if(!(fdcool = fopen(fname, "r"))) {
        printf(" Cannot read species cooling table in file `%s'\n", fname); endrun(456);}
    for(i=0;i<kspecies;i++) {
        for(j=0;j<i_nH;j++) {
            for(k=0;k<i_Temp;k++) {
                r=fread(&SpCoolTable0[i*i_nH*i_Temp + j*i_Temp + k],sizeof(float),1,fdcool);
                if(r!=1) {printf(" Reached Cooling EOF! \n");fflush(stdout);}
            }}}
    fclose(fdcool);
    /*
     GetMultiSpeciesFilename(iT,&fname,1);
     if(!(fdcool = fopen(fname, "r"))) {
     printf(" Cannot read species (He) cooling table in file `%s'\n", fname); endrun(456);}
     for(i=0;i<2;i++)
     for(j=0;j<i_nH;j++)
     for(k=0;k<i_Temp;k++)
     for(l=0;l<i_He;l++)
     fread(&SpCoolTable0_He[i][j][k][l],sizeof(float),1,fdcool);
     fclose(fdcool);
     */
    if (All.ComovingIntegrationOn && i<48) {
        fname=GetMultiSpeciesFilename(iT+1,0);
        if(ThisTask == 0) printf("Opening (z+) Cooling Table %s \n",fname);
        if(!(fdcool = fopen(fname, "r"))) {
            printf(" Cannot read species 1 cooling table in file `%s'\n", fname); endrun(456);}
        for(i=0;i<kspecies;i++) {
            for(j=0;j<i_nH;j++) {
                for(k=0;k<i_Temp;k++) {
                    r=fread(&SpCoolTable1[i*i_nH*i_Temp + j*i_Temp + k],sizeof(float),1,fdcool);
                    if(r!=1) {printf(" Reached Cooling EOF! \n");fflush(stdout);}
                }}}
        fclose(fdcool);
        /*
         GetMultiSpeciesFilename(iT+1,&fname,1);
         if(!(fdcool = fopen(fname, "r"))) {
         printf(" Cannot read species 1 (He) cooling table in file `%s'\n", fname); endrun(456);}
         for(i=0;i<2;i++)
         for(j=0;j<i_nH;j++)
         for(k=0;k<i_Temp;k++)
         for(l=0;l<i_He;l++)
         fread(&SpCoolTable1_He[i][j][k][l],sizeof(float),1,fdcool);
         fclose(fdcool);
         */
    }
}

char *GetMultiSpeciesFilename(int i, int hk)
{
    static char fname[100];
    if(i<0) i=0; if(i>48) i=48;
    if(hk==0) {
        sprintf(fname,"./spcool_tables/spcool_%d",i);
    } else {
        sprintf(fname,"./spcool_tables/spcool_He_%d",i);
    }
    return fname;
}

#endif



/* table input (from file TREECOOL) for ionizing parameters */
/* NOTE: we've switched to using the updated TREECOOL from CAFG, june11 version */

// #define JAMPL	1.0		/* amplitude factor relative to input table */
// #define TABLESIZE 250		/* Max # of lines in TREECOOL */

// static float inlogz[TABLESIZE];
// static float gH0[TABLESIZE], gHe[TABLESIZE], gHep[TABLESIZE];
// static float eH0[TABLESIZE], eHe[TABLESIZE], eHep[TABLESIZE];
// static int nheattab;		/* length of table */


void ReadIonizeParams(char *fname)
{
  int i;
  FILE *fdcool;

  if(!(fdcool = fopen(fname, "r")))
    {
      printf(" Cannot read ionization table in file `%s'\n", fname);
      endrun(456);
    }

  for(i = 0; i < TABLESIZE; i++)
    gH0[i] = 0;

  for(i = 0; i < TABLESIZE; i++)
    if(fscanf(fdcool, "%g %g %g %g %g %g %g",
	      &inlogz[i], &gH0[i], &gHe[i], &gHep[i], &eH0[i], &eHe[i], &eHep[i]) == EOF)
      break;

  fclose(fdcool);

  /*  nheattab is the number of entries in the table */

  for(i = 0, nheattab = 0; i < TABLESIZE; i++)
    if(gH0[i] != 0.0)
      nheattab++;
    else
      break;

  if(ThisTask == 0)
    printf("\n\nread ionization table with %d entries in file `%s'.\n\n", nheattab, fname);
}


void IonizeParams(void)
{
  IonizeParamsTable();

  /*
     IonizeParamsFunction();
   */
}



void IonizeParamsTable(void)
{
  int i, ilow;
  double logz, dzlow, dzhi;
  double redshift;

  if(All.ComovingIntegrationOn)
    redshift = 1 / All.Time - 1;
  else
    {
    /* in non-cosmological mode, still use, but adopt z=0 background */
    redshift = 0;
    /*
         gJHe0 = gJHep = gJH0 = 0;
         epsHe0 = epsHep = epsH0 = 0;
         J_UV = 0;
         return;
    */
    }

  logz = log10(redshift + 1.0);
  ilow = 0;
  for(i = 0; i < nheattab; i++)
    {
      if(inlogz[i] < logz)
	ilow = i;
      else
	break;
    }

  dzlow = logz - inlogz[ilow];
  dzhi = inlogz[ilow + 1] - logz;

  if(logz > inlogz[nheattab - 1] || gH0[ilow] == 0 || gH0[ilow + 1] == 0 || nheattab == 0)
    {
      gJHe0 = gJHep = gJH0 = 0;
      epsHe0 = epsHep = epsH0 = 0;
      J_UV = 0;
      return;
    }
  else
    J_UV = 1.e-21;		/* irrelevant as long as it's not 0 */

  gJH0 = JAMPL * pow(10., (dzhi * log10(gH0[ilow]) + dzlow * log10(gH0[ilow + 1])) / (dzlow + dzhi));
  gJHe0 = JAMPL * pow(10., (dzhi * log10(gHe[ilow]) + dzlow * log10(gHe[ilow + 1])) / (dzlow + dzhi));
  gJHep = JAMPL * pow(10., (dzhi * log10(gHep[ilow]) + dzlow * log10(gHep[ilow + 1])) / (dzlow + dzhi));
  epsH0 = JAMPL * pow(10., (dzhi * log10(eH0[ilow]) + dzlow * log10(eH0[ilow + 1])) / (dzlow + dzhi));
  epsHe0 = JAMPL * pow(10., (dzhi * log10(eHe[ilow]) + dzlow * log10(eHe[ilow + 1])) / (dzlow + dzhi));
  epsHep = JAMPL * pow(10., (dzhi * log10(eHep[ilow]) + dzlow * log10(eHep[ilow + 1])) / (dzlow + dzhi));

  return;
}


void SetZeroIonization(void)
{
  gJHe0 = gJHep = gJH0 = 0;
  epsHe0 = epsHep = epsH0 = 0;
  J_UV = 0;
}


void IonizeParamsFunction(void)
{
  int i, nint;
  double a0, planck, ev, e0_H, e0_He, e0_Hep;
  double gint, eint, t, tinv, fac, eps;
  double at, beta, s;
  double pi;

#define UVALPHA         1.0
  double Jold = -1.0;
  double redshift;

  J_UV = 0.;
  gJHe0 = gJHep = gJH0 = 0.;
  epsHe0 = epsHep = epsH0 = 0.;


  if(All.ComovingIntegrationOn)	/* analytically compute params from power law J_nu */
    {
      redshift = 1 / All.Time - 1;

      if(redshift >= 6)
	J_UV = 0.;
      else
	{
	  if(redshift >= 3)
	    J_UV = 4e-22 / (1 + redshift);
	  else
	    {
	      if(redshift >= 2)
		J_UV = 1e-22;
	      else
		J_UV = 1.e-22 * pow(3.0 / (1 + redshift), -3.0);
	    }
	}

      if(J_UV == Jold)
	return;


      Jold = J_UV;

      if(J_UV == 0)
	return;


      a0 = 6.30e-18;
      planck = 6.6262e-27;
      ev = 1.6022e-12;
      e0_H = 13.6058 * ev;
      e0_He = 24.59 * ev;
      e0_Hep = 54.4232 * ev;

      gint = 0.0;
      eint = 0.0;
      nint = 5000;
      at = 1. / ((double) nint);

      for(i = 1; i <= nint; i++)
	{
	  t = (double) i;
	  t = (t - 0.5) * at;
	  tinv = 1. / t;
	  eps = sqrt(tinv - 1.);
	  fac = exp(4. - 4. * atan(eps) / eps) / (1. - exp(-2. * M_PI / eps)) * pow(t, UVALPHA + 3.);
	  gint += fac * at;
	  eint += fac * (tinv - 1.) * at;
	}

      gJH0 = a0 * gint / planck;
      epsH0 = a0 * eint * (e0_H / planck);
      gJHep = gJH0 * pow(e0_H / e0_Hep, UVALPHA) / 4.0;
      epsHep = epsH0 * pow((e0_H / e0_Hep), UVALPHA - 1.) / 4.0;

      at = 7.83e-18;
      beta = 1.66;
      s = 2.05;

      gJHe0 = (at / planck) * pow((e0_H / e0_He), UVALPHA) *
	(beta / (UVALPHA + s) + (1. - beta) / (UVALPHA + s + 1));
      epsHe0 = (e0_He / planck) * at * pow(e0_H / e0_He, UVALPHA) *
	(beta / (UVALPHA + s - 1) + (1 - 2 * beta) / (UVALPHA + s) - (1 - beta) / (UVALPHA + s + 1));

      pi = M_PI;
      gJH0 *= 4. * pi * J_UV;
      gJHep *= 4. * pi * J_UV;
      gJHe0 *= 4. * pi * J_UV;
      epsH0 *= 4. * pi * J_UV;
      epsHep *= 4. * pi * J_UV;
      epsHe0 *= 4. * pi * J_UV;
    }
}





void InitCool(void)
{
    if(ThisTask == 0)
        printf("Initializing cooling ...\n");
    
    All.Time = All.TimeBegin;
    set_cosmo_factors_for_current_time();
    
#ifdef GRACKLE
    InitGrackle();
#endif

    InitCoolMemory();
    MakeCoolingTable();
    ReadIonizeParams("TREECOOL");
    IonizeParams();
#ifdef COOL_METAL_LINES_BY_SPECIES
    LoadMultiSpeciesTables();
#endif
}




#ifdef COOL_METAL_LINES_BY_SPECIES
double GetCoolingRateWSpecies(double nHcgs, double logT, double *Z)
{
    int k;
    double ne_over_nh_tbl=1, Lambda=0;
    int N_species_active = NUM_METAL_SPECIES-1;
#ifdef GALSF_FB_RPROCESS_ENRICHMENT
    //N_species_active -= 1;
    N_species_active -= NUM_RPROCESS_SPECIES;
#endif
    
    ne_over_nh_tbl = GetLambdaSpecies(0,nHcgs,logT,0);
    for (k=1;k<N_species_active;k++)
        Lambda += GetLambdaSpecies(k,nHcgs,logT,0) * Z[k+1]/(All.SolarAbundances[k+1]*0.0127/All.SolarAbundances[0]);
    
    if(ne_over_nh_tbl>0) Lambda /= ne_over_nh_tbl; else Lambda=0;
    return Lambda;
}

double getSpCoolTableVal(long i,long j,long k,long tblK)
{
    long i_nH=41; long i_T=176; long inHT=i_nH*i_T;
    if(tblK==0) {
        //return SpCoolTable0[i][j][k];
        return SpCoolTable0[i*inHT + j*i_T + k];
    }
    if(tblK==1) {
        //return SpCoolTable1[i][j][k];
        return SpCoolTable1[i*inHT + j*i_T + k];
    }
    return 0;
}

double GetLambdaSpecies(int k_species, double nHcgs, double logT, double fHe)
{
    int ix0,iy0,ix1,iy1,ixmax,iymax;/*,ih0,ih1,ihmax;*/
    double i1,i2,j1,j2,w1,w2,u1;
    double v000,v010,v100,v110,v001,v011,v101,v111;
    double dx,dy,dz,mdz;/*,dh,mdh;*/
    dx=dy=dz=0; ix0=ix1=iy0=iy1=0; ixmax=40; iymax=175; /*ihmax=6;*/
    
    dx = (log10(nHcgs)-(-8.0))/(0.0-(-8.0))*ixmax;
    if(dx<0) dx=0; if(dx>ixmax) dx=ixmax;
    ix0=(int)dx; ix1=ix0+1; if(ix1>ixmax) ix1=ixmax;
    dx=dx-ix0;
    dy = (logT-2.0)/(9.0-2.0)*iymax;
    if(dy<0) dy=0; if(dy>iymax) dy=iymax;
    iy0=(int)dy; iy1=iy0+1; if(iy1>iymax) iy1=iymax;
    dy=dy-iy0;
    /*
     if(k_species<=1) {
     dh = (fHe-0.238)/(0.298-0.238)*ihmax;
     if(dh<0) dh=0; if(dh>ihmax) dh=ihmax;
     ih0=(int)dh; ih1=ih0+1; if(ih1>ihmax) ih1=ihmax;
     dh=dh-ih0;
     } */
    
    v000=getSpCoolTableVal(k_species,ix0,iy0,0);
    v010=getSpCoolTableVal(k_species,ix0,iy1,0);
    v100=getSpCoolTableVal(k_species,ix1,iy0,0);
    v110=getSpCoolTableVal(k_species,ix1,iy1,0);
    /*if(k_species>1) {
     v000=SpCoolTable0[k_species-2][ix0][iy0];
     v010=SpCoolTable0[k_species-2][ix0][iy1];
     v100=SpCoolTable0[k_species-2][ix1][iy0];
     v110=SpCoolTable0[k_species-2][ix1][iy1];
     } else {
     v000=SpCoolTable0_He[k_species][ix0][iy0][ih0]*mdh + SpCoolTable0_He[k_species][ix0][iy0][ih1]*dh;
     v010=SpCoolTable0_He[k_species][ix0][iy1][ih0]*mdh + SpCoolTable0_He[k_species][ix0][iy1][ih1]*dh;
     v100=SpCoolTable0_He[k_species][ix1][iy0][ih0]*mdh + SpCoolTable0_He[k_species][ix1][iy0][ih1]*dh;
     v110=SpCoolTable0_He[k_species][ix1][iy1][ih0]*mdh + SpCoolTable0_He[k_species][ix1][iy1][ih1]*dh;
     }*/
    if(All.ComovingIntegrationOn) {
        if(All.SpeciesTableInUse<48) {
            dz=log10(1/All.Time)*48; dz=dz-(int)dz;
            v001=getSpCoolTableVal(k_species,ix0,iy0,1);
            v011=getSpCoolTableVal(k_species,ix0,iy1,1);
            v101=getSpCoolTableVal(k_species,ix1,iy0,1);
            v111=getSpCoolTableVal(k_species,ix1,iy1,1);
            /*if(k_species>1) {
             v001=SpCoolTable1[k_species-2][ix0][iy0];
             v011=SpCoolTable1[k_species-2][ix0][iy1];
             v101=SpCoolTable1[k_species-2][ix1][iy0];
             v111=SpCoolTable1[k_species-2][ix1][iy1];
             } else {
             v001=SpCoolTable1_He[k_species][ix0][iy0][ih0]*mdh + SpCoolTable1_He[k_species][ix0][iy0][ih1]*dh;
             v011=SpCoolTable1_He[k_species][ix0][iy1][ih0]*mdh + SpCoolTable1_He[k_species][ix0][iy1][ih1]*dh;
             v101=SpCoolTable1_He[k_species][ix1][iy0][ih0]*mdh + SpCoolTable1_He[k_species][ix1][iy0][ih1]*dh;
             v111=SpCoolTable1_He[k_species][ix1][iy1][ih0]*mdh + SpCoolTable1_He[k_species][ix1][iy1][ih1]*dh;
             }*/
        } else {
            v001=v011=v101=v111=dz=0;
        }
    } else {
        v001=v011=v101=v111=dz=0;
    }
    mdz=1-dz;
    i1=v000*mdz + v001*dz;
    i2=v010*mdz + v011*dz;
    j1=v100*mdz + v101*dz;
    j2=v110*mdz + v111*dz;
    w1=i1*(1-dy) + i2*dy;
    w2=j1*(1-dy) + j2*dy;
    u1=w1*(1-dx) + w2*dx;
    return u1;
}

#endif // COOL_METAL_LINES_BY_SPECIES



#ifdef GALSF_FB_LOCAL_UV_HEATING
void selfshield_local_incident_uv_flux(void)
{
    /* include local self-shielding with the following */
    int i;
    double GradRho,sigma_eff_0;
    
    sigma_eff_0 = 0.955*All.UnitMass_in_g*All.HubbleParam / (All.UnitLength_in_cm*All.UnitLength_in_cm);
    if(All.ComovingIntegrationOn) sigma_eff_0 /= All.Time*All.Time;
    
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].Type==0)
        {
            if((SphP[i].RadFluxUV>0) && (PPP[i].Hsml>0) && (SphP[i].Density>0) && (P[i].Mass>0) && (All.Time>0))
            {
                GradRho = sigma_eff_0 * evaluate_NH_from_GradRho(P[i].GradRho,PPP[i].Hsml,SphP[i].Density,PPP[i].NumNgb,1);
                SphP[i].RadFluxUV *= 1276.19 * sigma_eff_0 * exp(-KAPPA_UV*GradRho);
                GradRho *= 3.7e6; // 912 angstrom KAPPA_EUV //
                //SphP[i].RadFluxEUV *= 1276.19 * sigma_eff_0 * exp(-GradRho);
                SphP[i].RadFluxEUV *= 1276.19 * sigma_eff_0 * (0.01 + 0.99/(1.0+0.8*GradRho+0.85*GradRho*GradRho));
                // unit conversion and self-shielding (normalized to Habing (local MW) field)
            } else {
                SphP[i].RadFluxUV = 0;
                SphP[i].RadFluxEUV = 0;
            }}}
}
#endif // GALSF_FB_LOCAL_UV_HEATING

/* Params for computing cold gas fractions */
#define P0BLITZ 1.7e+04 // Leroy et al 2008, Fig 17 (THINGS) Table 6
#define ALPHA0BLITZ 0.8 // Leroy et al 2008, Fig 17 (THINGS) Table 6
#define NHILIM 1.73e18    // gives tau~1/fion(IGM)
#define FSHIELD 0.99    // maximum fraction of mass that is self-shielded

//#define KG11_H2	// if fH2 not tracked directly, use Krumholz+Gnedin 2011 Z-dependent H2 calculation

/* Computes cold gas fractions for particle i using auto-shielding approximation */
/* flag=0 returns total cold gas frac (HI+H2), flag=1 returns HI frac, flag=2 returns H2 frac */
double coldgasfrac(int i, double T, double ndens, int flag)
{
	double Nh,hsmooth;
        double fneut,fH2;
	int InitKernIntTable();
#ifdef LEROY_H2
	double Rmol;
#endif
#ifdef KG11_H2
        double tau_c,chi,sterm;
#endif

#ifdef GRACKLE
	Nh = SphP[i].grHI;	// neutral H fraction (assuming optically thin)
#else
	Nh = 1.+2*yhelium-SphP[i].Ne; if( Nh<0 ) Nh=0; // very approximate!  Better to use with grackle.
#endif
	if( ndens < 1.e-4 ) return Nh;	// Very low density gas can't be self-shielded
	if( T>1.e5 && SphP[i].Sfr==0 ) return Nh;	// only cold, dense gas considered
	if( SphP[i].Sfr>0 && flag==0 ) return FSHIELD;	// SF'ing gas always made up of cold gas
	if( gJH0 < 1.e-30 ) return Nh;	// Ig no ionizing background, can't self-shield against it

/* Rahmati+13 prescription for self-shielded gas */
  double nHss,fG,sigHI;
  double beta,Cp;
  double a=7.982e-11, b=0.7480, T0=3.148, T1=7.036e5;
	double z = 1./All.cf_atime-1.;
	nHcgs = XH * ndens;	/* hydrogen number dens in physical cgs units */
  sigHI = 3.27e-18 * pow(1+z,-0.2);
  nHss = 6.73e-3 * pow(sigHI/2.49e-18,-2./3) * pow(T*1.e-4,0.17) * pow(gJH0*1.e12,2./3) * pow((All.OmegaBaryon/All.Omega0)/0.17,-1./3);
  fG = 0.98 * pow(1+pow(nHcgs/nHss,1.64),-2.28) + 0.02 * pow(1+nHcgs/nHss,-0.84);
/* Calculate cold gas fraction using new GammaHI attenuated by fG; Popping+09 */
  beta = a / (sqrt(T/T0)* pow(1+sqrt(T/T0),1-b) * pow(1+sqrt(T/T1),1+b));	
  Cp = nHcgs*beta/(fG*gJH0);
  Nh = (2*Cp+1-sqrt((2*Cp+1)*(2*Cp+1)-4*Cp*Cp))/(2*Cp);
	//if(ThisTask==0 && Nh>0.1) printf("COLDGASFRAC: %d T=%g nH=%g nHss=%g fG=%g fneut=%g\n",i,log10(T),log10(nHcgs),log10(nHss),fG,Nh);

#ifdef GALSF_SFR_KMT
	if( SphP[i].fH2 > FSHIELD ) SphP[i].fH2 = FSHIELD;	// limit H2 frac to FSHIELD
	if( Nh < SphP[i].fH2 ) Nh = SphP[i].fH2;		// limit total cold frac to H2 frac
	if( flag == 1 ) return Nh-SphP[i].fH2;	
	if( flag == 2 ) return SphP[i].fH2; 		// with KMT, H2 frac is already calculated
#endif
	if( flag == 0 ) return Nh;	// all cold gas (neutral+molecular)

/* If not using KMT, then must calculate the molecular fraction */
        fneut = Nh;
        fH2 = 0.;
        if( SphP[i].Sfr>0 ) {
#ifdef LEROY_H2
	  // Empirical relation from Leroy+08
          Rmol = pow(ndens*T/P0BLITZ,ALPHA0BLITZ);
          fneut = FSHIELD*(1e+08-T)/1e+08/(1+Rmol);
          fH2 = FSHIELD-fneut;
#endif
#ifdef KG11_H2
          // Theoretical relation from Krumholz & Gnedin 2011
          tau_c = ndens*hsmooth*4*All.UnitLength_in_cm * 1./(All.Time*All.Time) * 1.e-21* (P[i].Metallicity[0]/All.SolarAbundances[0])/ 2.3e-24;
          chi = 3.1 * (1+3.1*pow(P[i].Metallicity[0]/All.SolarAbundances[0],0.365)) / 4.1;
          sterm = log(1+0.6*chi+0.01*chi*chi)/(0.6*tau_c);
          fH2 = 1 - 0.75*sterm / (1+0.25*sterm);
          if( fH2 > FSHIELD ) fH2 = FSHIELD; if( fH2 < 0 ) fH2 = 0;
          fneut = FSHIELD-fH2;
#endif
        }

/* Only get here if not using KMT */
	if( flag == 1 ) return fneut;	// using KG11 theoretical version here
	else if( flag == 2 ) return fH2;
	else
	  return 0.0;  //avoid compiler warning (SHOULD NOT HAPPEN!)
}


#endif
