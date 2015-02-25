#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "darkenergy.h"
#include "interface.h"


/* note the underscore! */
// int f77main_(double * hI, double * omegamI, double * omegavI, double * omegakI, double * omegaQI, double * wQI, double * wQpI, double * z1, double * z2, double * D1, double *D2, double *zmaxact, double *zminact, int *iwmodeI);
// The above declaration is included in interface.h.

// int main(int argc, char * argv[]){
// Arguments in order: redshift today, redshift at beginning of N-body sim, omega matter (total, i.e. CDM+baryons), omega dark energy, omega cosmological constant (not used), w_0, w_a, h (Hubble parameter).
double fvel_Interface(double z1, double z2, double omegamI, double omegaQI, double omegavI, double wQI, double wQpI, double hI, double d2, double d2minus,double delz,double zmaxact,double zminact,int iwmodeI,int *err)
{
  
  // double omegamI,omegavI,omegaQI,wQI,wQpI,hI, z1, z2;
  double omegakI;
  
  double OM, OL, OK, DE, Hz, H0, fvel, fspringel; // fspringel is fvel in Gadget units.
  int ret;

  
  /*
  omegamI=0.27;
  omegaQI=0.73;
  omegavI=0.0;
  wQI=-1.0;
  wQpI=0.0;
  hI=0.7;
  */
  
  /*
  // Redshifts between which growth factor calculated:
  z1=0.0;
  z2=60.0;

  // Cosmological parameters:
  omegamI=0.26;
  omegaQI=0.74;
  omegavI=0.0;
  wQI=-0.8;
  wQpI=0.0;
  hI=0.72;
  */
  
  // change seldom (determines span over which finite derivative is calculated for growth factor (should be very small but within numerical accuracy of code):
  // delz=0.000001; // 0.000001 seems a good number for accuracy before numerical artifacts start creeping in (some part of fortran code work probably just with floats).
  
  // do not change:
  //zmaxact=110.0; // make larger than any initial redshift of simulation needed (determines precalculated values, otherwise tables too short).
  //zminact=0.0;
  //iwmodeI=3;
  //na=5000; // direct in FORTRAN by means of parameters.inc file.
  //cspeed=2.99792458e5; // direct in FORTRAN
    
  // automatic initialization:
  DE=0.0;
  omegakI=1.0-omegamI-omegavI-omegaQI;
  OM=omegamI;
  OL=omegaQI;
  OK=omegakI;
  
  w0=wQI;
  wa=wQpI;
  ret = initialize_darkenergy();
  // w0 and wa above are global variables, need to be global for dark energy calculation functions.
  
  // Beware: Lam Hui's linear growth code does not seem to work well
  // (i.e. with dverk error) if wQpI differs from 0.0 by more than 0.15 or so.
  // He needs to revisit this later.

    
  // Doing a first test:  
  //return_code = f77main_(&hI, &omegamI, &omegavI, &omegakI, &omegaQI, &wQI, &wQpI, &z1, &z2, &D1, &D2, &zmaxact, &zminact, &iwmodeI);
  //printf("Redshifts and Growth Factors: %e %e %e %e \n", z1, z2, D1, D2);
  
  
  // Now doing finite difference between growth factors, to get time derivative:
  
  // printf("Now doing IC quantities:\n");
  
  //return_code = f77main_(&hI, &omegamI, &omegavI, &omegakI, &omegaQI, &wQI, &wQpI, &z1, &z2, &D1, &D2, &zmaxact, &zminact, &iwmodeI);
  // printf("Redshifts and Growth Factors: %e %e %e %e \n", z1, z2, D1, D2);

  //d2=D2;

  //return_code = f77main_(&hI, &omegamI, &omegavI, &omegakI, &omegaQI, &wQI, &wQpI, &z1, &z2minus, &D1, &D2, &zmaxact, &zminact, &iwmodeI);
  // printf("Redshifts and Growth Factors: %e %e %e %e \n", z1, z2, D1, D2);

  //d2minus=D2;
  
  
  H0=0.1; // in units km/s/kpc*h. H0=0.1 always in the internal units of Gadget, regardless of value of h (Hubble parameter).
  if (w0>-1.00001 && w0<-0.99999 && wa > -0.00001 && wa < 0.00001)
  {
	DE=1.0;
	fprintf(stderr,"[*]: Setting DE contribution hard to DE=1 (LCDM) to avoid wiggles in spline for insufficiently sampled constant integral by odeint. This is only an accuracy improvement measure.\n");
  }
  else DE=DarkEnergy(1.0/(z2+1.0));
  Hz=H0*sqrt( OM*(z2+1.0)*(z2+1.0)*(z2+1.0) + OK*(z2+1.0)*(z2+1.0) + OL*DE); // need to extend for w(z) by an integrator!
  
  // printf("Hz and Dark Energy Factor: %e %e\n", Hz, DE);
  
  int i;
  double test_z, testDE;
  for (i=0; i<1000; i++)
  {
	test_z=((double) i)/100.0;
	testDE=DarkEnergy(1.0/(test_z+1.0));
	// printf("z=%e, DE=%e\n", test_z, testDE);
  }
  
  
  fvel=((d2minus/d2-1.0)/delz)*Hz; 
    
  // Converting to Gadget-2 units as per Springel's N-GenIC:
  fspringel=fvel*sqrt(1+z2);
  
  free_darkenergy();
  *err = ret;
  
  return fspringel;

}
