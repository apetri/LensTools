/*
 *  hubble.c
 *  
 *  written by christian wagner
 *
 */

#include "coyote.h"


#define VC 299792.458
#define TWOPI 6.283185307179586
#define PI    3.141592653589793
#define Tcmb  2.725
#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))
#define pow4(x) ((x)*(x)*(x)*(x))

struct cosmo{
  double hub;
  double w0;
  double wa;
  double wm;
  double wb;
} mycosmo;


double hubble(double omega_m, double w0_de, double h_100, double a)
{
   double H0       = h_100 * 100.;
   double Omega_m  = omega_m / h_100 / h_100;
   double Omega_r  = 4.15e-5/pow2(h_100);
   double Omega_X  = 1.0 - Omega_m;
   double Omega_k  = 0.0;

   /* Assuming wa = 0 */
   return  H0*sqrt(Omega_r/pow4(a)+Omega_m/pow3(a)+Omega_k/pow2(a)+Omega_X*pow(a,-3*(1+w0_de)));
}


double co_distance_int(double a, void * params)
{
   struct cosmo mycosmo = *(struct cosmo *) params;

   return  1./(a*a*hubble(mycosmo.wm, mycosmo.w0, mycosmo.hub, a));
}

double co_distance(double omega_m, double w0_de, double h_100, double a)
{
   gsl_function  F;
   double result;
   double epsabs=0.001;
   double epsrel=0.001;
   double abserr;
   size_t neval;
   struct cosmo mycosmo;

   F.function  = &co_distance_int;
   mycosmo.wm  = omega_m;
   mycosmo.w0  = w0_de;
   mycosmo.hub = h_100;
   F.params = &mycosmo;

   gsl_integration_qng (&F, a, 1., epsabs, epsrel, &result, &abserr, &neval);

   return result*VC*h_100;
}

double z_lastscattering(double omega_m, double omega_b)
{
   /*
     Returns z_LS from Hu & White, DampingTail paper.
   */

   double wm = omega_m;
   double wb = omega_b;
   double b1 = 0.0783*pow(wb,(-0.238))/(1+39.5*pow(wb,0.763));
   double b2 = 0.560/(1+21.1*pow(wb,1.81));
   double zls= 1048.*(1+0.00124*pow(wb,(-0.738)))*(1+b1*pow(wm,b2));
   return zls;
}

double soundhorizon(double omega_m, double omega_b)
{
  /*
  A fit to the sound horizon, in Mpc, from Hu & Sugiyama (1995), Eq. B6
  */

   double wm = omega_m;
   double wb = omega_b;
   double wg = 2.4888e-5*pow4(Tcmb/2.73);
   double r0 = 0.75*wb/wg;
   double zeq= 5464.0*(wm/0.135)/pow4(Tcmb/2.725)/(1+0.6851)-1.0;
   double req= r0/(1.+zeq);
   double zLS= z_lastscattering(omega_m, omega_b);
   double rls= r0/(1+zLS);
   double tmp= (sqrt(1.+rls)+sqrt(rls+req))/(1.+sqrt(req));
   tmp= 3997.0*sqrt(wg/wm/wb)*log(tmp);

  //fprintf(stderr, "Sound horizon = %g Mpc (zLS = %g)\n", tmp, zLS);

  return(tmp);
}

double distls(double omega_m, double omega_b, double w0_de, double h_100)
{
   /*
     Returns the distance to LS, in Mpc.
   */

   double zLS = z_lastscattering(omega_m, omega_b);
   /* Used to be call to ang_distance (= co_distance for flat Universe) */
   double dLS = co_distance(omega_m, w0_de, h_100, 1.0/(1.0+zLS)) / h_100;

   return(dLS);
}

double set_omega(double omega, double h_100, int physical)
{
   double w;

   w = omega;
   if (physical == 0) w = w * h_100 * h_100;

   return w;
}

#define PREC 1.0e-3
double getH0fromCMB(double omega_m, double omega_b, double w0_de, int physical)
{
  /*
  solvehh(dLS,cosmo):
  Solves for h given the other cosmological parameters and the distance
  to last scattering (in Mpc).
  */

   double h_100_min, h_100_max, h_100_mid, ymin, ymax, ymid;
   double rs, dLS, wm, wb;

   dLS = 302.4 / pi;   /* Constraint from WMAP7 */

   h_100_min = 0.3;
   h_100_max = 1.0;

   wm = set_omega(omega_m, h_100_min, physical);
   wb = set_omega(omega_b, h_100_min, physical);
   rs   = soundhorizon(wm, wb);
   ymin = distls(wm, wb, w0_de, h_100_min) - dLS * rs;

   wm = set_omega(omega_m, h_100_max, physical);
   wb = set_omega(omega_b, h_100_max, physical);
   rs   = soundhorizon(wm, wb);
   ymax = distls(wm, wb, w0_de, h_100_max) - dLS * rs;

   while (fabs(h_100_max - h_100_min) > PREC) {

      h_100_mid = (h_100_max + h_100_min) / 2.0;

      wm = set_omega(omega_m, h_100_mid, physical);
      wb = set_omega(omega_b, h_100_mid, physical);
      rs   = soundhorizon(wm, wb);
      ymid = distls(wm, wb, w0_de, h_100_mid) - dLS * rs;

      if (ymin*ymid<0) {
	 h_100_max = h_100_mid;
	 ymax      = ymid;
      } else {
	 h_100_min = h_100_mid;
	 ymin      = ymid;
      }

   }

   return h_100_mid;
}
#undef PREC

