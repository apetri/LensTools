/* ============================================================ *
 * halomodel.h							*
 * Martin Kilbinger 2006-2009					*
 * ============================================================ */

#ifndef __HALOMODEL_H
#define __HALOMODEL_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <fftw3.h>


#include <gsl/gsl_sf_erf.h>

#include "io.h"
#include "errorlist.h"
#include "config.h"
#include "maths.h"

#include "cosmo.h"
#include "nofz.h"

#define hm_base     -1900
#define hm_hodtype   hm_base + 1
#define hm_Mmin      hm_base + 2
#define hm_pofk      hm_base + 3
#define hm_nfw       hm_base + 4
#define hm_par       hm_base + 5
#define hm_overflow  hm_base + 6
#define hm_io        hm_base + 7
#define hm_zbin      hm_base + 8
#define hm_alpha     hm_base + 9
#define hm_negative  hm_base + 10
#define hm_zmean_2h  hm_base + 11
#define hm_halo_bias hm_base + 12
#define hm_undef     hm_base + 13
#define hm_gsl_int   hm_base + 14

/* Ranges of interpolation tables */
#define k_max_HOD     3336.0

/* Present critical density [M_sol h^2 / Mpc^3] */
#define rho_c0  2.7754e11

/* Mass limits for integration over mass functions */
#define logMmin (3.0*log(10.0))
#define logMmax (16.0*log(10.0))

/* Number of steps for scale-factor-integration (redshift) */
#define Na_hm 20

/* Bit-coded power spectrum types */
typedef enum {pofk_undef=-1, pl=1, pnl=2, p1hdm=4, p2hdm=8, pthdm=16, p1hg=32, p2hg=64, pthg=128,
	      p1hgcs=256, p1hgss=512, pstellar=1024} pofk_t;

/* Halo mass function type */
typedef enum {ps, st, st2, j01} massfct_t;
#define smassfct_t(i) ( \
  i==ps  ? "ps" : \
  i==st  ? "st" : \
  i==st2 ? "st2" : \
  i==j01 ? "j01" : \
  "")
#define Nmassfct_t 4

/* Halo bias type */
typedef enum {halo_bias_sc, halo_bias_tinker05, halo_bias_tinker10} halo_bias_t;
#define shalo_bias_t(i) ( \
   i==halo_bias_sc       ? "halo_bias_sc" : \
   i==halo_bias_tinker05 ? "halo_bias_tinker05" : \
   i==halo_bias_tinker10 ? "halo_bias_tinker10" : \
   "")
#define Nhalo_bias_t 3

/* HOD (Halo occupation distribution) type */
#define Nhod_t 5
typedef enum {hod_none, hamana04, berwein02, berwein02_hexcl, leauthaud11} hod_t;
#define shod_t(i) (		\
 i==hod_none  ? "hod_none" :	\
 i==hamana04  ? "hamana04" :	\
 i==berwein02 ? "berwein02" :	\
 i==berwein02_hexcl ? "berwein02_hexcl" :\
 i==leauthaud11 ? "leauthaud11" :	\
 "")


/* ---------------------------------------------------------------- *
 * Global variables and functions                                   *
 * ---------------------------------------------------------------- */

typedef struct FFTLog_complex
{
  double re;
  double im;
  double amp;
  double arg;
}  FFTLog_complex;

typedef struct {
  int N;
  fftw_plan p_forward;
  fftw_plan p_backward;
  fftw_complex *an;
  fftw_complex *ak;
  fftw_complex *cm;
  fftw_complex *um;
  fftw_complex *cmum;
  double min;
  double max;
  double q;
  double mu;
  double kr;
} FFTLog_config;



typedef struct {
  
  cosmo *cosmo;
  redshift_t *redshift;
  double zmin, zmax;
  
  /* Dark matter halo profile                                                      */
  double c0;			/* concentration parameter                         */
  double alpha_NFW;		/* density slope                                   */
  double beta_NFW;		/* concentration slope as fct of mass              */
  massfct_t massfct;            /* halo mass function				   */
  halo_bias_t halo_bias;        /* Halo bias                                       */
  
  /* Mass function parameters (Sheth&Torman). Do not set manually, they are set    *
   * in set_massfct() according to enum massfct.				   */
  double nmz_a;                 /* Called q in CS02                                */
  double nmz_p;                 /* a=1, p=1/2 is Press-Schechter mass fct.         */
  
  /* HOD (halo occupation distribution) parameters				   */
  hod_t hod;                    /* HOD type                                        */
  double M1, M0, sigma_log_M;
  double M_min;
  double alpha;
  double pi_max;  
  double eta;                   /* central galaxy proportion                       */

  /* galaxy-galaxy lensing and wp(rp) */
  double log10Mhalo;
  double coord_phys;

  /* For Leauthaud11 model                                                         */
  double beta,delta,gamma,Mstar0;
  double beta_sat,B_sat,beta_cut,B_cut;
  double x;                     /* any parameter to propagate if needed            */
  double Mstellar_min, Mstellar_max;
  double fcen1, fcen2;
  
  /* Precomputed stuff */
  double A;			/* Mass function normalisation                     */
  double Mstar;                 /* M_*(a=1.0)					   */
  interTable2D *Pthdm;
  interTable *xir;
  interTable *xi_dm;
  interTable2D *rhohat;
  splineTable* sigRsqr;
  double a_xir;

  /* FFTLOG flag - OBSOLETE */
  int FFTLog;

} cosmo_hm;


typedef struct {
  
  cosmo    *cosmo;
  cosmo_hm *model;
  double   a, r, k, ng, ngp, eps, c;
  double  logMlim, bias_fac, Mh, Mstellar, Mstellar_min, Mstellar_max;
  double  M, r_vir, *kk;
  error    **err;
  
  double   logrmin, logrmax, rp, xi;
  gsl_interp_accel *acc;
  gsl_spline *spline;
  
  int i, j, type, asymptotic, logintegrate;
  
  double   (*bias_func)(double, void *);
} cosmo_hm_params;



typedef struct gsl_int_params
{
  void *params;
  funcwithpars func;
  error **err;

} gsl_int_params;



typedef struct { 
  double *z;
  double *fac;
  double *ypn; /* for spline interpolation */
  double zm;   /* average weighted redshift*/ 
  int nbins;
} nz_t;

cosmo_hm* init_parameters_hm(double OMEGAM, double OMEGADE, double W0_DE, double W1_DE,
			     double *W_POLY_DE, int N_POLY_DE,
			     double H100, double OMEGAB, double OMEGANUMASS, 
			     double NEFFNUMASS, double NORM, double NSPEC,
			     int Nzbin, const int *Nnz, const nofz_t *nofz, double *par_nz,
			     double zmin, double zmax,
			     nonlinear_t NONLINEAR, transfer_t TRANSFER, growth_t GROWTH,
			     de_param_t DEPARAM, norm_t normmode,
			     double C0, double ALPHANFW, double BETANFW, massfct_t MASSFCT, halo_bias_t HALO_BIAS,
			     double M_min, double M1, double M0, double sigma_log_M, double alpha,
			     double Mstar0, double beta, double delta, double gamma, double B_cut, double B_sat, 
			     double beta_cut, double beta_sat, double Mstellar_min, double Mstellar_max, double eta,
			     double fcen1, double fcen2,
			     hod_t HOD, double pi_max, error **err);

cosmo_hm* copy_parameters_hm_only(cosmo_hm* source, error **err);
cosmo_hm *copy_parameters_hm(cosmo_hm *source, error **err);
void read_cosmological_parameters_hm(cosmo_hm **model, FILE *F, error **err);
cosmo_hm *set_cosmological_parameters_to_default_hm(error **err);
void free_parameters_hm(cosmo_hm** model);

void set_massfct(massfct_t massfct, double *nmz_a, double *nmz_p, error **err);
void dump_param_only_hm(cosmo_hm* model, FILE *F);
void dump_param_hm(cosmo_hm* model, FILE *F, error **err);

double sm2_rtbis(double (*func)(double, void *, error **), double x1, double x2,
		 double xacc, void *param, error **err);

/* From nrcomplex.h,c */
#ifndef _DCOMPLEX_DECLARE_T_
typedef struct DCOMPLEX {double r,i;} dcomplex;
#define _DCOMPLEX_DECLARE_T_
#endif /* _DCOMPLEX_DECLARE_T_ */
dcomplex Complex(double re, double im);
dcomplex Cadd(dcomplex a, dcomplex b);
dcomplex Cmul(dcomplex a, dcomplex b);
dcomplex Cdiv(dcomplex a, dcomplex b);
dcomplex RCmul(double x, dcomplex a);
void sm2_cisi(double x, double *ci, double *si, error **err);


double delta_c(cosmo *model, double a, error **err);
double bis_Mstar(double logM, void *param, error **err);
double bis_Mstar_a(double logM, void *param, error **err);
double Mstar(cosmo_hm *model, error **err);
double Mstar_a(cosmo_hm *model, double a, error **err);
double concentration(cosmo_hm *model, double Mh, double a, error **err);
double Delta_vir(cosmo_hm *model, double a);
double dsigma_R_sqr_dR(cosmo_hm *model, double R, error **err);
double nufnu(cosmo_hm *model, double nu, int asymptotic, error **err);
double nufnu_j01(double x);
double sigma_R_sqr(cosmo_hm *model, double R, error **err);
double sigmasqr_M(cosmo_hm *model, double M, error **err);

double dsigma_m1_dlnM(cosmo_hm *model, double M, error **err);
double dnu_dlnM(cosmo_hm *model, double M, double a, error **err);
double dn_dlnM_lnM(double logM, void *intpar, error **err);
double dn_dlnM_uf(double M, cosmo_hm *model, double a, error **err);
double dn_dlnM(double M, void *intpar, error **err);

double r_vir(cosmo_hm *model, double M, double a, error **err);
double M_vir(cosmo_hm *model, double r_vir, double a, error **err);
double Delta_h(cosmo_hm *model, double a, error **err);
double rho_crit(cosmo_hm *model, double a, error **err);
double rho_crit_halo(cosmo_hm *model, double a, error **err);
double Omega_m_halo(cosmo_hm *model, double a, error **err);
double rho_halo(cosmo_hm *model, double r, double a, double Mh, double c, error **err);

double DeltaSigma_WB2000(cosmo_hm *model, double r, const double a, const double M,  double c,  double Delta, error **err);
double g_inf(double x, error **err);
double g_sup(double x, error **err);

double int_for_rhohat(double, void *, error **err);
double rhohat_halo(cosmo_hm *model, double k, double M, double a, double c, error **err);


double halo_bias(cosmo_hm *model, double M, double a, int k, error **err);
double bias(cosmo_hm *model, double M, double a, int k, error **err);
double bias_tinker(cosmo_hm *model, double M, double a, error **err);
double bias_tinker10(cosmo_hm *model, double M, double a, error **err);
double int_for_bias_norm(double logM, void *intpar, error **err);
double bias_norm(cosmo_hm *model, double a, error **err);

double int_for_M_ij(double, void *, error **);
double M_ij(cosmo_hm *model, int i, int j, double a, const double *k, error **err);
double P1h_dm(cosmo_hm *model, double a, double k, error **err);
double P2h_dm(cosmo_hm *model, double a, double k, error **err);

double xi_dm_NL_OBSOLETE(cosmo_hm *model, double a, double r, error **err);  // non-linear DM xi
double int_for_xi_dm_NL_OBSOLETE(double k, void *intpar,  error **err);  // non-linear DM xi

#define CHANGE(fct) int change_##fct(cosmo_hm*, cosmo_hm*)

/* ---------------------------------------------------------------- *
 * Utils                                                            *
 * ---------------------------------------------------------------- */

double int_gsl(funcwithpars func,void *params, double a, double b, double eps, error **err);
double integrand_gsl(double x,void *p);

CHANGE(massfct);
CHANGE(massfct_params);
CHANGE(halo_bias);
CHANGE(sigma_R_sqr);
CHANGE(Mstar);
CHANGE(rhohat_halo);
CHANGE(Pth);

#undef CHANGE


#endif

