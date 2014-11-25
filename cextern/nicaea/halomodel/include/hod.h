/* ---------------------------------------------------------------- *
 * hod.h						            *
 * Martin Kilbinger, Henry J. McCracken, Jean Coupon 2008-2013      *
 * ---------------------------------------------------------------- */

#ifndef __HOD_H
#define __HOD_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "io.h"
#include "errorlist.h"
#include "config.h"
#include "maths.h"

#include "cosmo.h"
#include "nofz.h"
#include "halomodel.h"

/* errors */
#define HOD_BASE        -50000
#define HOD_MSTAR_MAX  HOD_BASE  + 1
#define HOD_NGAL_FIT   HOD_BASE  + 2

#define GM_base      -400000
#define GM_nz_type   GM_base + 1
#define GM_pi_max    GM_base + 2
#define GM_SHMR      GM_base + 3
#define GM_SMF       GM_base + 4
#define GM_GGL       GM_base + 5
#define GM_MSTAR_MAX GM_base + 6
#define GM_ETA       GM_base + 7



/* Limits for xi(r), 1- and 2-halo, in Mpc.h */
#define RMIN1 0.001
#define RMAX1 5.0
#define RMIN2 0.1
#define RMAX2 400.0

#define RMIN  RMIN1
#define RMAX  RMAX2
#define NR    50

#define MAXCHAR 1024
#define NLINES 100

#define NFIELD 100
#define NCHAR 20

/* GG: galaxy-galaxy stuff. GM: galaxy-dar matter stuff (lensing) */
#define GG 0
#define GM 1

#define getDoubleValue(array,col)  atof(array+NCHAR*(col-1))
#define getIntValue(array,col)     atoi(array+NCHAR*(col-1))
#define getCharValue(array,col)    array+NCHAR*(col-1)
#define getLine(array,i)           array+NFIELD*NCHAR*i

#define Nhalodata_t 4
typedef enum {w_of_theta, wp_rp, deltaSigma, smf} halodata_t; 
#define shalodata_t(i) ( \
  i==w_of_theta    ? "w_of_theta" : \
  i==wp_rp         ? "wp_rp"      : \
  i==deltaSigma    ? "deltaSigma" : \
  i==smf           ? "smf"        : \
  "")

#define Nhalomode_t 3
typedef enum {galcorr_var, galcorr_cov, galcorr_log} halomode_t;
#define shalomode_t(i) (   \
 i==galcorr_var ? "galcorr_var" : \
 i==galcorr_cov ? "galcorr_cov" : \
 i==galcorr_log ? "galcorr_log" : \
 "")

#define Nngal_fit_t 5
typedef enum {ngal_log_fit, ngal_lin_fit, ngal_no_fit, ngal_lin_fit_only} ngal_fit_t;
#define sngal_fit_t(i) ( \
  i==ngal_log_fit  ? "ngal_log_fit"  : \
  i==ngal_lin_fit  ? "ngal_lin_fit"  : \
  i==ngal_no_fit   ? "ngal_no_fit"   : \
  i==ngal_lin_fit_only ? "ngal_lin_fit_only" : \
  "")

#define Nintconst_t 2
typedef enum {constant, random_file} intconst_t;
#define sintconst_t(i) ( \
  i==constant  ?    "constant"  : \
  i==random_file  ? "random_file"  : \
  "")


typedef struct { 
  char WTHETA[500];
  char COVNAME[500];
  char DNDZ[500];
  char OUTFILE[500];
  int nbins; 
  double ngal, ngalerr, area_deg, intconst, delta;
  double alpha_min,alpha_max,dalpha;
  double M1_min,M1_max,dM1;
  double Mmin_min,Mmin_max,dMmin;
  double Ngal_min,Ngal_max,dNgal;
  ngal_fit_t ngal_fit_type;
} config_hm_info;


/* ---------------------------------------------------------------- *
 * Global variables and functions                                   *
 * ---------------------------------------------------------------- */


#define ODD 0
#define EVEN 1

#define ABS(a) ((a) < 0 ? -(a) : (a))
#define PARITY(a) (a)%2 ? ODD : EVEN
#define FFTLog_SWAP(a,b) {FFTLog_TMP = (a); (a) = (b); (b) = FFTLog_TMP;}


#define CHANGE(fct) int change_##fct(cosmo_hm*, cosmo_hm*)
CHANGE(w_of_theta);
CHANGE(Pthg);
CHANGE(HOD);
CHANGE(ngd);
CHANGE(vc);
#undef CHANGE


/* ---------------------------------------------------------------- *
 * New types                                                        *
 * ---------------------------------------------------------------- */

/* This structure contains w_theta input data */
typedef struct {
  int nbins, nbins_RR;
  double *th;
  double *w;
  double *werr;   /* Error bars */
  double *wcov;   /* Covariance matrix */
  double *th_RR;  /* theta for random pairs (can be larger than for the data [even recommended]) */
  double *RR;     /* Random pairs (for the integral constraint) */
  double ngal;    /* Galaxy number density*/
  double ngal_err;
} wt_t;



typedef struct Nz_hist
{
  int nbins;
  double *z;
  double *n;
  double zmin;
  double zmax;
  double mean;
  double dz;
  double Norm;
  double max;
} Nz_hist;


/* ---------------------------------------------------------------- *
 * log-likelihood                                                   *
 * ---------------------------------------------------------------- */

double chi2_hm(cosmo_hm *model, halodata_t halodata, halomode_t halomode, const wt_t* data, ngal_fit_t ngal_fit_type,
	       double ngal, double ngalerr, intconst_t intconst_type, error **err);

/* ---------------------------------------------------------------- *
 * w(theta) and wp(rp)                                              *
 * ---------------------------------------------------------------- */

double *woftheta(cosmo_hm *model, pofk_t pofk, double *theta, int Ntheta, int i_bin, int j_bin, error **err);
double *woftheta_FFTLog_total(cosmo_hm *model, double ln_theta_min, double ln_delta_theta, int Ntheta,
			      double **theta, int i_bin, int j_bin, error **err);
double *wp(cosmo_hm *model, pofk_t pofk, const double *rp, int Nrp, double pi_max, int type, error **err);
double int_for_wp(double logr, void *params, error **err);

/* ---------------------------------------------------------------- *
 * Lensing functions (only for leauthaud11 model)                   *
 * ---------------------------------------------------------------- */

double *DeltaSigma(cosmo_hm *model, pofk_t pofk, const double *r, int N, error **err);
double int_for_Sigma(double logr, void *params, error **err);

/* ---------------------------------------------------------------- *
 * Stellar mass function (only for leauthaud11 model)               *
 * ---------------------------------------------------------------- */

double *dNdlogM10stellar(cosmo_hm *model, double *log10Mstellar, int N, error **err);
double *dNdlogM10stellar_c(cosmo_hm *model, double *log10Mstellar, int N, error **err);
double *dNdlogM10stellar_s(cosmo_hm *model, double *log10Mstellar, int N, error **err);

/* ---------------------------------------------------------------- *
 * HOD P(k) and xi(r)                                               *
 * ---------------------------------------------------------------- */

double *xiofr(cosmo_hm *model, pofk_t pofk, const double *r, int N, int type, error **err);
double *xi_1hcs(cosmo_hm *model,double a, const double *r, int N, int type, error **err);
double int_for_xi_1hcs(double logM, void *params, error **err);
double *xi_1hss(cosmo_hm *model,double a, const double *r, int N, int type, error **err);
double P1hss(double k, void *params);
double int_for_P1hss(double logM, void *params, error **err);
double concentration_sat(cosmo_hm *model, double Mh, double a, error **err);
double *xi_2h(cosmo_hm *model,double a, const double *r, int N, int type, error **err);
double *xi_P_NL(cosmo_hm *model, double a, const double *r, int N, error **err);
double FFTLog_P_NL(double k, void *params);
double P2h(double k, void *params);
double int_for_P2h(double logM, void *params, error **err);
double int_for_P2h_dm(double logM, void *params, error **err);
double bias_r(double M, void *params);

/* ---------------------------------------------------------------- *
 * HOD model functions                                              *
 * ---------------------------------------------------------------- */

double log10_fSHMR(double log10Mh, cosmo_hm *model);
double log10_fSHMR_inv_minus_x(double log10Mstar, void *p);
double log10_fSHMR_inv(double log10Mstar, cosmo_hm *model);
double Ngal_c(cosmo_hm *model, double Mh, double Mstellar_min, double Mstellar_max);
double int_for_Ngal_c(double Mstellar, void *params, error **err);
double phi_c_Mstellar(cosmo_hm *model, double log10Mstellar, double Mh, error **err);
double sigma_log_M(cosmo_hm *model, double log10Mstellar, error **err);
double eta_cen(cosmo_hm *model, double Mh, error **err);
double Ngal_s(cosmo_hm *model, double M, double Mstellar_min, double Mstellar_max);
double Ngal(cosmo_hm *model, double M, double Mstellar_min, double Mstellar_max);

double av_Mh_given_Mstar(cosmo_hm *model, double Mstellar, double a, error **err);
double int_for_phi_c_Mh(double log10Mh, void *params, error **err);
double int_for_phi_c_Mh_norm(double log10Mh, void *params, error **err); 

/* ---------------------------------------------------------------- *
 * Deduced quantities                                               *
 * ---------------------------------------------------------------- */
double Mstar_tot_c(cosmo_hm *model, double Mh, double Mstellar_min, double Mstellar_max, error **err);
double int_for_Mstar_tot_c(double logMstar, void *params, error **err);
double Mstar_tot_s(cosmo_hm *model, double Mh, double Mstellar_min, double Mstellar_max, error **err);
double int_for_Mstar_tot_s(double logMstar, void *params, error **err);
double av_gal_bias(cosmo_hm *model, double a, error **err);
double int_for_av_gal_bias(double logM, void *intpar, error **err);
double av_halo_mass(cosmo_hm *model, double a, error **err);
double int_for_av_halo_mass(double logM, void *intpar, error **err);
double mass_weighted_av_stellar_mass(cosmo_hm *model, double a, error **err);
double int_for_mass_weighted_av_stellar_mass(double logM, void *params, error **err);
double int_for_mass_denum_weighted_av_stellar_mass(double logM, void *params, error **err);
double av_stellar_mass(cosmo_hm *model, double a, error **err);
double int_for_av_stellar_mass(double logM, void *params, error **err);
double av_frsat(cosmo_hm *model, double a, error **err);
double int_for_av_frsat(double logM, void *intpar, error **err);
double av_halo_bias(cosmo_hm *model, double a, error **err);
double int_for_av_halo_bias(double logM, void *intpar, error **err);
double dn_dlnM_integrand(double logM, void *intpar, error **err);

/* ---------------------------------------------------------------- *
 * Physical quantities                                              *
 * ---------------------------------------------------------------- */

double vc(cosmo_hm *model,double zmin, double zmax, error **err);
double int_for_vc(double z, void *params, error **err);
double logM_lim(cosmo_hm *model, double a, double r, error **err);
double ngal_triax(cosmo_hm *model, double a, double r,  error **err);
double ngal_den(cosmo_hm *model, double a, double logMlim, double Mstellar_min, double Mstellar_max, error **err);
double int_for_ngal_den(double logM, void *params,  error **err);
double ngal_den_c(cosmo_hm *model, double a, double logMlim, double Mstellar_min, double Mstellar_max, error **err);
double int_for_ngal_den_c(double logM, void *params,  error **err);
double ngal_den_s(cosmo_hm *model, double a, double logMlim, double Mstellar_min, double Mstellar_max, error **err);
double int_for_ngal_den_s(double logM, void *params,  error **err);
double ngal_den_vol(cosmo_hm *model, error **err);
double ngal_weighted(cosmo_hm *model, error **err);

/* ---------------------------------------------------------------- *
 * FFTLog                                                           *
 * ---------------------------------------------------------------- */

double xi_from_Pkr(gsl_function *ak, double r_prime, FFTLog_config *fc, error **err);
void FFTLog(FFTLog_config *fc, const gsl_function *ar_in,double *k, double *ak_out, int dir, error **err);
FFTLog_config *FFTLog_init(int N, double min, double max, double q, double mu);
void FFTLog_free(FFTLog_config *fc);
FFTLog_complex FFTLog_U_mu(double mu, FFTLog_complex z);

/* ---------------------------------------------------------------- *
 * Utils                                                            *
 * ---------------------------------------------------------------- */

void updateFrom_hm(cosmo_hm* avant, cosmo_hm* apres, error **err);
int change_vc(cosmo_hm *avant, cosmo_hm *apres);
int change_ngd(cosmo_hm *avant, cosmo_hm *apres);
int change_HOD(cosmo_hm *avant, cosmo_hm *apres);
int change_Pthg(cosmo_hm* avant, cosmo_hm* apres);
int change_w_of_theta(cosmo_hm *avant, cosmo_hm *apres);
wt_t *init_wt_t(double nlines, int wcov, error **err);
void free_wt_t(wt_t **model);
wt_t *read_wtheta(const char *data_file_name, intconst_t intconst_type, double delta, double intconst, const char *ran_file_name, error **err);
int getStrings(char *line, char *strings, char *delimit, size_t *N, error **err);
void read_config_hm_file(config_hm_info *config, char cname[], error **err);
Nz_hist *get_nz(redshift_t *nz_mk, int n_bin, error **err);
void  free_nz( Nz_hist *nz);
double normalize(double *data, size_t Nbins, double binsize,error **err);

/* ---------------------------------------------------------------- *
 * Obsolete stuff                                                   *
 * ---------------------------------------------------------------- */

double *woftheta_total_OBSOLETE(cosmo_hm *model, double log_theta_min, double log_delta_theta, int Ntheta, error **err);
nz_t *read_dndz_OBSOLETE(char *infile,error **err);
double Ngal_mean(cosmo_hm *model, double a, error **err);

double xir_mwolk(cosmo_hm *model, double a, double r, error **err);
wt_t *read_wtheta_mwolk(const char *name, double delta, double intconst, error **err);
double chi2_mwolk(cosmo_hm *model, const wt_t* data, ngal_fit_t ngal_fit_type,
		  double ngal, double ngalerr, double area, error **err);
double wp_mwolk(cosmo_hm *model, double rp, error **err);
double compute_chisq_wp(cosmo_hm *model, const wt_t *wth, double ngd_obs, double ngd_err,
			ngal_fit_t ngal_fit_type, double *ngd, int dologw, error **err);


#endif
