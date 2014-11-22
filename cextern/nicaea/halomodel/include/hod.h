/* ============================================================ *
 * hod.h							*
 * Martin Kilbinger, Henry J. McCracken, Jean Coupon 2008-2012  *
 * ============================================================ */

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
#include <fftw3.h>

#include "io.h"
#include "errorlist.h"
#include "config.h"
#include "maths.h"

#include "cosmo.h"
#include "nofz.h"
#include "halomodel.h"

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

#define getDoubleValue(array,col)  atof(array+NCHAR*(col-1))
#define getIntValue(array,col)  atoi(array+NCHAR*(col-1))
#define getCharValue(array,col) array+NCHAR*(col-1)
#define getLine(array,i) array+NFIELD*NCHAR*i

#define Nhalodata_t 3
typedef enum {woftheta, wp_rp, deltaSigma} halodata_t; 
#define shalodata_t(i) (\
  i==woftheta      ? "woftheta" :\
  i==wp_rp         ? "wp_rp"    :\
  i==deltaSigma    ? "deltaSigma" :\
"")

#define Nhalomode_t 3
typedef enum {galcorr_var, galcorr_cov, galcorr_log} halomode_t;
#define shalomode_t(i) (   \
 i==galcorr_var ? "galcorr_var" : \
 i==galcorr_cov ? "galcorr_cov" : \
 i==galcorr_log ? "galcorr_log" : \
 "")

#define Nngal_fit_t 5
typedef enum {ngal_log_fit, ngal_lin_fit, ngal_no_fit, ngal_lin_fit_only, ngal_match_M1} ngal_fit_t;
#define sngal_fit_t(i) ( \
  i==ngal_log_fit  ? "ngal_log_fit"  : \
  i==ngal_lin_fit  ? "ngal_lin_fit"  : \
  i==ngal_no_fit   ? "ngal_no_fit"   : \
  i==ngal_lin_fit_only ? "ngal_lin_fit_only" : \
  i==ngal_match_M1 ? "ngal_match_M1" : \
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


/*OBSOLETE
typedef struct {
  double Mmin;
  double M1, M0, sigma_log_M;
  double alpha;
  double Ngal;
} haloparams;
*/

/*----------------------------------------------------------------*
 *Global variables and functions                                  *
 *----------------------------------------------------------------*/

double FFTLog_TMP;


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


/*----------------------------------------------------------------*
 *New types                                                       *
 *----------------------------------------------------------------*/

/* This structure contains w_theta input data */
typedef struct {
  int nbins;
  double *th;
  double *w;
  double *werr;   /* Error bars */
  double *wcov;   /* Covariance matrix */
  double *RR;     /* Random pairs */
  double ngal;    /* Galaxy number density*/
  double ngal_err;
} wt_t;


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


typedef struct gsl_int_params
{
  void *params;
  funcwithpars func;
  error **err;

} gsl_int_params;


/*----------------------------------------------------------------*
 *Melody's functions                                              *
 *----------------------------------------------------------------*/

double xir_mwolk(cosmo_hm *self, double a, double r, error **err);
wt_t *read_wtheta_mwolk(const char *name, double delta, double intconst, error **err);
double chi2_mwolk(cosmo_hm *model, const wt_t* data, ngal_fit_t ngal_fit_type,
		  double ngal, double ngalerr, double area, error **err);
double wp(cosmo_hm *self, double rp, error **err);
double compute_chisq_wp(cosmo_hm *model, const wt_t *wth, double ngd_obs, double ngd_err,
			ngal_fit_t ngal_fit_type, double *ngd, int dologw, error **err);
/*----------------------------------------------------------------*
 *log-likelihood                                                  *
 *----------------------------------------------------------------*/

double chi2_DeltaSigma(cosmo_hm *model, halomode_t halomode, const wt_t* data, error **err);
double chi2_hm(cosmo_hm *model, halomode_t halomode, const wt_t* data, ngal_fit_t ngal_fit_type,
	       double ngal, double ngalerr, double area, error **err);

/*----------------------------------------------------------------*
 *w(theta)                                                        *
 *----------------------------------------------------------------*/

double *woftheta_FFTLog(cosmo_hm *model, pofk_t pofk, double *theta, int Ntheta, error **err);
double *woftheta_FFTLog_total(cosmo_hm *model, double log_theta_min, double log_delta_theta, int Ntheta, error **err);

/*----------------------------------------------------------------*
 *HOD P(k) and xi(r)                                              *
 *----------------------------------------------------------------*/

double *xiofr(cosmo_hm *model, pofk_t pofk, double *r, int N, error **err);
double *xi_1hcs(cosmo_hm *model,double a,double *r, int N, error **err);
double int_for_xi_1hcs(double logM, void *params, error **err);
double *xi_1hss(cosmo_hm *model,double a,double *r, int N, error **err);
double P1hss(double k, void *params);
double int_for_P1hss(double logM, void *params, error **err);
double *xi_2h(cosmo_hm *model,double a,double *r, int N, error **err);
double *xi_P_NL(cosmo_hm *model, double a, double *r, int N, error **err);
double FFTLog_P_NL(double k, void *params);
double P2h(double k, void *params);
double int_for_P2h(double logM, void *params, error **err);
double bias_r(double M, void *params);

/*----------------------------------------------------------------*
 *HOD model functions                                             *
 *----------------------------------------------------------------*/

double log_fSHMR_inv_minus_Mh(double log10Mstar, void *p);
double log_fSHMR(double Mh, cosmo_hm *model);
double Ngal_c(cosmo_hm *model, double M);
double Ngal_s(cosmo_hm *model, double M);
double Ngal(cosmo_hm *model, double M);

/*----------------------------------------------------------------*
 *Deduced quantities                                              *
 *----------------------------------------------------------------*/

double av_gal_bias(cosmo_hm *model, double a, error **err);
double int_for_av_gal_bias(double logM, void *intpar, error **err);
double av_halo_mass(cosmo_hm *model, double a, error **err);
double int_for_av_halo_mass(double logM, void *intpar, error **err);
double av_frsat(cosmo_hm *model, double a, error **err);
double int_for_av_frsat(double logM, void *intpar, error **err);
double av_halo_bias(cosmo_hm *model, double a, error **err);
double int_for_av_halo_bias(double logM, void *intpar, error **err);
double Ngal_mean(cosmo_hm *model, double a, error **err);
double dn_dlnM_integrand(double logM, void *intpar, error **err);

/*----------------------------------------------------------------*
 *Physical quantities                                             *
 *----------------------------------------------------------------*/

double vc(cosmo_hm *model,double zmin, double zmax, error **err);
double int_for_vc(double z, void *params, error **err);
double M_vir(cosmo_hm *model, double r, double a);
double logM_lim(cosmo_hm *model, double a, double r, error **err);
double ngal_triax(cosmo_hm *model, double a, double r,  error **err);
double ngal_den(cosmo_hm *model, double a, double logMlim, error **err);
double int_for_ngal_den(double logM, void *intpar,  error **err);
double ngal_den_vol(cosmo_hm *model, error **err);
double ngal_weighted(cosmo_hm *model, error **err);

/*----------------------------------------------------------------*
 *FFTLog                                                          *
 *----------------------------------------------------------------*/

double xi_from_Pkr(gsl_function *ak, double r_prime, FFTLog_config *fc, error **err);
void FFTLog(FFTLog_config *fc, const gsl_function *ar_in,double *k, double *ak_out, int dir, error **err);
FFTLog_config *FFTLog_init(int N, double min, double max, double q, double mu);
void FFTLog_free(FFTLog_config *fc);
FFTLog_complex FFTLog_U_mu(double mu, FFTLog_complex z);

/*----------------------------------------------------------------*
 * Utils                                                          *
 *----------------------------------------------------------------*/

void updateFrom_hm(cosmo_hm* avant, cosmo_hm* apres, error **err);
int change_vc(cosmo_hm *avant, cosmo_hm *apres);
int change_ngd(cosmo_hm *avant, cosmo_hm *apres);
int change_HOD(cosmo_hm *avant, cosmo_hm *apres);
int change_Pthg(cosmo_hm* avant, cosmo_hm* apres);
int change_w_of_theta(cosmo_hm *avant, cosmo_hm *apres);
wt_t *init_wt_t(double nlines, int wcov, error **err);
void free_wt_t(wt_t **self);
wt_t *read_wtheta(const char *name, double delta, double intconst, error **err);
//wt_t *read_wtheta(config_hm_info *config, error **err);

int getStrings(char *line, char *strings, char *delimit, size_t *N, error **err);

nz_t *read_dndz(char *infile,error **err);
void read_config_hm_file(config_hm_info *config, char cname[], error **err);
Nz_hist *get_nz(redshift_t *nz_mk,error **err);
void  free_nz( Nz_hist *nz);
double normalize(double *data, size_t Nbins, double binsize,error **err);
double int_gsl(funcwithpars func,void *params, double a, double b, double eps, error **err);
double integrand_gsl(double x,void *p);

#endif
