/* ============================================================ *
 * sn1a.h							*
 * ============================================================ */


#ifndef __SN1A_H
#define __SN1A_H


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "cosmo.h"
#include "errorlist.h"
#include "io.h"
#include "maths.h"
#include "mvdens.h"


/* Number of derivatives */
#define NDER 8

/* Number of light curve parameters */
//#define NLCP 3
#define NLCP 4

/* Chi^2 methods */
typedef enum {chi2_simple, chi2_Theta2_denom_fixed, chi2_no_sc, chi2_betaz, chi2_dust, chi2_residual} chi2mode_t;
#define schi2mode_t(i) ( \
  i==chi2_simple ? "chi2_simple" : \
  i==chi2_Theta2_denom_fixed ? "chi2_Theta2_denom_fixed" : \
  i==chi2_no_sc ? "chi2_no_sc" : \
  i==chi2_betaz ? "chi2_betaz" : \
  i==chi2_dust ? "chi2_dust" : \
  i==chi2_residual ? "chi2_residual" : \
  "" )
#define Nchi2mode_t 6

/* Data formats */
typedef enum {SNLS_firstyear, SN_SALT} sndatformat_t;
#define ssndatformat_t(i) ( \
  i==SNLS_firstyear ? "SNLS_firstyear" : \
  i==SN_SALT ? "SN_SALT" : \
  "")
#define Nsndatformat_t 2

/* Cosmological and nuisance parameters */
typedef struct {

  cosmo *cosmo;                 /* Cosmological pararmeters (from cosmo.h) */
  double Theta1[NDER];          /* Systematic uncertainties, regarding photometric calibration and filter */
  double Theta2[NLCP];          /* 0=-M, 1=alpha, 2=-beta   */
  double Theta2_denom[NLCP];    /* fixed alppha, beta in sigma^2(musb), Theta2_denom[0] unused */
  double beta_d;                /* Dust absorption coefficient, = R_B     */

  double stretch, color;        /* True stretch, color */

  chi2mode_t chi2mode;

} cosmo_SN;

/* SNIa data type (an individual SN) */
typedef struct {

  double z;                     /* redshift				  */
  double musb;			/* brightness in B restframe		  */
  double s,c;			/* stretch color			  */

  double dust;                  /* Intergalactic dust absorption          */

  double derivative[NDER];      /* derivatives w.r.t systematics	  */
 
  double cov[NLCP][NLCP];       /* covariance matrix of (musb, 1, s-1, c) */
  double dl, mu_c;

} SnData;

/* SNIa sample type (contains all data) */
typedef struct {

  double int_disp;
  double sig_mu_pec_vel, logdetW1;
  double *W1;
  int W1dim, Nsample;

  SnData *data;

} SnSample;

/* Initialisation */
cosmo_SN *init_parameters_SN(double OMEGAM, double OMEGADE, double W0_DE, double W1_DE,
			     double *W_POLY_DE, int N_POLY_DE,
			     double H100, double OMEGAB, double OMEGANUMASS,
			     double NEFFNUMASS, double NORM, double NSPEC,
			     nonlinear_t NONLINEAR, transfer_t TRANSFER, growth_t GROWTH,
			     de_param_t DEPARAM, norm_t normmode,
			     chi2mode_t CHI2MODE, double THETA1[], double THETA2[], double BETA_D,
			     double AMIN, error **err);
cosmo_SN *copy_parameters_SN_only(cosmo_SN *source, error **err);
void read_cosmological_parameters_SN(cosmo_SN **self, FILE *F, error **err);
void updateFrom_SN(cosmo_SN* avant, cosmo_SN* apres, error **err);
cosmo_SN* set_cosmological_parameters_to_default_SN(error **err);
cosmo_SN* set_cosmological_parameters_to_best_fit_SNLS_WMAP5(error **err);
cosmo_SN* set_cosmological_parameters_to_best_fit_SNLS(error **err);
cosmo_SN* set_cosmological_parameters_to_best_fit_Union(error **err);
cosmo_SN* set_cosmological_parameters_to_EdS_SN(error **err);
void free_parameters_SN(cosmo_SN **self);

void dump_param_SN(cosmo_SN *self, FILE *F);

/* IO */
void readSnData(char *line, SnData *sndata, sndatformat_t sndatformat, error **err);
SnSample *SnSample_read(const char *FileName, sndatformat_t sndatformat, error **err);
void out_SnSample(const SnSample *sn, FILE *F);
void out_model(const cosmo_SN *cosmo, FILE *F, error **err);

void SetDl(cosmo_SN *self, SnSample *sn, error **err);
double DistModVariance(const SnData *snd, const double *Theta2);
double vect_scalar_product(const double *x, const double *A, int N);

double int_for_Nhalo_z(double z, void *intpar, error **err);
double Nhalo_z(cosmo *cosmo, double z, error **err);

double chi2_SN_residual(const cosmo_SN *cosmo, const SnSample *sn, error **err);
double chi2_SN(const cosmo_SN *cosmo, const SnSample *sn, mvdens *data_beta_d,
	       int wTheta1, int add_logdetCov, error **err);

double distance_module(cosmo *self, double dlum, error **err);

#endif /* __SN1A_H */
