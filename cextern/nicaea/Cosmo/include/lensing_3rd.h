/* ============================================================ *
 * lensing_3rd.h						*
 * Martin Kilbinger, Liping Fu 2010.				*
 * ============================================================ */

#ifndef __LENSING_3RD_H
#define __LENSING_3RD_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <fftw3.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>

#include "cosmo.h"
#include "nofz.h"
#include "lensing.h"
#include "errorlist.h"
#include "maths_base.h"


#define tenoverseven   1.428571428571
#define fouroverseven  0.571428571429


/* Minimum scale factor where a non-linear scale can be defined. *
 * See {a,b,c}scocou.						 */
#define A_NL_MIN 0.02

#define s2_min     0.1
#define s2_max     1.0e6
#define N_s2       50
#define epsilon0   1.0e-2
#define N_EFF_MIN -2.0


/* more ce_xyz in smith2.h */
#define ce_wrongmode    -100
#define ce_eof          -101
#define ce_incompatible -102
#define ce_rootbracket  -103

typedef enum {fgauss=0, fpoly=1, ftophat=2, fdelta=3, fxip=4, fxim=5, fall=6} filter_t;

typedef enum {PT=0, SCOCOU=1} bispmode_t;
#define sbispmode_t(i) (   \
  i==PT ? "PT" :           \
  i==SCOCOU ? "scocou01" : \
  "")
#define Nbispmode_t 2

typedef enum {kkk=0, kkg=1, kgg=2, ggg=3} bispfield_t;


typedef struct {

  /* Lensing, including basic cosmology and redshift distribution)(s) */
  cosmo_lens *lens;

  /* ============================================================ *
   * Precomputed stuff (at the moment only one redshift-bin).     *
   * ============================================================ */

  interTable *k_NL;
  double scale_NL_amin;
  interTable2D **B_kappa[3];

  bispmode_t bispmode;
} cosmo_3rd;


typedef struct {
  double r1, r2;
  cosmo_3rd *self;
} cosmo3ANDdoubleANDdouble;

typedef struct {
  double r1, r2;
  cosmo_3rd *self;
  int i, n_bin[3];
} cosmo3ANDtomo;

typedef struct {
  cosmo_3rd *self;
  double R[3];
  filter_t wfilter;
  int n_bin[3];
  int m;
  error **err;
} cosmo3ANDmorestuff;


cosmo_3rd *init_parameters_3rd(double OMEGAM, double OMEGAV, double W0_DE, double W1_DE,
			       double *W_POLY_DE, int N_POLY_DE,
			       double H100, double OMEGAB, double OMEGANUMASS,
			       double NEFFNUMASS, double NORM, double NSPEC,
			       int Nzbin, const int *Nnz, const nofz_t *nofz, double *par_nz,
			       nonlinear_t NONLINEAR, transfer_t TRANSFER,
			       growth_t GROWTH, de_param_t DEPARAM,
			       norm_t normmode, bispmode_t BISPMODE, error **err);
cosmo_3rd *copy_parameters_3rd_only(cosmo_3rd *source, error **err);
void updateFrom_3rd(cosmo_3rd *avant, cosmo_3rd *apres, error **err);
cosmo_3rd *set_cosmological_parameters_to_default_lens_3rd(error **err);
void read_cosmological_parameters_lens_3rd(cosmo_3rd **self, FILE *F, error **err);
void free_parameters_3rd(cosmo_3rd **self);
void dump_param_3rd(cosmo_3rd* self, FILE *F, error **err);

double dcub(double a);
double n_eff(cosmo *, double, int, error **);
double Q3(double);
double temp_NL(double, double, cosmo *, error **);
double scale_NL(cosmo_3rd *, double, error **);

double ascocou(cosmo_3rd *, double, double, double *, error **);
double bscocou(cosmo_3rd *, double, double, double *, error **);
double cscocou(cosmo_3rd *, double, double, double *, error **);

double int_for_B_kappa_bar0(double a, void *intpar, error **err);
double int_for_B_kappa_bar1(double a, void *intpar, error **err);
double int_for_B_kappa_bar2(double a, void *intpar, error **err);
double int_for_B_kappa_bar3(double a, void *intpar, error **err);

double F2eff(cosmo_3rd *, double a, double k1, double k2, double cosphi, error **);
double F2bar(int, double, double, error **);
double F2(int, double, double, error **);
double F2cos(int, double, error **err);

double Q_123(cosmo_3rd *, double, double, double, double, error **);
double bb(cosmo_3rd *, double, double, double, int i_bin, int j_bin, int k_bin, error **);
double B_kappa_bar(cosmo_3rd *self, double s1, double s2, int abc, int i_bin, int j_bin, int k_bin, error **err);

double hept_rtbis(double (*func)(double,double,cosmo*,error**), double, double,
		  double, double, cosmo *, error **);
double B_delta(cosmo_3rd *self, double k1, double k2, double cosbeta, double a,
	       error **err);

double Uhat_one(double x, filter_t wfilter);
double Uhat(double eta, filter_t wfilter);
void permute3(double *x, int offset);
double int_for_map3_3d(double x[], size_t dim, void *intpar);
double map3_perm(cosmo_3rd *self, double R[3], int i_bin, int j_bin, int k_bin, filter_t wfilter, error **err);
double map3(cosmo_3rd *self, double R[3], int i_bin, int j_bin, int k_bin, filter_t wfilter, error **err);


double lensing_signal_3rd(cosmo_3rd *self, double theta[3], int i_bin, int j_bin, int k_bin, error **err);
double chi2_lensing_3rd(cosmo_3rd *self, datcov *dc, error **err);



#endif

