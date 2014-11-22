/* ============================================================ *
 * decomp_eb.h							*
 * Martin Kilbinger, Liping Fu 2008, 2009			*
 * ============================================================ */

#ifndef __DECOMP_EB_H
#define __DECOMP_EB_H


#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_combination.h>

#include "maths.h"
#include "io.h"
#include "errorlist.h"


/* Error codes */
#define mr_base		-1500
#define mr_range	-1 + mr_base
#define mr_poly		-2 + mr_base
#define mr_func		-3 + mr_base
#define mr_type		-4 + mr_base
#define mr_incompatible -5 + mr_base
#define mr_null         -6 + mr_base
#define mr_dimension    -7 + mr_base
#define mr_file         -8 + mr_base


/* Maximum COSEBI mode. The code is accurate up *
 * to Nmax_cosebi = 13			        */
#define NMAX_COSEBI     20


typedef enum {cheby, cheby2, legen, cosebi} poly_t;
#define spoly_t(i) ( \
  i==cheby  ? "cheby" : \
  i==cheby2 ? "cheby2" : \
  i==legen  ? "legen" : \
  i==cosebi  ? "cosebi" : \
  "")
#define Npoly_t 4

/*
typedef enum {cosebi} filter_t;
#define sfilter_t(i) ( \
  i==cosebu  ? "cosebi" : \
  "")
#define Nfilter_t 1
*/


typedef enum {comb_none, all_sc} mring_combinef_t;
#define smring_combinef_t(i) ( \
   i==comb_none    ? "comb_none" : \
   i==all_sc       ? "all_sc" : \
   "")
#define Nmring_combinef_t 2

typedef enum {inner, outer, fromZ} cov_mode_t;

/* === CFHTLS Wide 3rd data release, Fu&Kilbinger (2010) === */
#define N_FK10 6
#define eta_FK10_SN 0.02
#define eta_FK10_FoM_eta10 0.1
#define eta_FK10_FoM_eta50 0.02


double Cheby(double x, int n, error **err);
double Cheby2(double x, int n, error **err);
double Legen(double x, int n);
double C(double x, int n, poly_t poly, error **err);
double Tp(double x, const double *a, int N, poly_t poly, error **err);

double Fn0(double x, int n, poly_t poly, error **err);
void Fnnu(double x, int n,  poly_t poly, double Fn[], error **err);
double alpha_explicit(double x, int n, double R, poly_t poly, error **err);
double Tm(double x, const double *a, int N, poly_t poly, double R, error **err);

double RR_data(const double *xip, const double *xim, const double *th, const int Nxi,
	       double THETA_MIN, double THETA_MAX, const double *a, int N,
	       poly_t poly, int pm, error **err);

double cov_RR(const double *THETA_MIN, const double *THETA_MAX, const double *a, int N, poly_t poly,
	      const double *theta, const double *cov_xi, int Ntheta, 
	      cov_mode_t cov_mode, int n, int m, double fac, error **err);
double cov_RR_diag_xi(const double *THETA_MIN, const double *THETA_MAX, const double *a, int N, poly_t poly,
		      const double *theta, const double *var_xi, int Ntheta, int islog,
		      error **err);

double chi2_RB_null(const double *RB, const double *covRB, int NRB);

double *read_zeros_norm_cosebi_auto_check(double Psimin, double Psimax, const char *path, error **err);
double *read_zeros_norm_cosebi(const char *rname, double *psimin, double *psimax, error **err);
double Tplog(double x, const double *r, int n, error **err);
double Tplog_c(double z, const double *c, int n, error **err);
double Tmlog(double x, const double *c, int n, error **err);
double an2(int n,  const double *c);
double an4(int n, const double *c);
double dnm(int n, const int m, const double *c);
double sum_combinations(int j, int n, const double *r, error **err);

#endif
