/* ============================================================ *
 * lensing.h							*
 *								*
 * Martin Kilbinger, Karim Benabed 2006-2012			*
 * ============================================================ */

#ifndef __LENSING_H
#define __LENSING_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_rng.h>

#include "errorlist.h"
#include "maths.h"
#include "io.h"
#include "mvdens.h"
#include "par.h"

#include "hod.h"
#include "halomodel.h"
#include "cosmo.h"
#include "nofz.h"
#include "decomp_eb.h"
#include "reduced_fit.h"


/* Dimensions of interpolation tables */
/* N_s was increased from 200 to 400, for linear tabulation of P_kappa */
#define N_s     400
#define N_theta 100

/* Ranges of interpolation table for convergence power spectrum. *
 * Power-law extrapolation outside these ranges.		 */
#define s_min     1.0e-2
#define s_max     1.0e6

/* Ranges of interpolation table for reduced-shear correction    *
 * power spectrum. No extrapolation outside these ranges.	 */
#define ELL_MIN_REDUCED            0.1
#define ELL_MAX_REDUCED            2.0e5
#define THETA_P_MIN_REDUCED        (0.1*arcmin)
#define THETA_M_MIN_REDUCED        (0.5*arcmin)
#define THETA_MAP_MIN_REDUCED      (0.2*arcmin)
#define THETA_GSQR_MIN_REDUCED     (0.1*arcmin)
#define THETA_MAPGAUSS_MIN_REDUCED (0.1*arcmin)
#define THETA_MAX_REDUCED          (1000.0*arcmin)

#define NELL_REDUCED 50


#define lensing_base                     -1400
#define lensing_inconsistent             -1 + lensing_base
#define lensing_baryon_fraction          -2 + lensing_base
#define lensing_tomoij			 -3 + lensing_base
#define lensing_initialised              -4 + lensing_base
#define lensing_unknown                  -5 + lensing_base
#define lensing_pm                       -6 + lensing_base
#define lensing_type                     -7 + lensing_base
#define lensing_fastxi                   -8 + lensing_base
#define lensing_nperm                    -9 + lensing_base
#define lensing_range                   -10 + lensing_base
#define lensing_cosebi_n_max            -11 + lensing_base
#define lensing_ia                      -12 + lensing_base
#define lensing_angle_format            -13 + lensing_base
#define lensing_nzbin                   -14 + lensing_base

/* If Ob/Oc > BARYON_FRAC, chi2 produces an error */
#define BARYON_FRAC 0.75


/* Intrinsic alignment, constant amplitude C_1 * rho_crit, *
 * with C_1 = 5e-14 h^2 Mpc^3/M_sol.                       */
#define ia_c1_rho_crit     0.0134


typedef enum {xipm, xip, xim, map2poly, map2gauss, gsqr, decomp_eb, nofz, pkappa, map3gauss,
	      map3gauss_diag, map2gauss_map3gauss_diag, map2gauss_map3gauss,
	      decomp_eb_map3gauss_diag, decomp_eb_map3gauss}
  lensdata_t;
#define slensdata_t(i) ( \
 i==xipm      ? "xipm" : \
 i==xip       ? "xip"  : \
 i==xim       ? "xim"  : \
 i==map2poly  ? "map2poly" : \
 i==map2gauss ? "map2gauss" : \
 i==gsqr      ? "gsqr" : \
 i==decomp_eb ? "decomp_eb" : \
 i==nofz      ? "nofz" : \
 i==pkappa    ? "pkappa" : \
 i==map3gauss ? "map3gauss" : \
 i==map3gauss_diag ? "map3gauss_diag" : \
 i==map2gauss_map3gauss_diag ? "map2gauss_map3gauss_diag" : \
 i==map2gauss_map3gauss ? "map2gauss_map3gauss" : \
 i==decomp_eb_map3gauss_diag ? "decomp_eb_map3gauss_diag" : \
 i==decomp_eb_map3gauss ? "decomp_eb_map3gauss" : \
 "")
#define Nlensdata_t 15

typedef enum {decomp_eb_none, FK10_SN, FK10_FoM_eta10, FK10_FoM_eta50, COSEBIs_log} decomp_eb_filter_t;
#define sdecomp_eb_filter_t(i) (		\
 i==decomp_eb_none ? "none" : \
 i==FK10_SN        ? "FK10_SN" : \
 i==FK10_FoM_eta10 ? "FK10_FoM_eta10" : \
 i==FK10_FoM_eta50 ? "FK10_FoM_eta50" : \
 i==COSEBIs_log    ? "COSEBIs_log" : \
 "")
#define Ndecomp_eb_filter_t 5

/* The following arrays are defined in decomp_eb.c */
extern const double a_FK10_SN[], a_FK10_FoM_eta10[], a_FK10_FoM_eta50[];
// r_COSEB[];

typedef enum {angle_center, angle_mean, angle_wlinear, angle_wquadr} lensformat_t;
#define slensformat_t(i) ( \
 i==angle_center  ? "angle_center" : \
 i==angle_mean    ? "angle_mean" : \
 i==angle_wlinear ? "angle_wlinear" : \
 i==angle_wquadr  ? "angle_wquadr" : \
 "")
#define Nlensformat_t 4

typedef enum {cov_const, cov_ESH09} cov_scaling_t;
#define scov_scaling_t(i) ( \
  i==cov_const ? "cov_const" : \
  i==cov_ESH09 ? "cov_ESH09" : \
  "__undef__")
#define Ncov_scaling_t 2

typedef enum {reduced_none, reduced_K10} reduced_t;
#define sreduced_t(i) ( \
  i==reduced_none ? "none" : \
  i==reduced_K10  ? "K10" : \
  "")
#define Nreduced_t 2

/* Intrinsic alignment model */
typedef enum {ia_none, ia_HS04} ia_t;
#define sia_t(i) ( \
  i==ia_none ? "none" :	\
  i==ia_HS04 ? "HS04" : \
  "")
#define Nia_t 2

/* Bit-coded IA terms */
typedef enum {ia_undef, ia_GI_II, ia_only_GI, ia_only_II} ia_terms_t;
#define sia_terms_t(i) (      \
  i==ia_undef   ? "undef" :   \
  i==ia_GI_II   ? "GI_II" :   \
  i==ia_only_GI ? "only_GI" : \
  i==ia_only_II ? "only_II" : \
  "")
#define Nia_terms_t 4
 
typedef enum {second_order=2, third_order=3} order_t;

typedef enum {tomo_all, tomo_auto_only, tomo_cross_only} tomo_t;
#define stomo_t(i) ( \
 i==tomo_all        ? "tomo_all" : \
 i==tomo_auto_only  ? "tomo_auto_only" : \
 i==tomo_cross_only ? "tomo_cross_only" : \
 "")
#define Ntomo_t 3

typedef struct {
  int n_max;
  double th_min, th_max;
  char path[1024];
} cosebi_info_t;

typedef struct {

  /* Basic cosmology */
  cosmo *cosmo;

  /* Redshift distribution(s) */
  redshift_t *redshift;

  /* Tomography type */
  tomo_t tomo;

  /* Reduced-shear correction */
  reduced_t reduced;
  double q_mag_size;  /* q_mag_size = 2(alpha+beta-1),         *
		       * alpha, beta: slopes of number density *
		       * with flux (alpha), size (beta)        */

  /* Intrinsic aligmnent */
  ia_t ia;
  double A_ia;           /* IA amplitude */
  ia_terms_t ia_terms;  /* Bit-coded terms, GG=1, GI=2, II=4 */

  /* Halomodel stuff (only initialised if cosmo->nonlinear=halodm) */
  cosmo_hm *hm;

  /* ============================================================ *
   * Precomputed stuff.						  *
   * ============================================================ */

  interTable **g_source;
  interTable **Pshear, **Pg1;

  /* Shear second-order functions */
  interTable **xiP, **xiM, **gamma, **map_gauss, **map_poly;
  double *c_cosebi, psimin_cosebi, psimax_cosebi;
  int N_cosebi;

} cosmo_lens;

typedef struct {
  double r;
  cosmo_lens* self;
} cosmo_lensANDdouble;

typedef struct {
  int i;
  double r;
  cosmo_lens *self;
} cosmo_lensANDintANDdouble;

typedef struct {
  int i, j;
  double r;
  cosmo_lens *self;
} cosmo_lensANDiid;

typedef struct {
  int i, j, t;
  double r;
  cosmo_lens *self;
} cosmo_lensANDiiid;

typedef struct {
  int i_bin, j_bin, pm, n;
  const double *c;
  double thmin;
  cosmo_lens *self;
  error **err;
} cosmo_lensANDextra;

typedef struct {
   int Ntheta, Nzbin;  /* Number of angular and redshift bins */
   int Ntheta2;        /* For combined 2nd and 3rd-order */
   int Nzcorr;         /* Number of z-correlations, Nzcorr=Nzbin*(Nzbin+1)/2 */
   int n;              /* Number of total entries in data vector, n=Ntheta*Nzcorr */
   double *theta;      /* n-dimensional vector of angular scales */
   double *theta2;     /* For angle_range lensformats: (theta,theta2) = (lower,upper) bin limits */ 
   double *data;       /* n-dimensional data vector */
   double *var;        /* n-dimensional vector with variance */
   double *cov[3];     /* Maximum three nxn-dimensional covariance matrix */
   double a1, a2;      /* Coefficients for 'angle_wquadr' */
   double lndetC;
   int usecov;
   lensdata_t type;
   lensformat_t format;
   order_t order;
   decomp_eb_filter_t decomp_eb_filter;
   cov_scaling_t cov_scaling;
   cosmo_lens *fiducial;   /* Needed for ESH09 cov scaling */
} datcov;


/* ============================================================ *
 * Initialisation.						*
 * ============================================================ */

cosmo_lens *init_parameters_lens(double OMEGAM, double OMEGAV, double W0_DE, double W1_DE,
				 double *W_POLY_DE, int N_POLY_DE,
				 double H100, double OMEGAB, double OMEGANUMASS,
				 double NEFFNUMASS, double NORM, double NSPEC,
				 int Nzbin, const int *Nnz, const nofz_t *nofz, double *par_nz,
				 nonlinear_t NONLINEAR, transfer_t TRANSFER,
				 growth_t GROWTH, de_param_t DEPARAM,
				 norm_t normmode, tomo_t TOMO, reduced_t REDUCED, double Q_MAG_SIZE,
				 ia_t IA, ia_terms_t ia_terms, double A_IA, error **err);

void consistency_parameters_lens(const cosmo_lens *self, error **err);
cosmo_lens* copy_parameters_lens_only(cosmo_lens* source, error **err);
cosmo_lens* copy_parameters_lens(cosmo_lens* source, sm2_error **err);
void updateFrom_lens(cosmo_lens* avant, cosmo_lens* apres, error **err);
void copy_parameters_lenshm_cosmo(cosmo_lens *model, error **err);
void read_cosmological_parameters_lens(cosmo_lens **self, FILE *F, error **err);
cosmo_lens* set_cosmological_parameters_to_default_lens(error **err);
void free_parameters_lens(cosmo_lens** self);
void dump_param_lens(cosmo_lens* self, FILE *F, int wnofz, error **err);

/* ============================================================ *
 * Lensing functions.						*
 * ============================================================ */

/* Projection */
double int_for_g(double aprime, void *intpar, error **err);
double g_source(cosmo_lens*, double a, int n_bin, error **err);
double G(cosmo_lens* self, double a, int n_bin, error **err);
double int_for_p_2(double a, void *intpar,error **err);
double P_NL_tot(cosmo_lens *self, double a, double k, error **err);
double Pshear(cosmo_lens *self, double a, int i_bin, int j_bin, error **err);
double P_projected_kappa(void *self, double l, int i_bin, int j_bin, error **err);
double int_over_P_kappa(cosmo_lens *self, funcwithpars int_for_p, void *intpar, error **err);

double int_for_p_GI(double a, void *intpar, error **err);
double int_for_p_II(double a, void *intpar, error **err);


/* Reduced-shear correction (K10) */
extern const int parameter[M_PAR];
cosmo *set_cosmological_parameters_to_WMAP7(const redshift_t *nofz, tomo_t tomo, error **err);
double *par_to_pointer(cosmo *self, par_t par, error **err);
void fill_dpar(cosmo *model, cosmo *wmap7, double *dpar, error **err);
double Fbar(cosmo_lens *self, double a, int m_bin, int n_bin, error **err);
void fill_Fbar_array(cosmo_lens *self, double *fbar, int m_bin, int n_bin, double amin, int N_a,
		     double da, error **err);
void fill_dFbar_dp_array(cosmo_lens *self, par_t par, double *dfbar_dp, int m_bin, int n_bin, double amin,
			 int N_a, double da, error **err);
double Pg1(cosmo_lens *self, double s, int i_bin, int j_bin, error **err);

/* Second-order shear functions */
double xi(cosmo_lens*, int pm, double theta, int i_bin, int j_bin, error **err);
double gamma2(cosmo_lens*, double theta, int i_bin, int j_bin, error **err);
double map2_poly(cosmo_lens*, double theta, int i_bin, int j_bin, error **err);
double map2_gauss(cosmo_lens*, double theta, int i_bin, int j_bin, error **err);
double RR(cosmo_lens *lens, double THETA_MIN, double THETA_MAX, const double *a, int N,
	  poly_t poly, int pm, error **err);
double E_cosebi(cosmo_lens *lens, int n, double Psimin, double Psimax, int i_bin, int j_bin,
		const char *path, double *B_cosebi, error **err);
double RR_cosebi(cosmo_lens *lens, double THETA_MIN, double THETA_MAX, int i_bin, int j_bin,
		 int n, int pm, error **err);
double dRR_cosebi_dz_MC(double *z, int ndim, void *intpar);
double dRR_cosebi_dz(double z, void *intpar, error **err);
double int_for_map2_slow(double ell, void *intpar, error **err);
double map2_slow(cosmo_lens *self, double theta, tpstat_t tpstat, int i_bin, int j_bin, error **err);

/* Reading data files */
datcov *init_data_cov_tomo(char* dataname, char *dataname2, char** covname_ptr, lensdata_t type,
			   decomp_eb_filter_t decomp_eb_filter, 
			   lensformat_t format, double corr_invcov,
			   double a1, double a2, order_t order,
			   cov_scaling_t cov_scaling, error **err);
datcov *init_datcov_for_cov_only(int Nzbin, int Ntheta, error **err);

void del_data_cov(datcov** dc);
void read_data_tomo(datcov *dc, char data_name[], int Nzbin, order_t order, error **err);
void read_cov_tomo(datcov* dc, char cov_name[], int icov, error **err);
void datcov2xipm(const datcov *dc, int i_bin, int j_bin, double **xip, double **xim, double **theta,
		 double **theta2, int *N, error **err);

//void read_cov(datcov* dc, char cov_name[], error **err);
void read_cov_col(datcov *dc, char cov_name[], error **err);

int get_pm(lensdata_t type, int i, int Ntheta, error **err);
int find_bin(double x, const double *list, int N, int prev, error **err);
void scale_cosmic_variance_ESH09(cosmo_lens *model, gsl_matrix *cov, const datcov *dc, error **err);
void scale_mixed_ESH09(const cosmo_lens *model, gsl_matrix *cov, const datcov *dc, error **err);
double lensing_signal(cosmo_lens *model, double theta, int i_bin, int j_bin, lensdata_t type,
		      decomp_eb_filter_t decomp_eb_filter, const cosebi_info_t *cosebi_info, error **err);
double chi2_lensing(cosmo_lens* csm, datcov* dc, int return_model, double **model_array, int *Nmodel,
		    const cosebi_info_t *cosebi_info, error **err);

/* Some third-order stuff which is called from lensing.c */
int Nperm_to_Ntheta(int Nperm, error **err);
void read_data_3rd(datcov *dc, char data_name[], int Nzbin, error **err);
void read_data_2nd_3rd(datcov *res, char *dataname, char *dataname2, error **err);


#define CHANGE(fct) int change_##fct(cosmo_lens*, cosmo_lens*)
CHANGE(g_source);
CHANGE(Pshear);
CHANGE(xi);
CHANGE(gamma2);
CHANGE(map2);
#undef CHANGE


#endif /* __LENSING_H */

