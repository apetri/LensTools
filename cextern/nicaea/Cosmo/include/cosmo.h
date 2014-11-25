/* ============================================================ *
 * cosmo.h							*
 * Martin Kilbinger, Karim Benabed 2006-2009		        *
 * ============================================================ */

#ifndef __COSMO_H
#define __COSMO_H


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "errorlist.h"
#include "maths.h"
#include "io.h"
#include "config.h"
#include "coyote.h"


/* If fastxi is defined, the integration over comoving distances for the lensing *
 * power spectrum P_kappa done without the need of interpolating P_NL.  	 */
#define fastxi


/* Hubble radius c/H_0 in units of h^{-1} Mpc */
#define R_HUBBLE 2997.92458

/* Physical photon density omega_gamma = Omega_gamma h^2, for T_cmb = 2.725 */
#define omega_gamma 2.469e-5

/* Effective number of massless neutrino species */
#define NEFF 3.04

/* CMB temperature */
#define T_CMB  2.725

/* Dimensions of interpolation tables */
#define _N_a    100 // 50 // 200
/* New! Changed from 100 (nicaea v2.2) -> 250 */
#define N_k     250

/* Ranges of interpolation tables */
#define a_minmin  0.0009
/* [k] = h/Mpc */
#define k_min     3.336e-6
#define k_max     333.6
//ATTENTION k_max is set to 3336.0 for HOD (see halomodel/include/halomodel.h
//#define k_max     3336.0

typedef enum {linear, pd96, smith03, smith03_de, coyote10, coyote13, halodm, smith03_revised} nonlinear_t;
#define snonlinear_t(i) ( \
  i==linear          ? "linear" : \
  i==pd96            ? "pd96"   : \
  i==smith03         ? "smith03" : \
  i==smith03_de      ? "smith03_de" : \
  i==coyote10        ? "coyote10" : \
  i==coyote13        ? "coyote13" : \
  i==halodm          ? "halodm" : \
  i==smith03_revised ? "smith03_revised" : \
  "")
#define Nnonlinear_t 8

typedef enum {bbks, eisenhu, eisenhu_osc, Class, be84} transfer_t;
#define stransfer_t(i) ( \
  i==bbks           ? "bbks" :           \
  i==eisenhu        ? "eisenhu" :	 \
  i==eisenhu_osc    ? "eisenhu_osc" :    \
  i==Class          ? "class" :          \
  i==be84           ? "be84" :		 \
  "")
#define Ntransfer_t 4

typedef enum {heath, growth_de} growth_t;
#define sgrowth_t(i) ( \
  i==heath             ? "heath" : \
  i==growth_de         ? "growth_de" : \
  "")
#define Ngrowth_t 2

typedef enum {jassal, linder, earlyDE, poly_DE} de_param_t;
#define sde_param_t(i) ( \
  i==jassal  ? "jassal" : \
  i==linder  ? "linder" : \
  i==earlyDE ? "earlyDE" : \
  i==poly_DE ? "poly_DE" : \
  "")
#define Nde_param_t 4

typedef enum {norm_s8=0, norm_as=1} norm_t;

/* Cosmology-related errors */
#define ce_none                       0
#define ce_spcpd                    -10
#define ce_nonlin                   -11
#define ce_noknl                    -12
#define ce_zmax                     -13
#define ce_interpol2small           -15
#define ce_interpol2big             -16
#define ce_interpoloutofrangemax    -17
#define ce_alloc                    -18
#define ce_wrongValue               -19
#define ce_tooManySteps             -20
#define ce_badFormat		    -21
#define ce_singularValue            -22
#define ce_interpoloutofrange       -23
#define ce_underflow                -24
#define ce_unknown                  -25
#define ce_negative                 -26
#define ce_overflow                 -27
#define ce_infnan		    -28
#define ce_file			    -29
#define ce_transfer                 -30
#define ce_no_init                  -31
#define ce_nonsquare                -32
#define ce_npar                     -33
#define ce_range                    -34
#define ce_de                       -35
#define ce_omega                    -36

typedef struct {

  /* ============================================================ *
   * Cosmological parameters					  *
   * ============================================================ */

  double Omega_m;		/* matter density parameter                       */
  double Omega_de;		/* dark energy parameter                          */
  double w0_de, w1_de;          /* dark energy eos parametrization		  */
  double *w_poly_de;            /* for w(z) = sum_i w_i a^i                       */
  int    N_poly_de;
  double h_100;                 /* Hubble parameter: H_0 = 100 h_100 km/s/Mpc     */
  double Omega_b;               /* baryon density parameter			  */
  double Omega_nu_mass;         /* density of massive neutrinos                   */
  double Neff_nu_mass;          /* Effective number of massive neutrinos          */  
  double normalization;         /* normalization                                  */
  double sigma_8;		/* power spectrum normalization                   */
  double As;                    /* As						  */
  double n_spec;		/* spectral index of initial power spectrum       */

  /* ============================================================ *
   * Flags.							  *
   * ============================================================ */

  nonlinear_t nonlinear;  /* linear, pd96, smith03, smith03_de, coyote10, coyote13, halodm, smith03_revised */
  transfer_t transfer;    /* bbks, eisenhu, eisenhu_osc	                    */
  growth_t growth;        /* heath, growth_de				    */
  de_param_t de_param;    /* jassal, linder, earlyDE, poly_DE               */
  norm_t normmode;        /* 0 sigma8, 1 As                                 */

  double a_min;           /* Minimum scale factor  */
  int N_a;
 
  /* ============================================================ *
   * Precomputed stuff.						  *
   * ============================================================ */

  interTable* linearGrowth;
  double growth_delta0;

  interTable* transferFct;
  double transfer_alpha_Gamma;
  double transfer_s;
  interTable* transferBE;

  double cmp_sigma8;                  /* sigma8 computed */
  interTable2D* P_NL;
  interTable* slope;
  interTable* w;
   //interTable *k_NL;
  double *ystar_allz; /* For nonlinear=coyote13 */

} cosmo;


typedef struct {
  double r;
  cosmo* self;
} cosmoANDdouble;

typedef struct {
  double r, a;
  cosmo *self;
} cosmoAND2double;

typedef struct {
  int i;
  double r;
  cosmo *self;
} cosmoANDintANDdouble;

typedef struct {
  double r1, r2;
  cosmo *self;
  int i[2];
} cosmoANDdoubleANDdoubleAND2int;

typedef struct {
  int i;
  cosmo *self;
} cosmoANDint;

#define NCLOSE(avant,apres) (fabs(avant-apres)>epsilon1)
#define NCOCLOSE(avant,apres,field) (fabs(avant->field-apres->field)>EPSILON1)
#define NCOEQ(avant,apres,field) (avant->field!=apres->field)


/* ============================================================ *
 * Initialisation.						*
 * ============================================================ */

cosmo* init_parameters(double OMEGAM, double OMEGAV, double W0_DE, double W1_DE,
		       double *W_POLY_DE, int N_POLY_DE,
		       double H100, double OMEGAB, double OMEGANUMASS,
		       double NEFFNUMASS, double NORM, double NSPEC,
		       nonlinear_t NONLINEAR, transfer_t TRANSFER, growth_t GROWTH,
		       de_param_t DEPARAM, norm_t normmode, double AMIN, error **err);

void consistency_parameters(const cosmo *self, error **err);
cosmo* copy_parameters_only(cosmo* source,error **err);
cosmo* copy_parameters(cosmo* source,error **err);
void read_cosmological_parameters(cosmo **self, FILE *F, error **err);
cosmo* set_cosmological_parameters_to_default(error **err);
cosmo *set_cosmological_parameters_to_default2(error **err);

void free_parameters(cosmo** self);
void updateFrom(cosmo* avant, cosmo* apres, error **err);
void updateParameters(cosmo* model, double OMEGAM, double OMEGAV, double W0_DE, double W1_DE,
		      double *W_POLY_DE, int N_POLY_DE,
		      double H100, double OMEGAB, double OMEGANUMASS, double NEFFNUMASS,
		      double NORM, double NSPEC, nonlinear_t NONLINEAR, transfer_t TRANSFER,
		      growth_t GROWTH, de_param_t DEPARAM, norm_t normmode, double AMIN,
		      error **err);

void set_norm(cosmo*, error **err);
void set_w_poly_de(double **w_target, int *N_targer, const double *w_source, int N_source, int check, error **err);

void dump_param(cosmo*, FILE *F);
void dump_param2(cosmo*, FILE *F, char*);

/* ============================================================ *
 * Cosmology.							*
 * ============================================================ */

/* Homogeneous Universe */
double da_dtau(cosmo*,double a,error **err);
double da_dtau_m3(double a, void *intpar,error **err);
double b_early(double w0, double Omegam, double Omegadeinf, error **err);
double w_de(cosmo*, double a, error **err);
double f_de(cosmo*, double a, error **err);
double Esqr(cosmo*, double a, int wOmegar, error **err);
double Omega_m_a(cosmo*,double a, double Esqrpre, error **err);
double Omega_de_a(cosmo*,double a, double Esqrpre, error **err);
void Omega_a(cosmo*,double a, double *omega_m, double *omega_v);
double w_nu_mass(cosmo *self, double a);

/* Geometry, distances */
double int_for_w(double, void *intpar, error **err);
double w(cosmo *, double a, int wOmegar, error **err);
double dwoverda(cosmo *self, double a, error **err);
double drdz(cosmo *self, double a, error **err);
double dvdz(cosmo *self, double a, error **err);
double volume_ana(cosmo *self, double a, error **err);
double f_K(cosmo*, double w, error **err);
double D_lum(cosmo *self, double a, error **err);

/* Growth factor */
void D_plus_derivs(double a, double *y, double *yp, void* extra, error **err);
double D_plus(cosmo*, double a, int normalised, error **err);
double g(cosmo*, double a);

/* Transfer function */
double r_sound(cosmo *model);
double int_for_r_sound(double a, void *intpar, error **err);
double r_sound_integral(cosmo *self, double a, error **err);
double r_sound_drag_fit(cosmo *model, error **err);
double r_sound_drag_analytical(cosmo *self, error **err);
double k_silk(const cosmo *model);
double ratio_b_gamma(cosmo *self, double a);
double z_drag(cosmo *self);
double G_EH98(double y);
double T_tilde(const cosmo *self, double k, double alpha_c, double beta_c);
double Tsqr_one(cosmo*,double k,double Gamma_eff,error **err);
double Tsqr(cosmo*,double k,error **err);

/* Linear power spectrum */
double W_tophat(double x);
double int_for_sigma_R(double logk, void *intpar, error **err);
double sigma_8_sqr(cosmo*, error **err);
double P_L(cosmo* ,double a, double k, error **err);
double P_L_nonorm(cosmo* self, double a, double k, error **err);
double sm2_dfridr(double (*func)(cosmo*,double,double,error **), double x, double h,
                  double *err, double aa, cosmo* self,error **);


/* Peacock&Dodds non-linear power spectrum */
double n_L(cosmo*,double a, double k, error **err);
double f_NL(cosmo*,double x, double a, double k, error **err);

/* Smith et al. non-linear power spectrum (halofit) */
double sm2_transfer(cosmo*, double k, error **err);
double Delta_L_BE2(cosmo*, double k, error **err);
double int_for_wint2_knl(double logk, void *intpar, error **err);
double int_for_wint2_neff(double logk, void *intpar, error **err);
double int_for_wint2_ncur(double logk, void *intpar, error **err);
void wint2(cosmo*,double r,double *sig,double *d1,double *d2, double a, int onlysig,
	   error **err, double precision);
double slope_NL(double rn, double rncur, double om_m, double om_v);
void halofit(double rk, double rn, double rncur, double rknl, double plin,
	     double om_m, double om_v, double *pnl, nonlinear_t nonlinear, double aa, cosmo *self, error **err);
double dlog(double x);
double P_NL(cosmo *self, double a, double k, error **err);
double P_NL_fitting(cosmo*, double a, double k, error **err);
void set_H0_Coyote(cosmo *self, error **err);

/* Coyote emulator interface */
double P_NL_coyote(cosmo *self, double a, double k, error **err);
void set_H0_for_Coyote(cosmo *self, int iterative, error **err);

int test_range_de_conservative(cosmo *model, error **err);

/* Index functions for tomography. Used in lensing & halomodel */
int idx_zz(int i_bin, int j_bin, int Nzbin, error **err);
int idx_zzz(int i_bin, int j_bin, int k_bin, int Nzbin);


/* ============================================================ *
 * Quick and slow direction support				*
 * ============================================================ */
#define CHANGE(fct) int change_##fct(cosmo*, cosmo*)
CHANGE(Esqr);
CHANGE(D_plus);
CHANGE(Tsqr);
CHANGE(sigma_8_sqr);
CHANGE(norm);
CHANGE(Delta_L_BE2);
CHANGE(P_NL);
CHANGE(w);
CHANGE(w_de);
#undef CHANGE

#endif /* __COSMO_H */
