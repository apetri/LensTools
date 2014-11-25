#ifndef __COYOTE_H
#define __COYOTE_H




#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

#include "errorlist.h"
#include "maths.h"


/* Error codes */
#define coyote_base -2200
#define coyote_h        -1 + coyote_base
#define coyote_flat     -2 + coyote_base
#define coyote_range    -3 + coyote_base


/* =============================================== *
 * From other header files which have been removed *
 * Constants which are used also outside of Coyote *
 * (e.g. in cosmo.c) 				   *
 * =============================================== */

/* Coyote I */

/* Number of k values */
#define nsim 1995
extern const double ksim[nsim];

/* Coyote II */

#define fr_nsim 582
#define fr_rs   11
extern const double fr_ksim[fr_nsim];



/* v1: Functions from hubble.c and emu.c */
double getH0fromCMB(double omega_m, double omega_b, double w0_de, int physical);
void fill_xstar_wo_z(double omega_m, double omega_b, double n_spec, double sigma_8, double w0_de, double xstar[]);
double P_NL_coyote5(double omega_m, double omega_b, double n_spec, double sigma_8, double w0_de,
		    double h_100, double a, double k, error **err);
void emu(double *xstar, double *ystar, error **err);

/* v2: Functions from fr_emu.c */
double P_NL_coyote6(double omega_m, double omega_b, double n_spec, double sigma_8, double w0_de,
		    double h_100, double a, double k, double **ystar_allz, error **err);
void fr_check_range(const double *xstar, error **err);
void fr_fill_ystar_allz(double *ystar_allz, const double *xstar, error **err);
void fr_emu(const double *xstar, double *ystar, const double *ystar_allz, error **err);
void fill_xstar6_wo_z(double omega_m, double omega_b, double n_spec, double sigma_8, double w0_de, double h_100,
                      double xstar[]);


#endif
