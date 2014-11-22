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
 * =============================================== */

/* Number of k values */
#define nsim 1995
#define sd 1.4613227729637035e-01

/* Constants which are used also outide of Coyote (e.g. in cosmo.c) */
const double ksim[nsim];


/* Functions from hubble.c and emu.c */
double getH0fromCMB(double omega_m, double omega_b, double w0_de, int physical);
void fill_xstar_wo_z(double omega_m, double omega_b, double n_spec, double sigma_8, double w0_de, double xstar[]);
double P_NL_coyote5(double omega_m, double omega_b, double n_spec, double sigma_8, double w0_de,
		    double h_100, double a, double k, error **err);
void emu(double *xstar, double *ystar, error **err);


#endif
