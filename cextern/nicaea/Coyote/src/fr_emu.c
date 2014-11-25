/*
 *  emu.c
 *  Edited and renamed to fr_emu.c for cosmo_pmc by Martin Kilbinger 2013.
 *  
 *
 *  Created by Earl Lawrence on 9/17/09.
 *  Update 10/9/2012
 *
 *  This program was prepared by Los Alamos National Security, LLC at Los Alamos National Laboratory (LANL) 
 *  under contract No. DE-AC52-06NA25396 with the U.S. Department of Energy (DOE). All rights in the program 
 *  are reserved by the DOE and Los Alamos National Security, LLC.  Permission is granted to the public to 
 *  copy and use this software without charge, provided that this Notice and any statement of authorship are 
 *  reproduced on all copies.  Neither the U.S. Government nor LANS makes any warranty, express or implied, 
 *  or assumes any liability or responsibility for the use of this software.  
 *
 *  
 */


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

#include "coyote.h"
#include "fr_constants.h"

/* Copy parameters to emulator input array */
void fill_xstar6_wo_z(double omega_m, double omega_b, double n_spec, double sigma_8, double w0_de, double h_100,
                      double xstar[])
{
   /* Orders reversed:  (omega_b <-> omega_m), (w0_de <-> sigma_8), compared to v1 */
   xstar[0] = omega_b;
   xstar[1] = omega_m;
   xstar[2] = n_spec;
   xstar[3] = h_100 * 100;
   xstar[4] = w0_de;
   xstar[5] = sigma_8;
}


/* ============================================================ *
 * Returns the Coyote non-linear power spectrum. k is in units  *
 * of h/Mpc, h_100 is now an independent parameter (Coyote v2). *
 * It is also used to transform it to [1/Mpc] as                *
 * needed by the Coyote emulator.                               *
 * Six parameters, for coyote13.				*
 * ============================================================ */
double P_NL_coyote6(double omega_m, double omega_b, double n_spec, double sigma_8, double w0_de,
		    double h_100, double a, double k, double **ystar_allz, error **err)
{
   double xstar[p], ystar[2*fr_nsim], val, k_invMpc;
   int ihi, ilo, i;

   /* Copy cosmological parameters */
   fill_xstar6_wo_z(omega_m, omega_b, n_spec, sigma_8, w0_de, h_100, xstar);

   /* Redshift */
   xstar[p] = 1.0/a - 1.0;

   /* Input k is [h/Mpc], Coyote interpolates on k_invMpc in [1/Mpc] */
   k_invMpc  = k * h_100;

   /* Check k-range */
   testErrorRetVA(k_invMpc < fr_ksim[0] || k_invMpc > fr_ksim[fr_nsim-1], coyote_range,
           "Fourier mode k=%g 1/Mpc out of range [%g;%g]",
           *err, __LINE__, -1.0, k_invMpc, fr_ksim[0], fr_ksim[fr_nsim-1]);

   if (*ystar_allz == NULL) {
      *ystar_allz = malloc_err(sizeof(double) * fr_rs * fr_nsim, err);
      forwardError(*err, __LINE__, -1.0);

      /* The actual emulator, filling the array in (k, z) */
      //testErrorRet(ystar_allz == NULL, io_null, "Coyote13: array for (k,z) 'ystar_allz' not initialised",
      //	   *err, __LINE__,);
      fr_fill_ystar_allz(*ystar_allz, xstar, err);
      forwardError(*err, __LINE__, -1.0);
   }

   /* ystar: first half = k, second half = P(k) */
   fr_emu(xstar, ystar, *ystar_allz, err);
   forwardError(*err, __LINE__, -1.0);

   /* Looking for right k-index with bisection */
   ilo = 0;
   ihi = fr_nsim-1;
   while (ihi - ilo > 1) {
      i = (ihi + ilo) >> 1;
      if (ystar[i] > k_invMpc) ihi = i;
      else ilo = i;
   }
   testErrorRetVA(ihi == ilo, math_wrongValue, "Bisection failed, both indices equal (%d)", *err, __LINE__, -1.0, ihi);

   /* Linear interpolation */
   val = (k_invMpc - ystar[ilo]) / (ystar[ihi] - ystar[ilo]) * (ystar[fr_nsim + ihi]
         - ystar[fr_nsim + ilo]) + ystar[fr_nsim + ilo];

   /* Coyote P(k) [Mpc^3] -> output P(k) [(Mpc/h)^3] */
   val *= h_100 * h_100 * h_100;
 
   return val;
}


void fr_check_range(const double *xstar, error **err)
{
   int i;

   for(i=0; i<p; i++) {
      testErrorRetVA((xstar[i] < fr_xmin[i]) || (xstar[i] > fr_xmin[i]+fr_xrange[i]), coyote_range,
		     "Coyote FrankenEmu parameter #%d (%g) out of range [%g; %g]",
		     *err, __LINE__,, i, xstar[i], fr_xmin[i], fr_xmin[i] + fr_xrange[i]);
   } // for(i=0; i<p; i++)

   // Check redshift to make sure we're interpolating
   testErrorRetVA((xstar[p] < 0) || (xstar[p] > 4), coyote_range,
		  "Redshift %g out of range [%g;%g]", *err, __LINE__,, xstar[p], 0.0, 4.0);
}

/* ============================================================ *
 * The actual emulation. Fills the array ystar_allz for (k, z). *
 * Cosmological parameters, placeholder for the output, type of *
 * output.
 * xstar[p] = redshift is not used here, only the p=6 cosmo.    *
 * parameters.							*
 * ============================================================ */
void fr_fill_ystar_allz(double *ystar_allz, const double *xstar, error **err)
{
   const double sd = 0.16002;
   int i, j, k;
   double xstarstd[p], wstar[peta], Sigmastar[peta][m], ystaremu[neta];
   double logc;

   // Interpolation stuff for k
   gsl_spline *lininterp_k = gsl_spline_alloc(gsl_interp_linear, neta/rs);
   gsl_interp_accel *accel = gsl_interp_accel_alloc();


   //fprintf(stderr, "MKDEBUG: Filling ystar_allz for (%g, %g, %g, %g, %g, %g)\n",
   //	   xstar[0], xstar[1], xstar[2], xstar[3], xstar[4], xstar[5]);

    // Standardize the inputs
    for(i=0; i<p; i++) {
        xstarstd[i] = (xstar[i] - fr_xmin[i]) / fr_xrange[i];
    }

    // Compute the covariances between the new input and sims for all PCs
    for(i=0; i<peta; i++) {
        for(j=0; j<m; j++) {
            logc = 0.0;
            for(k=0; k<p; k++) {
                logc -= fr_beta[i][k]*dsqr(fr_x[j][k]-xstarstd[k]);
            }
            Sigmastar[i][j] = exp(logc)/fr_lamz[i];
        }
    }

    // Compute wstar, the predicted PC weights for the new input
    for(i=0; i<peta; i++) {
        wstar[i]=0.0;
        for(j=0; j<m; j++) {
            wstar[i] += Sigmastar[i][j] * fr_KrigBasis[i][j];
        }
    }

    // Compute ystar, the new output
    for(i=0; i<neta; i++) {
        ystaremu[i] = 0.0;
        for(j=0; j<peta; j++) {
            ystaremu[i] += fr_K[i][j]*wstar[j];
        }
        ystaremu[i] = ystaremu[i]*sd + fr_mean[i];
    }

    // Interpolate the emulated output onto the original domain.
    for(i=0; i<rs; i++) {
        gsl_spline_init(lininterp_k, fr_kemu, &ystaremu[i*neta/rs], neta/rs);
        for(j=0; j<fr_nsim; j++) {
            ystar_allz[i*fr_nsim+j] = gsl_spline_eval(lininterp_k, fr_ksim[j], accel);
        }
        gsl_interp_accel_reset(accel);
    }

    gsl_spline_free(lininterp_k);
    gsl_interp_accel_free(accel);
}

/* ============================================================ *
 * This used to be the actual emulation, now fr_fill_ystar_allz *
 * is called from here.						*
 * xstar[p+1] contains the six cosmological parameters plus the *
 * redshift.							*
 * ============================================================ */
void fr_emu(const double *xstar, double *ystar, const double *ystar_allz, error **err)
{
    int i, j;
    //double ystar_allz[rs*fr_nsim];
    double zemu[rs], ybyz[rs];

    /* CosmoPMC version of Coyote: no initialisation, no static variables */

    // Interpolation stuff z
    gsl_spline *lininterp_z = gsl_spline_alloc(gsl_interp_linear, rs);
    gsl_interp_accel *accel = gsl_interp_accel_alloc();
    

    // Check the inputs to make sure we're interpolating.
    fr_check_range(xstar, err);
    forwardError(*err, __LINE__,);

    // Fill in the k values for the final output; first half of ystar
    for(i=0; i<fr_nsim; i++) {
        ystar[i] = fr_ksim[i];
    }

    // Interpolate on to the desired redshift
    // The order needs to be reversed here.
    for(i=0; i<rs; i++) {
        zemu[i] = (1/fr_aemu[rs-i-1]) - 1.0;
    }
    for(i=0; i<fr_nsim; i++) {
        // Build an array with the values of y for a given value of z
        // Reverse the order
        for(j=0; j<rs; j++) {
            ybyz[rs-j-1] = ystar_allz[j*fr_nsim+i];
        }
        gsl_spline_init(lininterp_z, zemu, ybyz, rs);
        ystar[fr_nsim+i] = gsl_spline_eval(lininterp_z, xstar[p], accel);
        gsl_interp_accel_reset(accel);
    }
    
    //printf("z interped\n");
    
    // Transform to P(k); second half of ystar
    for(j=0; j<fr_nsim; j++) {
       ystar[fr_nsim+j] = ystar[fr_nsim+j] - 1.5*log10(fr_ksim[j]);
       ystar[fr_nsim+j] = pow(10.0,ystar[fr_nsim+j])*2.0*M_PI*M_PI;
    }

    // Free some stuff.  I always forget this.  
    // Thanks to Tim Eifler for discovering it (apparently the hard way).
    gsl_spline_free(lininterp_z);
    gsl_interp_accel_free(accel);
}

/* Undefine constants from fr_constants.h */
#undef m
#undef neta
#undef p
#undef peta
#undef rs

