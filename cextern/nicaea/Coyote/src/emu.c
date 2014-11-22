/*
 *  emu.c
 *  
 *
 *  Created by Earl Lawrence on 9/17/09.
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

// Sizes of stuff and number of redshifts

#include "coyote.h"
#include "constants.h"


void KrigBasisOut()
{
   FILE *F;
   int i, j;

   F = fopen("KrigBasis.dat", "w");
   fprintf(F, "const double KrigBasis[peta][m] = {\n");
   for (i=0; i<peta; i++) {
      for (j=0; j<m; j++) {
	 fprintf(F, "% .16f,   ", KrigBasis[i][j]);
      }
      fprintf(F, "\n");
   }
   fprintf(F, "};\n");
   fclose(F);
}

/* Copy parameters to emulator input array */
void fill_xstar_wo_z(double omega_m, double omega_b, double n_spec, double sigma_8, double w0_de, double xstar[])
{
   xstar[0] = omega_m;
   xstar[1] = omega_b;
   xstar[2] = n_spec;
   xstar[3] = sigma_8;
   xstar[4] = w0_de;
}

/* ============================================================ *
 * Returns the Coyote non-linear power spectrum. k is in units  *
 * of h/Mpc, h_100 is needed to transform it to [1/Mpc] as      *
 * needed by the Coyote emulator.				*
 * ============================================================ */
double P_NL_coyote5(double omega_m, double omega_b, double n_spec, double sigma_8, double w0_de,
		    double h_100, double a, double k, error **err)
{
   double xstar[6], ystar[2*nsim], val, k_invMpc;
   int ihi, ilo, i;

   /* Copy cosmological parameters to xstar */
   fill_xstar_wo_z(omega_m, omega_b, n_spec, sigma_8, w0_de, xstar);

   /* Copy redshift */
   xstar[5] = 1.0/a - 1.0;

   /* Input k is [h/Mpc], Coyote interpolates on k_invMpc in [1/Mpc] */
   k_invMpc  = k * h_100;

   /* Check k-range */
   testErrorRetVA(k_invMpc < ksim[0] || k_invMpc > ksim[nsim-1], coyote_range,
   		  "Fourier mode k=%g 1/Mpc out of range [%g;%g]",
   		  *err, __LINE__, -1.0, k_invMpc, ksim[0], ksim[nsim-1]);

   emu(xstar, ystar, err);
   forwardError(*err, __LINE__, -1.0);

   /* Looking for right k-index with bisection */
   ilo = 0;
   ihi = nsim-1;
   while (ihi - ilo > 1) {
      i = (ihi + ilo) >> 1;
      if (ystar[i] > k_invMpc) ihi = i;
      else ilo = i;
   }
   testErrorRetVA(ihi == ilo, math_wrongValue, "Bisection failed, both indices equal (%d)", *err, __LINE__, -1.0, ihi);

   /* Linear interpolation */
   val = (k_invMpc - ystar[ilo]) / (ystar[ihi] - ystar[ilo]) * (ystar[nsim + ihi] - ystar[nsim + ilo]) + ystar[nsim + ilo];

   /* Coyote P(k) [Mpc^3] -> output P(k) [(Mpc/h)^3] */
   val *= h_100 * h_100 * h_100;

   return val;
}

// The actual emulation
// Cosmological parameters, placeholder for the output, type of output
void emu(double *xstar, double *ystar, error **err)
{
    int i, j, k;
    double wstar[peta], Sigmastar[peta][m], ystaremu[neta], ystar_allz[rs*nsim], logc;
    double xstarstd[p];
    double zemu[rs], ybyz[rs];
    // Interpolation stuff for k and then z
    gsl_spline *lininterp_k = gsl_spline_alloc(gsl_interp_linear, neta/rs);
    gsl_spline *lininterp_z = gsl_spline_alloc(gsl_interp_linear, rs);
    gsl_interp_accel *accel = gsl_interp_accel_alloc();
    
    /* CosmoPMC version of Coyote: no initialisation, no static variables */
    
    // Check the inputs to make sure we're interpolating.
    for(i=0; i<p; i++) {
       testErrorRetVA((xstar[i] < xmin[i]) || (xstar[i] > xmin[i]+xrange[i]), coyote_range,
       		      "Coyote emulator parameter #%d (%g) out of range [%g; %g]",
		      *err, __LINE__,, i, xstar[i], xmin[i], xmin[i] + xrange[i]);
    }
    
    // Check redshift to make sure we're interpolating
    testErrorRetVA((xstar[5] < 0.0) || (xstar[5] > 1.0), coyote_range,
		   "Redshift %g out of range [%g;%g]", *err, __LINE__,, xstar[5], 0.0, 1.0);
    
    // Standardize the inputs
    for(i=0; i<p; i++) {
        xstarstd[i] = (xstar[i] - xmin[i]) / xrange[i];
    }
    
    // Compute the covariances between the new input and sims for all PCs
    for(i=0; i<peta; i++) {
        for(j=0; j<m; j++) {
            logc = 0.0;
            for(k=0; k<p; k++) {
                logc -= beta[i][k]*pow(x[j][k]-xstarstd[k], 2.0);
            }
            Sigmastar[i][j] = exp(logc)/lamz[i];
        }
    }
    
    // Compute wstar, the predicted PC weights for the new input
    for(i=0; i<peta; i++) {
        wstar[i]=0.0;
        for(j=0; j<m; j++) {
            wstar[i] += Sigmastar[i][j] * KrigBasis[i][j];
        }
    }
    
    // Compute ystar, the new output
    for(i=0; i<neta; i++) {
        ystaremu[i] = 0.0;
        for(j=0; j<peta; j++) {
            ystaremu[i] += K[i][j]*wstar[j];
        }
        ystaremu[i] = ystaremu[i]*sd + mean[i];
    }
    
    
    // Interpolate the emulated output onto the original domain.
    for(i=0; i<rs; i++) {
        gsl_spline_init(lininterp_k, kemu, &ystaremu[i*neta/rs], neta/rs);
        for(j=0; j<nsim; j++) {
            ystar_allz[i*nsim+j] = gsl_spline_eval(lininterp_k, ksim[j], accel);
        }
        gsl_interp_accel_reset(accel);
    }
    
    // Fill in the k values for the final output
    for(i=0; i<nsim; i++) {
       ystar[i] = ksim[i];
    }
    
    // Interpolate on to the desired redshift
    // The order needs to be reversed here.
    for(i=0; i<rs; i++) {
        zemu[i] = (1/aemu[rs-i-1]) - 1.0;
    }
    for(i=0; i<nsim; i++) {
        // Build an array with the values of y for a given value of z
        // Reverse the order
        for(j=0; j<rs; j++) {
            ybyz[rs-j-1] = ystar_allz[j*nsim+i];
        }
        gsl_spline_init(lininterp_z, zemu, ybyz, rs);
        ystar[nsim+i] = gsl_spline_eval(lininterp_z, xstar[5], accel);
        gsl_interp_accel_reset(accel);
    }
    
    // Transform to P(k)
    for(j=0; j<nsim; j++) {
       ystar[nsim+j] = ystar[nsim+j] - 1.5*log10(ksim[j]);
       ystar[nsim+j] = pow(10.0,ystar[nsim+j])*2.0*M_PI*M_PI;
    }

    // Free some stuff. 
    // Thanks to Tim Eifler for discovering it (apparently the hard way).
    gsl_spline_free(lininterp_z);
    gsl_spline_free(lininterp_k);
    gsl_interp_accel_free(accel);

    //printf("The spectrum is computed using the CMB derived value for h.\n");
    
    //inited=1;
}


#undef m
#undef neta
#undef p
#undef peta
#undef rs
