/* ---------------------------------------------------------------- *
 * hod.c						            *
 * Martin Kilbinger, Henry J. McCracken, Jean Coupon 2008-2013      *
 * Refs:						            *
 *   - Leauthaud et al. 2011 (ApJ, 738, 45)                         *
 *   - Coupon et al. 2012 (A&A 542A, 5)                             *
 *   - Brown et al. 2008 (ApJ 682, 937)		                    * 
 *   - Zheng et al. 2005 (ApJ, 633, 791)         		    *
 *   - Tinker et al. 2005 (ApJ, 631, 41)	                    *
 *   - Kravstov et al. 2004 (ApJ, 609, 35)    	                    *
 *   - Hamana et al 2004 (MNRAS 347, 813)		            *
 *   - Zheng et al. 2005 (ApJ, 633, 791)            	            *
 *   - Berlind & Weinberg 2002 (ApJ 575, 587)			    *
 *                                                                  *
 * WARNINGS:                                                        *
 * - The Hamana model is currently NOT supported                    *
 * - the FFTLog routines have not been tested in case               *
 *   of a change of cosmological parameters during a pmc run.       *
 *   Only tested for changes in HOD parameters.                     *
 * - sm2_romberg should never be used with eps > 1e-4               *
 *   (use int_gsl instead)                                          *
 * ---------------------------------------------------------------- */

/* ---------------------------------------------------------------- *
 * To add a new fitting parameter newly created in this file,
 * here is the procedure:
 *
 * 1. first add parameter in structure "model_hm" in 
 * halomodel/include/halomodel.h
 * 2. then 
 * > make clean
 * > make
 * 3. add parameter in halomodel/src/halomodel.c: 
 * init_parameters_hm(...) 2 places, 
 * copy_parameters_hm_only(...), copy_parameters_hm(...) 
 * and read_cosmological_parameters_hm(...)
 * > make clean
 * > make
 * 4. Cosmo/src/lensing.c: init_parameters_hm(...)
 * 5. update exec/getHODModel.c: initPara(...)
 * 6. add parameter in 
 * tools/include/par.h: 2 places + update #define Npar_t
 * > make clean
 * > make
 * 7. add parameter in wrappers/src/halo.c: fill_parameters_hm(...)
 * > make clean
 * > make
 *
 *                                 and good luck....(`~')
 *
 * ---------------------------------------------------------------- */



#include "hod.h"
double FFTLog_TMP;

double chi2_hm(cosmo_hm *model, halodata_t halodata, halomode_t halomode, const wt_t* data, ngal_fit_t ngal_fit_type,
	       double ngal, double ngal_err, intconst_t intconst_type, error **err)
/* ---------------------------------------------------------------- *
 * Returns the chi^2 for w(theta), wp(rp) and the 
 * galaxy number density.
 * ---------------------------------------------------------------- */
{
  
  int i, j;
  double chi2, ngal_model, zm;
    
  /* to do */
  testErrorRet(model->redshift->Nzbin!=1, ce_overflow,
	       "More than one redshift bin not yet supported in likelihood",
	       *err, __LINE__, -1.0);
  
  /* ngal fit not used in leauthaud11 model */
  testErrorExit(model->hod == leauthaud11 && ngal_fit_type != ngal_no_fit, HOD_NGAL_FIT,
	       "ngal fit cannot be used with leauthaud11 mode; set \"sngal_fit_type ngal_no_fit\" and use stellar mass function instead",
	       *err, __LINE__);
  

  /* mean redshift */
  zm = zmean(model->redshift,0, err);
  forwardError(*err, __LINE__, -1.0);
  
  /* galaxy number density in h^3 Mpc^{-3} */
  ngal_model = ngal_den(model,1.0/(1.0+zm), logMmax, model->Mstellar_min, model->Mstellar_max, err);   
  forwardError(*err, __LINE__, -1.0);
  
  /* if ngal only is used in chi2 */
  if(ngal_fit_type == ngal_lin_fit_only){
    chi2 = dsqr((ngal - ngal_model)/ngal_err);
    testErrorRetVA(chi2<0.0, math_negative, "Negative chi^2 %g",*err, __LINE__, -1.0, chi2);
    testErrorRet(!finite(chi2), ce_infnan, "inf or nan chi^2", *err, __LINE__, -1.0);
    return -0.5*chi2;
  }
  
  /* w(theta), wp(rp) */
  double *x, *y_model, IC, RR_total;
  double *ystellar, *y1hgcs, *y1hgss, *y2hg; 
  switch (halodata) {
  case w_of_theta :
    y_model = malloc_err(sizeof(double)*data->nbins, err);  forwardError(*err, __LINE__, -1.0);
    y1hgcs  = woftheta(model, p1hgcs, data->th, data->nbins, 0, 0, err); forwardError(*err, __LINE__, -1.0);
    y1hgss  = woftheta(model, p1hgss, data->th, data->nbins, 0, 0, err); forwardError(*err, __LINE__, -1.0);
    y2hg    = woftheta(model, p2hg,   data->th, data->nbins, 0, 0, err); forwardError(*err, __LINE__, -1.0);
    for (i=0; i< data->nbins; i++) y_model[i] = y1hgcs[i] + y1hgss[i] + y2hg[i];
    free(y1hgcs); free(y1hgss); free(y2hg);
    
    if(intconst_type == random_file){
      y1hgcs  = woftheta(model, p1hgcs, data->th_RR, data->nbins_RR, 0, 0, err); forwardError(*err, __LINE__, -1.0);
      y1hgss  = woftheta(model, p1hgss, data->th_RR, data->nbins_RR, 0, 0, err); forwardError(*err, __LINE__, -1.0);
      y2hg    = woftheta(model, p2hg,   data->th_RR, data->nbins_RR, 0, 0, err); forwardError(*err, __LINE__, -1.0);
      IC       = 0.0;
      RR_total = 0.0;
      for (i=0; i< data->nbins_RR; i++){
	IC       += (y1hgcs[i] + y1hgss[i] + y2hg[i])*data->RR[i];
	RR_total += data->RR[i];
      }
      for (i=0; i< data->nbins; i++) y_model[i] -= IC/RR_total;
      free(y1hgcs); free(y1hgss); free(y2hg);
    }
    break;
  case wp_rp :
    y_model = malloc_err(sizeof(double)*data->nbins, err);  forwardError(*err, __LINE__, -1.0);
    y1hgcs = wp(model, p1hgcs, data->th, data->nbins, model->pi_max, GG, err);  forwardError(*err, __LINE__, -1.0);
    y1hgss = wp(model, p1hgss, data->th, data->nbins, model->pi_max, GG, err);  forwardError(*err, __LINE__, -1.0);
    y2hg   = wp(model, p2hg,   data->th, data->nbins, model->pi_max, GG, err);  forwardError(*err, __LINE__, -1.0);
    for (i=0; i< data->nbins; i++) y_model[i] = y1hgcs[i] + y1hgss[i] + y2hg[i];
    free(y1hgcs); free(y1hgss); free(y2hg);
    break;
  case deltaSigma :
    y_model = malloc_err(sizeof(double)*data->nbins, err);  forwardError(*err, __LINE__, -1.0);
    ystellar = DeltaSigma(model, pstellar, data->th, data->nbins, err);  forwardError(*err, __LINE__, -1.0);
    y1hgcs   = DeltaSigma(model, p1hgcs,   data->th, data->nbins, err);  forwardError(*err, __LINE__, -1.0);
    y1hgss   = DeltaSigma(model, p1hgss,   data->th, data->nbins, err);  forwardError(*err, __LINE__, -1.0);
    y2hg     = DeltaSigma(model, p2hg,     data->th, data->nbins, err);  forwardError(*err, __LINE__, -1.0);
    for (i=0; i< data->nbins; i++) y_model[i] = ystellar[i] + y1hgcs[i] + y1hgss[i] + y2hg[i];
    free(ystellar); free(y1hgcs); free(y1hgss); free(y2hg);
    break;
  case smf :
    x = malloc_err(sizeof(double)*data->nbins, err);  forwardError(*err, __LINE__, -1.0);
    for(i=0;i<data->nbins;i++) x[i]  = log10(data->th[i]);
    y_model = dNdlogM10stellar(model, x, data->nbins, err);  forwardError(*err, __LINE__, -1.0);
    free(x);
    break;
  }
  
  /* log(w(theta)) */
  /* [Jean]: useless. To be deprecated */
  if (halomode == galcorr_log) for (i=0; i< data->nbins; i++) y_model[i] = log(y_model[i]);
  
  /* chi2 */
  chi2 = 0.0;
  if (data->wcov == NULL){
    for (i=0; i < data->nbins; i++) chi2 += dsqr((data->w[i] - y_model[i])/(data->werr[i]));
  }else{
    for (i=0; i < data->nbins; i++) {
      for (j=0; j < data->nbins; j++) {
	chi2 += (data->w[i] - y_model[i])*data->wcov[i*data->nbins+j]*(data->w[j] - y_model[j]);
	testErrorRet(!finite(data->wcov[i*data->nbins+j]), ce_infnan, "inf or nan in logl", *err, __LINE__, -1.0);
      }
    }
  }
  free(y_model);
  
  /* Number density of galaxies */
  switch (ngal_fit_type) {
  case ngal_log_fit : /* "Number density varies logarithmically with M1" [Hamana 2004] */
    chi2 += dsqr((log(ngal) - log(ngal_model))/(ngal_err/ngal)); /* Delta ln x = Delta x / x */
    break;
  case ngal_lin_fit :
    chi2 += dsqr((ngal - ngal_model)/ngal_err);
    break;
  case ngal_no_fit :
    break;
    /* DEBUGGING
       JEAN: removed (useless) */
    /* case ngal_match_M1 : */
    /* MKDEBUG: TODO */
  default :
    *err = addErrorVA(ce_unknown, "Wrong ngal_fit_type %d", *err, __LINE__, (int)ngal_fit_type);
    return -1.0;
  }

  /* DEBUGGING [Jean]: useless
  if(ngal_fit_type != ngal_lin_fit_only){
    fprintf(stderr, "chi2_hm: ngden(mod,obs,err) = (%5.2e,%5.2e,%5.2e) ln(ngal_model)(mod,obs,err) = (%g,%g,%g, %g)\n",
	    ngal_model, ngal, ngal_err, log(ngal_model), log(ngal), ngal_err/ngal, log(ngal_err+ngal) - log(ngal-ngal_err));
  }
  */
  
  testErrorRetVA(chi2<0.0, math_negative, 
		 "Negative chi^2 %g. Maybe the covariance matrix is not positive",*err, __LINE__, -1.0, chi2);
  testErrorRet(!finite(chi2), ce_infnan, "inf or nan chi2", *err, __LINE__, -1.0);

  printf("chi2 = %f\n", chi2);
  
  /* det C ... */
  return -0.5*chi2;
}

/* ---------------------------------------------------------------- *
 * w(theta)                                                         *
 * ---------------------------------------------------------------- */

double *woftheta(cosmo_hm *model, pofk_t pofk, double *theta, int Ntheta, int i_bin, int j_bin, error **err)

/* ---------------------------------------------------------------- *
 * First compute xi, project using Limber equation at mean z        *
 * and return w(theta).                                             *
 * See Bartelmann & Schneider (2001) eq. (2.79),                    *
 * Tinker et al. 2010 eq. (3), Ross et al. 2009 eq (27), ...        *
 * Theta in degree.                                                 *
 * N is the number of points that sample xi.The speed of the code   *
 * is basically inversely proportional to this number...choose      *
 * it carefuly!! You can also play with N in FFTLog routines.       *
 * ---------------------------------------------------------------- */
{
  testErrorRetVA(MIN(i_bin, j_bin)<0 || MAX(i_bin, j_bin)>=model->redshift->Nzbin, redshift_Nzbin,
		 "Requested redshift bins (%d, %d) out of range [0; %d]",
		 *err, __LINE__, NULL, i_bin, j_bin, model->redshift->Nzbin-1);
  
  double nsqr_dzdr;
  double *result  = malloc_err(Ntheta*sizeof(double),err); 
  forwardError(*err, __LINE__, NULL);
  
  /* tabulate xi(r) */
  int i,j,k,N      = 40;
  double *u        = malloc_err(N*sizeof(double), err);  forwardError(*err, __LINE__, NULL);
  double *logu     = malloc_err(N*sizeof(double), err);  forwardError(*err, __LINE__, NULL);
  double umin      = 0.001, umax = 800;
  double dlogu     = log(umax/umin)/(double)N;
  
  for(i=0;i<N;i++){
    logu[i] = log(umin)+dlogu*(double)i;
    u[i]    = exp(logu[i]);
  }
  
  double *xi = xiofr(model, pofk, u, N, GG, err);
  forwardError(*err, __LINE__, NULL);
  
  /* interpolate xi(r) */
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline    = gsl_spline_alloc (gsl_interp_cspline,N);
  gsl_spline_init(spline,logu,xi,N);
  
  Nz_hist *nz1, *nz2;
  nz1 = get_nz(model->redshift, i_bin, err);   forwardError(*err, __LINE__, NULL);
  nz2 = get_nz(model->redshift, j_bin, err);   forwardError(*err, __LINE__, NULL);

  testErrorRetVA(nz1->nbins != nz2->nbins, redshift_Nzbin,
		 "Number of histogram bins has to be equal for all redshift bins. Found: (%d, %d) for zbins (%d, %d)",
		 *err, __LINE__, NULL, nz1->nbins, nz2->nbins, i_bin, j_bin);

  testErrorRetVA(nz1->dz != nz2->dz, redshift_dz,
		 "Histogram dz-values have to be the same for all redshift bins. Found (%g, %g) for zbins (%d, %d)",
		 *err, __LINE__, NULL, j, nz1->dz, nz2->dz, i_bin, j_bin);

  int wOmegar           = 0;
  double deg_to_rad_sqr = dsqr(pi/180.0);
  double ww, r, x, sum, a;
  
  /* Limber equation - project xi to get w */
  for(i=0;i<Ntheta;i++){ /* loop over theta */
    result[i] = 0.0;
    
    /* loop over z */
    for (j=0; j<nz1->nbins; j++) { 
      testErrorRetVA(nz1->z[j] != nz2->z[j], redshift_dz,
		     "Histogram z-values have to be the same for all redshift bins. Found z[%d]=(%g, %g) for zbins (%d, %d)",
		     *err, __LINE__, NULL, j, nz1->z[j], nz2->z[j], i_bin, j_bin);
      
      a   = 1.0/(1.0 + nz1->z[j]);
      ww  = w(model->cosmo, a, wOmegar, err);                 forwardError(*err, __LINE__, NULL);
      x   = f_K(model->cosmo, ww, err);                       forwardError(*err, __LINE__, NULL);
      sum = 0.0;
      
      /* loop over u */
      for (k=0;k<N;k++) { 
	r    = sqrt(u[k]*u[k] + x*x*theta[i]*theta[i]*deg_to_rad_sqr);
	if (log(r) < logu[N-1]) {
	  /* to unsure log(r) lies within interpolation limits */
	  sum += u[k]*gsl_spline_eval(spline,log(r),acc);
	}
      }
      nsqr_dzdr = nz1->n[j] * nz2->n[j] / drdz(model->cosmo, a, err);
      forwardError(*err, __LINE__, NULL);
      
      result[i] += nsqr_dzdr * sum;
      
    }

    result[i] *= 2.0 * nz1->dz * dlogu;
    testErrorRetVA(!finite(result[i]), ce_infnan, "inf or nan in w(theta_%d=%g)",
		   *err, __LINE__, NULL, i, theta[i]/arcmin);
  }
  
  free(xi);
  free(u);
  free(logu);
  free_nz(nz1); free_nz(nz2);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  
  return result;
}

/* Returns an array of w(theta) for the sum of the three HOD terms (1hcs + 1hss + 2h) */
double* woftheta_FFTLog_total(cosmo_hm *model, double ln_theta_min, double ln_delta_theta, int Ntheta,
			      double **theta, int i_bin, int j_bin, error **err)
{
   int i;
   double *wh, *w1hgcs, *w1hgss, *w2hg;

   wh    = malloc_err(sizeof(double)*Ntheta, err);  forwardError(*err, __LINE__, NULL);

   if (*theta == NULL) {
      printf("Alloc\n");
      *theta = malloc_err(sizeof(double)*Ntheta, err);  forwardError(*err, __LINE__, NULL);
      for (i=0; i<Ntheta; i++) {
         (*theta)[i] = exp(ln_theta_min + i*ln_delta_theta);
      }
   } else {
      printf("No alloc\n");
   }

   w1hgcs = woftheta(model, p1hgcs, *theta, Ntheta, i_bin, j_bin, err);
   forwardError(*err, __LINE__, NULL);

   w1hgss = woftheta(model, p1hgss, *theta, Ntheta, i_bin, j_bin, err);
   forwardError(*err, __LINE__, NULL);

   w2hg   = woftheta(model, p2hg, *theta, Ntheta, i_bin, j_bin, err);
   forwardError(*err, __LINE__, NULL);

  for (i=0; i<Ntheta; i++) wh[i] = w1hgcs[i] + w1hgss[i] + w2hg[i];

  free(w1hgcs); free(w1hgss); free(w2hg);

  return wh;
}


/*----------------------------------------------------------------*
 *wp(rp)                                                          *
 *----------------------------------------------------------------*/
double *wp(cosmo_hm *model, pofk_t pofk, const double *rp, int Nrp, double pi_max, int type, error **err)
/* ---------------------------------------------------------------- *
 * Computes xi, projects it along pi and returns wp.                * 
 * See e.g. Zehavi et al. (2005) Eq. (3).                           *
 * rp in mpc.                                                       *
 * ---------------------------------------------------------------- */
{
  
  /* Coordinates (como or phys) TO DO: add this as parameter */
  double fac = 1.0, zm = zmean(model->redshift, 0, err); forwardError(*err, __LINE__, NULL);
  
  model->coord_phys = 1;
  if(model->coord_phys) fac = 1.0 + zm;

  /* rmax should be smaller than pi_max */
  testErrorExit(rp[Nrp-1] > pi_max, GM_pi_max,  
		"rmax(physical) should be smaller than pi_max", *err, __LINE__);
  
  
  /* interpolate to speed up integration  */
  int i, Ninter      = 40;
  double *logrinter  = malloc_err(Ninter*sizeof(double), err); forwardError(*err, __LINE__, NULL);
  double *rinter     = malloc_err(Ninter*sizeof(double), err); forwardError(*err, __LINE__, NULL);
  double rinter_max  = RMAX;
  double rinter_min  = RMIN;
  double dlogrinter  = log(rinter_max/rinter_min)/(double)Ninter;

  for(i=0;i<Ninter;i++){
    logrinter[i] = log(rinter_min)+dlogrinter*(double)i;
    rinter[i]    = exp(logrinter[i]);
  }  
  
  double *xi = xiofr(model, pofk, rinter, Ninter, type, err);
  forwardError(*err, __LINE__, NULL);
  
  /* interpolate xi(r) */
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline    = gsl_spline_alloc (gsl_interp_cspline, Ninter);
  gsl_spline_init(spline, logrinter, xi, Ninter);
  
  cosmo_hm_params params;
  params.acc       = acc;
  params.spline    = spline;
  params.eps       = 1.0e-4;
  params.logrmin   = logrinter[0];
  params.logrmax   = logrinter[Ninter-1];
  
  double *result   = malloc_err(Nrp*sizeof(double), err); forwardError(*err, __LINE__, NULL);
  
  for(i=0;i<Nrp;i++){
    params.rp = rp[i]*fac;
    if(params.logrmin < log(rp[i]*fac) && log(rp[i]*fac) < params.logrmax){
      result[i] = 2.0*int_gsl(int_for_wp, (void*)&params, log(rp[i]*fac), log(pi_max*fac), params.eps, err)/fac; 
      forwardError(*err, __LINE__, NULL);
    }else{
      result[i] = 0.0;
    }      
  }
  

  free(xi);
  free(rinter);
  free(logrinter);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);

  return result;
}

double int_for_wp(double logr, void *params, error **err)
{
  double result         = 0.0;
  double r              = exp(logr);
  double logrmin        = ((cosmo_hm_params *)params)->logrmin;
  double logrmax        = ((cosmo_hm_params *)params)->logrmax;
  double rp             = ((cosmo_hm_params *)params)->rp;
  gsl_interp_accel *acc = ((cosmo_hm_params *)params)->acc;
  gsl_spline *spline    = ((cosmo_hm_params *)params)->spline;
  
  if(logrmin < logr && logr < logrmax){
    result = r*r*gsl_spline_eval(spline, logr, acc)/sqrt(r*r - rp*rp);
  }
  
  return result;
}


/* ---------------------------------------------------------------- *
 * Lensing functions (only for leauthaud11 model)                   *
 * ---------------------------------------------------------------- */

double *DeltaSigma(cosmo_hm *model, pofk_t pofk, const double *r, int N, error **err)

/* ---------------------------------------------------------------- *
 * Computes DeltaSigma = Sigma( < r) - Sigma(r)                     *
 * See Yoo et al. (2006), Leauthaud et al. (2011)                   *
 * In [h M_sun pc^-2]                                               *
 * This is gamma X Sigma_crit X 1e-12                               *
 * ---------------------------------------------------------------- */
{
  
  int i;
  
  /* DeltaSigma in physical coordinates */
  double zm = zmean(model->redshift, 0, err); forwardError(*err, __LINE__, NULL);
  double  Mstellar_mean, fac = 1.0;
  
  /* Coordinates (como or phys) TO DO: add this as parameter */
  model->coord_phys = 1;
  if(model->coord_phys) fac = 1.0 + zm;
  
  /* rmax should be smaller than pi_max */
  testErrorExit(r[N-1] > model->pi_max, GM_pi_max,  
		"rmax(physical) should be smaller than pi_max", *err, __LINE__);
 
  double *result = malloc_err(N*sizeof(double), err); forwardError(*err, __LINE__, NULL);
  
  if(pofk == pstellar){
    if (model->hod == leauthaud11){
      Mstellar_mean = mass_weighted_av_stellar_mass(model, 1.0/(1.0+zm), err);
    }else{
      /* "berwein" models don't know about SMF, so we take the mean Mstellar.
	 If the sample only has a Mstallar threshold (Mstellar_max = -1),
	 calculate <Mstar> as Mstellar_min.
	 If both min and max undefind (-1), DeltaSigma_stellar = 0.
      */
      if (model->Mstellar_min < 0 || model->Mstellar_max < 0) {
	Mstellar_mean = 0;
      } else {
	Mstellar_mean = (MAX(model->Mstellar_max, model->Mstellar_min) + model->Mstellar_min)/2.0;
      }
    }
    
    /* Stellar mass as point source, Leauthaud et al. (2011) eq. (37) */
    for(i=0; i<N; i++) result[i] = 1.0e-12 * Mstellar_mean/(pi*dsqr(r[i]));
    
  }else{
    

    /* interpolate to speed up integration  */
    int Ninter         = 40;
    double *logrinter  = malloc_err(Ninter*sizeof(double), err); forwardError(*err, __LINE__, NULL);
    double *rinter     = malloc_err(Ninter*sizeof(double), err); forwardError(*err, __LINE__, NULL);
    double rinter_max  = model->pi_max;
    double rinter_min  = MIN(RMIN1, r[0]);
    double dlogrinter  = log(rinter_max/rinter_min)/(double)Ninter;
    
    for(i=0;i<Ninter;i++){
      logrinter[i] = log(rinter_min)+dlogrinter*(double)i;
      rinter[i]    = exp(logrinter[i]);
    }
    
    /* projected density contrast (see Yoo et al. (2006) Eq. 2 */
    double *Sigma = wp(model, pofk, rinter, Ninter, model->pi_max,  GM, err); forwardError(*err, __LINE__, NULL);
    
    /* Interpolation */
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline    = gsl_spline_alloc (gsl_interp_cspline, Ninter);
    gsl_spline_init(spline, logrinter, Sigma, Ninter);
    
    cosmo_hm_params params;
    params.acc       = acc;
    params.spline    = spline;
    params.eps       = 1.0e-4;
    params.logrmin   = logrinter[0];
    params.logrmax   = logrinter[Ninter-1];
    
    /* Those are Omega_m_0 and rho_c0. Functions used only for debugging */
    double Omega_m = Omega_m_halo(model, 1.0/(1.0+zm), err);  forwardError(*err, __LINE__, NULL);
    double rhocrit = rho_crit_halo(model, 1.0/(1.0+zm), err); forwardError(*err, __LINE__, NULL);
    /* Loop over input scales */
    for(i=0; i<N; i++){
      if(params.logrmin < log(r[i]) && log(r[i]) < params.logrmax){
	result[i] = 1.0e-12*(2.0/dsqr(r[i])
			     * int_gsl(int_for_Sigma, (void*)&params, log(1.0e-6), log(r[i]), params.eps, err)
			     - gsl_spline_eval(spline, log(r[i]), acc));
	result[i] *= rhocrit*Omega_m*fac*fac*fac;
	forwardError(*err, __LINE__, NULL);
      }else{
	result[i] = 0.0;
      }      
    }
    
    free(Sigma);
    free(logrinter);
    free(rinter);
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
  }
  
  return result;
}

double int_for_Sigma(double logr, void *params, error **err)
{
  double result         = 0.0;
  double r              = exp(logr);
  double logrmin        = ((cosmo_hm_params *)params)->logrmin;
  double logrmax        = ((cosmo_hm_params *)params)->logrmax;
  gsl_interp_accel *acc = ((cosmo_hm_params *)params)->acc;
  gsl_spline *spline    = ((cosmo_hm_params *)params)->spline;
  
  if(logrmin < logr && logr < logrmax){
    result = r*r*gsl_spline_eval(spline, logr, acc);
  }
  
  return result;
}


/* ---------------------------------------------------------------- *
 * Stellar mass function (only for leauthaud11 model)               *
 * ---------------------------------------------------------------- */

double *dNdlogM10stellar(cosmo_hm *model, double *log10Mstellar, int N, error **err)
{
  /* Computes the stellar mass function. See Leauthaud et al. (2011) eq. 17 
   * and Coupon et al. (2012) Eq. A. 21.
   * In h^3 Mpc^-3 dex^-1
   */
  
  int i;
  
  double *result = malloc_err(N*sizeof(double), err); 
  forwardError(*err, __LINE__, NULL);
  
  double zm = zmean(model->redshift, 0, err);
  forwardError(*err, __LINE__, NULL);
  
  double dlog10Mstellar = 0.01;
  
  for(i=0; i<N; i++){
    
    result[i] = ngal_den(model, 1.0/(1.0+zm), logMmax, 
			 pow(10.0, log10Mstellar[i] - dlog10Mstellar/2.0), 
			 pow(10.0, log10Mstellar[i] + dlog10Mstellar/2.0), 
			 err)/dlog10Mstellar;
    forwardError(*err, __LINE__, NULL);

  }
  
  return result;
}


double *dNdlogM10stellar_c(cosmo_hm *model, double *log10Mstellar, int N, error **err)
{
  /* Computes the stellar mass function. See Leauthaud et al. (2011) eq. 17 
   * and Coupon et al. (2012) Eq. A. 21.
   * In h^3 Mpc^-3 dex^-1
   */
  
  int i;
  
  double *result = malloc_err(N*sizeof(double), err); 
  forwardError(*err, __LINE__, NULL);
  
  double zm = zmean(model->redshift, 0, err);
  forwardError(*err, __LINE__, NULL);
  
  double dlog10Mstellar = 0.01;
  
  for(i=0; i<N; i++){
    result[i] = ngal_den_c(model, 1.0/(1.0+zm), logMmax, 
			 pow(10.0, log10Mstellar[i] - dlog10Mstellar/2.0), 
			 pow(10.0, log10Mstellar[i] + dlog10Mstellar/2.0), 
			   err)/dlog10Mstellar;
  }
  
  return result;
}


double *dNdlogM10stellar_s(cosmo_hm *model, double *log10Mstellar, int N, error **err)
{
  /* Computes the stellar mass function. See Leauthaud et al. (2011) eq. 17 
   * and Coupon et al. (2012) Eq. A. 21.
   * In h^3 Mpc^-3 dex^-1
   */
  
  int i;
  
  double *result = malloc_err(N*sizeof(double), err); 
  forwardError(*err, __LINE__, NULL);
  
  double zm = zmean(model->redshift, 0, err);
  forwardError(*err, __LINE__, NULL);
  
  double dlog10Mstellar = 0.01;
  
  for(i=0; i<N; i++){
    result[i] = ngal_den_s(model, 1.0/(1.0+zm), logMmax, 
			 pow(10.0, log10Mstellar[i] - dlog10Mstellar/2.0), 
			 pow(10.0, log10Mstellar[i] + dlog10Mstellar/2.0), 
			 err)/dlog10Mstellar;
  }
  
  return result;
}


/*----------------------------------------------------------------*
 *HOD P(k) and xi(r)                                              *
 *----------------------------------------------------------------*/

double *xiofr(cosmo_hm *model, pofk_t pofk, const double *r, int N, int type, error **err){
  /*Compute xi(r) at mean redshift*/
  

   double *res;
   double a = 1.0/(1.0 + zmean(model->redshift, 0, err));
   forwardError(*err, __LINE__, NULL);
   
   /* q = 0.0, m = 1/2 */
   switch (pofk) {
   case p1hgcs : /* 1-halo term: central-satellite term */
     res = xi_1hcs(model, a, r, N, type, err);
     forwardError(*err, __LINE__, NULL);
     break;
   case p1hgss : /* 1-halo term: satellite-satellite term */
     res = xi_1hss(model, a, r, N, type, err);
     forwardError(*err, __LINE__, NULL);
     break;
   case p2hg   : /* 2-halo term */ 
     res = xi_2h(model, a, r, N, type, err);
     forwardError(*err, __LINE__, NULL);
     break;
   case pnl    : /* This is xi dark matter. 'type' is not taken into account */
     res = xi_P_NL(model, a, r, N, err);
     forwardError(*err, __LINE__, NULL);
     break;
   default:
     *err = addErrorVA(hm_pofk, "Wrong pofk type %, has to be one of %d (p1hgcs), %d (p1hgss) or %d (p2hg)",
		       *err, __LINE__, pofk, p1hgcs, p1hgss, p2hg);
     return NULL;
   }
   
   return res;
}

double *xi_1hcs(cosmo_hm *model,double a, const double *r, int N, int type, error **err){
  /* 1-halo central-satellite */

  int i;
  double *result = malloc_err(N*sizeof(double), err); forwardError(*err, __LINE__, NULL);
  
  cosmo_hm_params params;
  params.cosmo  = model->cosmo;
  params.model  = model;
  params.a      = a;
  params.ng     = ngal_den(model, a, logMmax, model->Mstellar_min, model->Mstellar_max, err);
  forwardError(*err, __LINE__, NULL);
  params.err    = err;
  params.eps    = 1.0e-4;
  
  for(i=0;i<N;i++){
    if(r[i] > RMAX1){
      result[i] = 0.0;
    }else{
      params.r = r[i];
      params.type = type;
      result[i] = int_gsl(int_for_xi_1hcs, (void*)&params, log(M_vir(model, r[i], a, err)),
			  logMmax, params.eps, err);
      forwardError(*err, __LINE__, NULL);
    }
  }
  return result;
}

double int_for_xi_1hcs(double logM, void *params, error **err)
/* Type can be:
 * GG: galaxy-galaxy stuff,
 * GM: galaxy-Dark matter stuff (for lensing).
 */
{
  double c_sat, rho, res, Omega_m, rhocrit;
  double M         = exp(logM);
  cosmo_hm *model  =  ((cosmo_hm_params *)params)->model;
  double a         =  ((cosmo_hm_params *)params)->a; 
  double r         =  ((cosmo_hm_params *)params)->r;
  double ng        =  ((cosmo_hm_params *)params)->ng;
  int type         =  ((cosmo_hm_params *)params)->type;
  
  
  /* Hala mass function parameters */
  ((cosmo_hm_params *)params)->asymptotic = 0;
  
  res  = dn_dlnM(M, params, err);  forwardError(*err, __LINE__, 0.0);
  
  switch(type){
  case GG:
    
    c_sat = concentration_sat(model, M, a, err);
    forwardError(*err, __LINE__, 0.0);
    
    rho = rho_halo(model, r, a, M, c_sat, err);
    forwardError(*err, __LINE__, 0.0);
    
    res *= Ngal_c(model, M, model->Mstellar_min, model->Mstellar_max)*
      Ngal_s(model, M, model->Mstellar_min, model->Mstellar_max)/(0.5*ng*ng)*rho/M;
    break;
  case GM:
    
    rho = rho_halo(model, r, a, M, -1.0, err);
    forwardError(*err, __LINE__, 0.0);
    
    Omega_m = Omega_m_halo(model, a, err);  forwardError(*err, __LINE__, 0.0);
    rhocrit = rho_crit_halo(model, a, err); forwardError(*err, __LINE__, 0.0);
    
    res  *= Ngal_c(model, M, model->Mstellar_min, model->Mstellar_max)/ng * rho/(rhocrit*Omega_m);
    forwardError(*err, __LINE__, 0.0);
    break;
  }
  
  return res;
}

double *xi_1hss(cosmo_hm *model,double a, const double *r, int N, int type, error **err){
  /* 1-halo satellite-satellite */

  int i,j;
  double *result = malloc_err(N*sizeof(double), err);  forwardError(*err, __LINE__, NULL);
  
  cosmo_hm_params params;
  params.cosmo = model->cosmo;
  params.model = model;
  params.a     = a;
  params.ng    = ngal_den(model, a, logMmax, model->Mstellar_min, model->Mstellar_max, err);
  forwardError(*err, __LINE__, NULL);
  params.err   = err;
  params.eps   = 1.0e-4;
  params.type  = type;
  
  FFTLog_config *fc = FFTLog_init(128, k_min, k_max_HOD, 0.0, 0.5);
  
  double *r_FFT    = malloc_err(fc->N*sizeof(double), err);  forwardError(*err, __LINE__, NULL);
  double *ar       = malloc_err(fc->N*sizeof(double), err);  forwardError(*err, __LINE__, NULL);
  double *logr_FFT = malloc_err(fc->N*sizeof(double), err);  forwardError(*err, __LINE__, NULL);
  
  gsl_function Pk;
  Pk.function = &P1hss;
  Pk.params   = &params;
  
  /* FFTlog */
  FFTLog(fc, &Pk, r_FFT, ar, -1, err); forwardError(*err, __LINE__, NULL);
 
  /* Interpolation */
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline    = gsl_spline_alloc (gsl_interp_cspline,fc->N);
  
  /* Attention: N and fc->N are different */
  for(j=0;j<fc->N;j++) logr_FFT[j] = log(r_FFT[j]);
  gsl_spline_init (spline,logr_FFT, ar, fc->N);
  
  for(i=0;i<N;i++){
    if(r[i] > RMAX1){
      result[i] = 0.0;
    }else{
      result[i] = gsl_spline_eval (spline,log(r[i]), acc)*pow(2.0*pi*r[i],-1.5);
    }
  }
  
  free(r_FFT);
  free(ar);
  free(logr_FFT);
  FFTLog_free(fc);
  
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);

  return result;
}

double P1hss(double k, void *params)
{
  double res;
  ((cosmo_hm_params *)params)->k = k;
  double eps       = ((cosmo_hm_params *)params)->eps;
  error **err      = ((cosmo_hm_params *)params)->err;
  
  res = pow(k,1.5)*int_gsl(int_for_P1hss, params, logMmin, logMmax, eps, err);
  forwardError(*err, __LINE__, 0.0);

  
  return res;
}

double int_for_P1hss(double logM, void *params, error **err)
{
  double res, Omega_m, rhocrit;
  cosmo_hm *model  = ((cosmo_hm_params *)params)->model; 
  double a         = ((cosmo_hm_params *)params)->a; 
  double k         = ((cosmo_hm_params *)params)->k;
  double ng        = ((cosmo_hm_params *)params)->ng;
  int type         = ((cosmo_hm_params *)params)->type;
  double M         = exp(logM);

  double Ngal_sat = Ngal_s(model, M, model->Mstellar_min, model->Mstellar_max);
  
  /* Fourrier transform of halo profile */
  double c_sat = concentration_sat(model, M, a, err);
  forwardError(*err, __LINE__, 0.0);
  double rhohat =  rhohat_halo(model, k, M, a, c_sat, err);
  forwardError(*err, __LINE__, 0.0);

  /* Hala mass function parameters */
  ((cosmo_hm_params *)params)->asymptotic = 0;
  
  res  = dn_dlnM(M, params, err);
  forwardError(*err, __LINE__, 0.0);

  switch(type){
  case GG:
    res *= dsqr(Ngal_sat/ng)*rhohat*rhohat;
    break;
  case GM:

    Omega_m = Omega_m_halo(model, a, err);  forwardError(*err, __LINE__, 0.0);
    rhocrit = rho_crit_halo(model, a, err); forwardError(*err, __LINE__, 0.0);
    
    res *= Ngal_sat/ng * rhohat*rhohat*M/(rhocrit*Omega_m);
    forwardError(*err, __LINE__, 0.0);
    break;
  }
  
  return res;
}

double concentration_sat(cosmo_hm *model, double Mh, double a, error **err)
{
  double c, fac;
   /* to implement someday */

  //fac = 3.0*pow(model->Mstellar_min/1.0e10,-0.3);
  fac = 1.0;
  
  c = fac*concentration(model, Mh, a, err);
  forwardError(*err, __LINE__, 0.0);
  return c;
   
}


#define EPS 1.0e-5
double *xi_2h(cosmo_hm *model, double a, const double *r, int N, int type, error **err)
{
  int i,j;
  double *result = malloc_err(N*sizeof(double), err);  forwardError(*err, __LINE__, NULL);
  
  double xi_dm, amean;
  cosmo_hm_params params;
  params.cosmo  = model->cosmo;
  params.model  = model;
  params.a      = a;
  params.ng     = ngal_den(model, a, logMmax, model->Mstellar_min, model->Mstellar_max, err);
  forwardError(*err, __LINE__, NULL);
  params.err    = err;
  params.eps    = 1.0e-3;
  params.type   = type;

  FFTLog_config *fc;
  gsl_function Pk;

  /* New: consistency check */
  amean = 1.0 / (1.0 + zmean(model->redshift, 0, err));
  testErrorRetVA(fabs(amean - a) > EPS, hm_zmean_2h, "Precalculated xi_dm at a=%g not consistent with a=%g",
		 *err, __LINE__, NULL, amean, a);

  /* Tabulate xi_dm */
  if (model->xi_dm == NULL) {
    
    /* NEW: xi(r) now computed in a consistent way everywhere */
    int N_inter = 128;
    
    double *r_inter    = malloc_err(N_inter*sizeof(double), err); forwardError(*err, __LINE__, NULL);
    double dlogr_inter = (log(RMAX2) - log(RMIN1))/(double)N_inter;
    
    for(j=0;j<N_inter;j++) r_inter[j] = exp(log(RMIN1) + dlogr_inter*(double)j);
    
    /* initialize table */
    model->xi_dm  = init_interTable(N_inter, log(r_inter[0]), log(r_inter[N_inter-1]), dlogr_inter, 0.0, 0.0, err);
    forwardError(*err, __LINE__, NULL);
    
    /* compute xi(r) */
    double *xi = xi_P_NL(model, a, r_inter, N_inter, err);
    forwardError(*err, __LINE__, NULL);
    
    /* fill in the table */
    for(j=0;j<N_inter;j++) model->xi_dm->table[j] = xi[j];
    
    /* free memory */
    free(r_inter);
    free(xi);
    
  
  }

  /* xi 2h */
  Pk.function = &P2h;
  Pk.params   = &params;
  fc          = FFTLog_init(64, k_min, k_max_HOD, 0.0, 0.5);
  
  for(i=0;i<N;i++){
    if(r[i] < RMIN2){
      result[i] = 0.0;
    }else{
      if (model->hod == berwein02_hexcl || model->hod == leauthaud11){
	params.logMlim    = logM_lim(model,a,r[i],err);
	forwardError(*err, __LINE__, NULL);
	params.ngp        = ngal_den(model, a, params.logMlim, model->Mstellar_min, model->Mstellar_max, err);
	forwardError(*err, __LINE__, NULL);
	//debugging
	//params.logMlim   = logMmax;
	//params.ngp       = params.ng;

	params.bias_func  = &bias_r;
	/* Debugging */
	/* xi_dm             = xi_dm_NL(model,a,r[i],err); */
	xi_dm             = interpol_wr(model->xi_dm, log(r[i]), err);
	forwardError(*err, __LINE__, NULL);
	params.bias_fac   = pow(1.0+1.17*xi_dm,1.49)/pow(1.0+0.69*xi_dm,2.09);
      }else if (model->hod == berwein02){
	params.logMlim   = logMmax;
	params.ngp       = params.ng;
	params.bias_func = &bias_r;
	params.bias_fac  = 1.0;    /* No scale-dependent bias */
      }      
      
      if(params.ng < 1.0e-14 || params.ngp < 1.0e-14){
	/* NEW: was 1e-8 before and gave 0 2h-term for very massive halos */
	result[i] = 0.0;
      }else{
	result[i] = dsqr(params.ngp/params.ng)*xi_from_Pkr(&Pk, r[i], fc, err);
	forwardError(*err, __LINE__, NULL);
      }
      
    }
  }
  FFTLog_free(fc);
  
  return result;
}
#undef EPS

/* Angular correlation function of dark matter */
double *xi_P_NL(cosmo_hm *model, double a, const double *r, int N, error **err)
{
  int i;
  double *result = malloc_err(N*sizeof(double), err);  forwardError(*err, __LINE__, NULL);
  
  cosmo_hm_params params;
  params.cosmo  = model->cosmo;
  params.model  = model;
  params.a      = a;
  params.ng     = -1.0;
  params.err    = err;
  params.eps    = 1.0e-3;
  
  FFTLog_config *fc;
  gsl_function Pk;

  /* 3-d correlation function xi */
  Pk.function = &FFTLog_P_NL;
  Pk.params   = &params;
  fc          = FFTLog_init(64, k_min, k_max_HOD, 0.0, 0.5);
  
  for(i=0;i<N;i++){
    result[i] = xi_from_Pkr(&Pk, r[i], fc, err);
    forwardError(*err, __LINE__, NULL);
  }
  
  FFTLog_free(fc);
  return result;
}

double FFTLog_P_NL(double k, void *params)
{
  double res;
  cosmo_hm *model  = ((cosmo_hm_params *)params)->model; 
  double a         = ((cosmo_hm_params *)params)->a;
  error **err      = ((cosmo_hm_params *)params)->err;
  
  ((cosmo_hm_params *)params)->k = k;
  res = pow(k,1.5)*P_NL(model->cosmo, a, k, err);
  
  forwardError(*err, __LINE__, 0.0);
  return res;
}

double P2h(double k, void *params)
{
  double res;  
  cosmo_hm *model  = ((cosmo_hm_params *)params)->model; 
  double a         = ((cosmo_hm_params *)params)->a;
  double ng        = ((cosmo_hm_params *)params)->ngp;
  double eps       = ((cosmo_hm_params *)params)->eps;
  double logMlim   = ((cosmo_hm_params *)params)->logMlim;
  error **err      = ((cosmo_hm_params *)params)->err;
  int type         = ((cosmo_hm_params *)params)->type;

  ((cosmo_hm_params *)params)->k = k;
  
  res  = pow(k,1.5)*P_NL(model->cosmo, a, k, err);
  forwardError(*err, __LINE__, 0.0);
 
  switch(type){
  case GG:
    res *= dsqr(int_gsl(int_for_P2h, params, logMmin, logMlim, eps, err)/ng);
    forwardError(*err, __LINE__, 0.0);
    break;
  case GM:
    res *= int_gsl(int_for_P2h, params, logMmin, logMlim, eps, err)*int_gsl(int_for_P2h_dm, params, logMmin, logMlim, eps, err)/ng;
    forwardError(*err, __LINE__, 0.0);
    break;
  }
  
  return res;
}

double int_for_P2h(double logM, void *params, error **err)
{
  double res;
  double M         = exp(logM);
  cosmo_hm *model  = ((cosmo_hm_params *)params)->model; 
  double a         = ((cosmo_hm_params *)params)->a; 
  double k         = ((cosmo_hm_params *)params)->k;
  double (*bias_func)(double, void *) =  ((cosmo_hm_params *)params)->bias_func; 
  
  /* Fourrier transform of halo profile */
  double rhohat =  rhohat_halo(model, k, M, a, -1, err);
  forwardError(*err, __LINE__, 0.0);
  
  /* Hala mass function parameters */
  ((cosmo_hm_params *)params)->asymptotic = 0; 
  
  double b = bias_func(M, params);
  
  res = dn_dlnM(M, params,err)*b*Ngal(model, M, model->Mstellar_min, model->Mstellar_max)*rhohat;
  forwardError(*err, __LINE__, 0.0);
  
  return res;
}

double int_for_P2h_dm(double logM, void *params, error **err)
{
  double res, Omega_m, rhocrit;
  double M         = exp(logM);
  cosmo_hm *model  = ((cosmo_hm_params *)params)->model; 
  double a         = ((cosmo_hm_params *)params)->a; 
  double k         = ((cosmo_hm_params *)params)->k;
  double (*bias_func)(double, void *) =  ((cosmo_hm_params *)params)->bias_func; 
  
  /* Fourrier transform of halo profile */
  double rhohat =  rhohat_halo(model, k, M, a, -1, err);
  forwardError(*err, __LINE__, 0.0);
  
  /* Hala mass function parameters */
  ((cosmo_hm_params *)params)->asymptotic = 0; 
  
  double b = bias_func(M, params);
  
  Omega_m = Omega_m_halo(model, a, err);  forwardError(*err, __LINE__, 0.0);
  rhocrit = rho_crit_halo(model, a, err); forwardError(*err, __LINE__, 0.0);
  
  res = dn_dlnM(M, params,err)*b*(M/(rhocrit*Omega_m))*rhohat;
  forwardError(*err, __LINE__, 0.0);
  
  return res;
}


double bias_r(double M, void *params)
/* ---------------------------------------------------------------- *
 * Tinker 2005 eq B7. Depends on r via bias_fac.                    *
 * Used to be called "bias_tinker_r". Now call to bias_sc is        *
 * allowed depending on model->halo_bias.		            *
 * ---------------------------------------------------------------- */
{
  
  double b;
  cosmo_hm *model  = ((cosmo_hm_params *)params)->model; 
  double a         = ((cosmo_hm_params *)params)->a;
  double bias_fac  = ((cosmo_hm_params *)params)->bias_fac;
  error **err      = ((cosmo_hm_params *)params)->err;
  
  b = halo_bias(model, M, a, 1, err);
  forwardError(*err, __LINE__, 0.0);
  
  b *= sqrt(bias_fac);
  
  return b;
}

/* ---------------------------------------------------------------- *
 * HOD model functions                                              *
 * ---------------------------------------------------------------- */

double log10_fSHMR(double log10Mh, cosmo_hm *model){
  /* TO DO (?): interpolation + limits */
  
  // DEBUGGING to fix later
  if(model->beta < 0.1) return 1000000.0;

  int status;
  int iter = 0, max_iter = 100;
  double x, x_lo = logMmin - log10Mh, x_hi = logMmax - log10Mh;
  
  /* ATTENTION: do not choose eps larger than 
   * eps set for int_gsl in P1hss(..)
   */
  double eps    = 1.0e-5;
  
  gsl_function F;
  F.function = &log10_fSHMR_inv_minus_x;
  F.params   = model;
  model->x   = log10Mh;
  
  const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
  gsl_root_fsolver *s            = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);
  
  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      x      = gsl_root_fsolver_root (s);
      x_lo   = gsl_root_fsolver_x_lower (s);
      x_hi   = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi, 0, eps);
    }
  while (status == GSL_CONTINUE && iter < max_iter);
  
  gsl_root_fsolver_free (s);
  
  return x;
}

double log10_fSHMR_inv_minus_x(double log10Mstellar, void *p){
  cosmo_hm *model = (cosmo_hm *)p;
  return log10_fSHMR_inv(log10Mstellar, model) - model->x;
}


double log10_fSHMR_inv(double log10Mstellar, cosmo_hm *model){
  double result;  

  double A = pow(10.0, log10Mstellar)/model->Mstar0;
  result   = log10(model->M1) + model->beta*log10(A);
  result  += pow(A, model->delta)/(1.0 + pow(A, -model->gamma)) - 0.5;

  return result;
}

double Ngal_c(cosmo_hm *model, double Mh, double Mstellar_min, double Mstellar_max){
  /* Number of central galaxies per halo between 
   * Mstellar_min and Mstellar_max. Set Mstellar_max = -1 to get threshold number */
  double arg, result;
  
  // To fix later
  error *myerr = NULL, **err;
  err = &myerr;
  
  if(Mh < 1.0e8) return 0.0;   // MKDEBUG: Martin -> Jean: ??? -> [Jean changed from 10e10 to 1e8]
  
  /* /\* This assumes sigma_log_M is constant within a stellar mass bin *\/ */
  if (model->hod == leauthaud11){

    //DEBUGGING
    arg    = (log10(Mstellar_min) - log10_fSHMR(log10(Mh), model))/(sqrt(2.0)*sigma_log_M(model,log10(Mstellar_min), err));
    result = eta_cen(model, Mh, err)*0.5*(1.0-gsl_sf_erf(arg));
    if(Mstellar_max > 0){
      arg    = (log10(Mstellar_max) - log10_fSHMR(log10(Mh), model))/(sqrt(2.0)*sigma_log_M(model,log10(Mstellar_max), err));
      result -= eta_cen(model, Mh, err)*0.5*(1.0-gsl_sf_erf(arg));
    }
    
    
    /* if no approx on varying sigma_log... but very slow.*/
    
  /* cosmo_hm_params params; */
  /* params.model  = model; */
  /* params.Mh     = Mh; */

  /* result = int_gsl(int_for_Ngal_c, &params, log10(Mstellar_min), 13.00, 1.0e-3, err); */
  /* if(Mstellar_max > 0){ */
  /*   result -= int_gsl(int_for_Ngal_c, &params, log10(Mstellar_max), 13.00, 1.0e-3, err); */
  /* } */
    
  
    return result;
  }else{
    arg    = (log10(Mh/model->M_min)/model->sigma_log_M);
    result =  0.5*(1.0+gsl_sf_erf(arg));
    
    /* Ncen X central galaxy proportion */
    return model->eta*result;
  }
  
  /* eta must be <= 1 */
  testErrorExit(model->eta > 1.0, GM_ETA, "eta should be smaller than 1.0", *err, __LINE__);
}

double eta_cen(cosmo_hm *model, double Mh, error **err){
   /* fraction of centrals as function of halo mass:
    * = 1 if fcen1 = -1 and fcen2 = -2 (total sample)
    * increasing form if fcen1 > 0.0 (passive) 
    * decreasing form if fcen1 < 0.0 (star-forming)
    */

   double result;

   if(model->fcen1 < 0.0 && model->fcen2 < 0.0){
      result = 1.0;
   }else if(model->fcen1 > 0.0){
      result = 0.5*(1.0-gsl_sf_erf((model->fcen1-log10(Mh))/model->fcen2));
   }else{
      result = 0.5*(1.0-gsl_sf_erf((model->fcen1+log10(Mh))/model->fcen2));
  }
  
  forwardError(*err, __LINE__, 0.0);
  return result;
}

double int_for_Ngal_c(double log10Mstellar, void *params, error **err)
{
  double res;
  cosmo_hm *model  = ((cosmo_hm_params *)params)->model; 
  double Mh        = ((cosmo_hm_params *)params)->Mh;
  
  res = phi_c_Mstellar(model, log10Mstellar, Mh, err);
  forwardError(*err, __LINE__, 0.0);

  return res;
}

double phi_c_Mstellar(cosmo_hm *model, double log10Mstellar, double Mh, error **err){
  /* Returns P(M*|Mh). 
   * Leauthaud et al. 2011, Eq. 1. This is not exactly phi_c but phi_c * ln(10)
   * (to be integrated with variable change x -> ln x)
   */
  
  double A =  sigma_log_M(model, log10Mstellar,  err)*sqrt(2.0*pi);
  return exp(-dsqr(log10Mstellar - log10_fSHMR(log10(Mh), model))/(2.0*dsqr(sigma_log_M(model, log10Mstellar, err))))/A;
}

double av_Mh_given_Mstar(cosmo_hm *model, double Mstellar, double a, error **err){
  /* Returns <Mh|M*> */
  double norm, res;
  
  cosmo_hm_params params;
  params.model      = model;
  params.Mstellar   = Mstellar;
  params.a          = a;
  params.asymptotic = 0;
  
  res  = int_gsl(int_for_phi_c_Mh,      &params, logMmin, logMmax, 1.0e-3, err);
  norm = int_gsl(int_for_phi_c_Mh_norm, &params, logMmin, logMmax, 1.0e-3, err);
  
  return res/norm;

}

double int_for_phi_c_Mh(double logMh, void *params, error **err){
  
  cosmo_hm *model  = ((cosmo_hm_params *)params)->model; 
  double Mstellar  = ((cosmo_hm_params *)params)->Mstellar;
  double Mh       = exp(logMh);
  
  double res = phi_c_Mstellar(model, log10(Mstellar), Mh, err)*Mh*dn_dlnM(Mh, params, err);
  forwardError(*err, __LINE__, 0.0);

  return res;
  
}


double int_for_phi_c_Mh_norm(double logMh, void *params, error **err){
  
  cosmo_hm *model  = ((cosmo_hm_params *)params)->model; 
  double Mstellar  = ((cosmo_hm_params *)params)->Mstellar;
  double Mh       = exp(logMh);
  
  double res =  phi_c_Mstellar(model, log10(Mstellar), Mh, err)*dn_dlnM(Mh, params, err);

  forwardError(*err, __LINE__, 0.0);
  return res;
}


double sigma_log_M(cosmo_hm *model, double log10Mstellar, error **err){
  /* varying sigma_log_M */
  
  return model->sigma_log_M * pow(pow(10.0, log10Mstellar)/1.0e10, -model->beta_cut);
  
  /* this matches Leauthaud et al. (2012) */
  /* return model->sigma_log_M; */
  
}

double Ngal_s(cosmo_hm *model, double Mh, double Mstellar_min, double Mstellar_max){
  /*
    Number of satellite galaxies per halo mass.
  */
  
  double M0, log10Msat,log10Mcut,  result;
  
  if(Mh < 1.0e8) return 0.0;
  
  if (model->hod == leauthaud11){
    
  
    /* Stellar mass threshold */
    /* log10Msat = log10(model->B_sat) + model->beta_sat*log10_fSHMR_inv(log10(Mstellar_min), model) + 12.0*(1.0-model->beta_sat); */
    /* M0 = pow(10.0, log10_fSHMR_inv(log10(Mstellar_min), model)); */
    /* if(Mh - M0 > 0.0){ */
    /*   result = pow((Mh - M0)/pow(10.0, log10Msat), model->alpha); */
    /* }else{ */
    /*   result = 0.0; */
    /* } */
    /* /\* Stellar bin *\/ */
    /* if(Mstellar_max > 0){ */
    /*   log10Msat = log10(model->B_sat) + model->beta_sat*log10_fSHMR_inv(log10(Mstellar_max), model) + 12.0*(1.0-model->beta_sat); */
    /*   M0 = pow(10.0, log10_fSHMR_inv(log10(Mstellar_max), model)); */
    /*   if(Mh - M0 > 0.0){ */
    /* 	result -= pow((Mh - M0)/pow(10.0, log10Msat), model->alpha); */
    /*   } */
    /* } */


    // DEBUGGING
    log10Msat = log10(model->B_sat) + model->beta_sat*log10_fSHMR_inv(log10(Mstellar_min), model) + 12.0*(1.0-model->beta_sat);
    //log10Mcut = log10(model->B_cut) + 0.5*(log10_fSHMR_inv(log10(Mstellar_min), model)-12.0) + 12.0;
    log10Mcut = log10_fSHMR_inv(log10(Mstellar_min), model) - 0.5;
    M0        = pow(10.0, log10Mcut);
    if(Mh - M0 > 0.0){
      result = pow((Mh - M0)/pow(10.0, log10Msat), model->alpha);
    }else{
      result = 0.0;
    }
    /* Stellar bin */
    if(Mstellar_max > 0){
      log10Msat = log10(model->B_sat) + model->beta_sat*log10_fSHMR_inv(log10(Mstellar_max), model) + 12.0*(1.0-model->beta_sat);
      //log10Mcut = log10(model->B_cut) + 0.5*(log10_fSHMR_inv(log10(Mstellar_max), model)-12.0) + 12.0;
      log10Mcut = log10_fSHMR_inv(log10(Mstellar_max),  model) - 0.5;
      M0        = pow(10.0, log10Mcut);
      if(Mh - M0 > 0.0){
    	result -= pow((Mh - M0)/pow(10.0, log10Msat), model->alpha);
      }
    }
    
    
    /* /\* Leauthaud et al. (2011). Eqs. 13 & 14. *\/ */
    /* log10Msat = log10(model->B_sat) + model->beta_sat*(log10_fSHMR_inv(log10(Mstellar_min), model)-12.0) + 12.0; */
    /* log10Mcut = log10(model->B_cut) + model->beta_cut*(log10_fSHMR_inv(log10(Mstellar_min), model)-12.0) + 12.0; */
    /* result =  Ngal_c(model, Mh, Mstellar_min, -1) * pow(Mh/pow(10.0, log10Msat), model->alpha) * exp(-pow(10.0, log10Mcut)/Mh); */
    /* if(Mstellar_max > 0){ */
    /*   log10Msat = log10(model->B_sat) + model->beta_sat*(log10_fSHMR_inv(log10(Mstellar_max), model)-12.0) + 12.0; */
    /*   log10Mcut = log10(model->B_cut) + model->beta_cut*(log10_fSHMR_inv(log10(Mstellar_max), model)-12.0) + 12.0; */
    /*   result   -= Ngal_c(model, Mh, Mstellar_max, -1) * pow(Mh/pow(10.0, log10Msat), model->alpha) * exp(-pow(10.0, log10Mcut)/Mh); */
    /* } */
    
    /* tinker et al. (2013). Eqs. 13 & 14. */
    
    
    /* log10Msat = log10(model->B_sat) + model->beta_sat*log10_fSHMR_inv(log10(Mstellar_min), model) + 12.0*(1.0-model->beta_sat); */
    /* log10Mcut = log10(model->B_cut) + model->beta_cut*log10_fSHMR_inv(log10(Mstellar_min), model) + 12.0*(1.0-model->beta_cut); */
    /* result = pow(Mh/pow(10.0, log10Msat), model->alpha) * exp(-(pow(10.0, log10Mcut)+pow(10.0, log10_fSHMR_inv(log10(Mstellar_min), model)))/Mh); */
    /* if(Mstellar_max > 0){ */
    /*   log10Msat = log10(model->B_sat) + model->beta_sat*log10_fSHMR_inv(log10(Mstellar_max), model) + 12.0*(1.0-model->beta_sat); */
    /*   log10Mcut = log10(model->B_cut) + model->beta_cut*log10_fSHMR_inv(log10(Mstellar_max), model) + 12.0*(1.0-model->beta_cut); */
    /*   result   -= pow(Mh/pow(10.0, log10Msat), model->alpha) * exp(-(pow(10.0, log10Mcut)+pow(10.0, log10_fSHMR_inv(log10(Mstellar_max), model)))/Mh); */
    /* } */
	 
    
  }else{
    /* Two solutions here:
       1. N(M) = Ncen(1+Nsat) with Nsat \prop M^alpha and
       2. N(M) = Ncen + Nsat  with Nsat \prop Ncen * M^alpha 
       These are equivalent but <N(N-1)> must be changed accordingly. 
       It seems to matter for unusual sigma_log_M
       
       As in Brown et al. (2008): Ngal_c(model, M)*pow((M-M0)/model->M1, model->alpha);
       As in Blake et al. (2008): pow(M/model->M1,model->alpha);
       As in Zheng et al. (2005): pow((M - M0)/model->M1,model->alpha);    
    */     
    if (model->M0 < 0.0) {     /* Negative M0: set to M_min */
      M0 = model->M_min;
    } else {
      M0 = model->M0;
    }

    if (Mh - M0 < 0.0) return 0.0;
    
    /* [Jean] here we relax the constraint on ncen to allow a != 1 galaxy occupation function,
       but the equations need to be checked.
       It might change HOD parameters, but I wonder if this has any consequence on the 
       final conclusions... Anyway THIS HAS TO BE CHECKED as well as the consistency of the 
       equations for derived parameters. */
    
    result = pow((Mh - M0)/model->M1, model->alpha);
    
    // DEBUGGING
    //result = Ngal_c(model, Mh, Mstellar_min, -1) * pow((Mh - M0)/model->M1, model->alpha);
  }
  
  return MAX(0.0, result);
}

double Ngal(cosmo_hm *model, double Mh, double Mstellar_min, double Mstellar_max){
  /* Total number of galaxies per halo.
   */
  double result;
  
  result = Ngal_c(model, Mh, Mstellar_min, Mstellar_max) + Ngal_s(model, Mh, Mstellar_min, Mstellar_max);
  
  return result;
}

/*----------------------------------------------------------------*
 *Deduced quantities                                              *
 *----------------------------------------------------------------*/


double Mstar_tot_c(cosmo_hm *model, double Mh, double Mstellar_min, double Mstellar_max, error **err)
{
  /* returns the total stellar mass from central galaxies as function of halo mass and 
   * stellar bin. 
   * Currently only supported for leauthaud11's model. */

  if(model->hod == leauthaud11){
    double res, log_Mstellar_max;
    cosmo_hm_params params;
    params.model = model;
    params.Mh    = Mh;


    /* New: to be set in config file, -1 no longer allowed
    if(Mstellar_max < 0){
      log_Mstellar_max = 12.5;
    }else{
      log_Mstellar_max = log(Mstellar_max);
    }
    */
    
    log_Mstellar_max = log(Mstellar_max);

    res  = int_gsl(int_for_Mstar_tot_c, (void*)&params, log(Mstellar_min), log_Mstellar_max, 1.e-3, err);
    forwardError(*err, __LINE__, 0.0);
    res -= Ngal_c(model, Mh, Mstellar_max, -1)*Mstellar_max - Ngal_c(model, Mh, Mstellar_min, -1)*Mstellar_min;
    
    return res;
  }else{
    return -99.0;
  }
}

double int_for_Mstar_tot_c(double logMstar, void *params, error **err)
{
  
  double res, Mstar = exp(logMstar);
  cosmo_hm *model   = ((cosmo_hm_params *)params)->model;
  double Mh         = ((cosmo_hm_params *)params)->Mh; 
  
  res = Ngal_c(model, Mh, Mstar, -1)*Mstar;
  //forwardError(*err, __LINE__, 0.0);
  
  return res;
}

double Mstar_tot_s(cosmo_hm *model, double Mh, double Mstellar_min, double Mstellar_max, error **err)
{
  /* returns the total stellar mass from satellite galaxies as function of halo mass and 
   * stellar bin.
   * Currently only supported for leauthaud11's model. */
  
  if(model->hod == leauthaud11){
    double res, log_Mstellar_max;
    cosmo_hm_params params;
    params.model = model;
    params.Mh    = Mh;
    
    if(Mstellar_max < 0){
      log_Mstellar_max = 12.5;
    }else{
      log_Mstellar_max = log(Mstellar_max);
    }
    
    res  = int_gsl(int_for_Mstar_tot_s, (void*)&params, log(Mstellar_min), log_Mstellar_max, 1.e-3, err);
    forwardError(*err, __LINE__, 0.0);
    res -= Ngal_s(model, Mh, Mstellar_max, -1)*Mstellar_max - Ngal_s(model, Mh, Mstellar_min, -1)*Mstellar_min;
    
    return res;
  }else{

    return -99.0;
    
  }
}

double int_for_Mstar_tot_s(double logMstar, void *params, error **err)
{

  double res, Mstar = exp(logMstar);
  cosmo_hm *model   = ((cosmo_hm_params *)params)->model;
  double Mh         = ((cosmo_hm_params *)params)->Mh; 
  
  res = Ngal_s(model, Mh, Mstar, -1)*Mstar;
  //forwardError(*err, __LINE__, 0.0);
  
  return res;
}


double av_gal_bias(cosmo_hm *model, double a, error **err)
{
  /* returns the average galaxy bias wrt to redshift (depends on HOD). */
 
  double ng,res;
  cosmo_hm_params params;
  params.model = model;
  params.a    = a;
  params.asymptotic = 0;

  double min = logMmin;

  ng = ngal_den(model, a, logMmax, model->Mstellar_min, model->Mstellar_max, err);
  forwardError(*err, __LINE__, 0.0);
  res = int_gsl(int_for_av_gal_bias,(void*)&params,min,logMmax,1.e-3,err);
  forwardError(*err, __LINE__, 0.0);
  
  return res/ng;
}
double int_for_av_gal_bias(double logM, void *params, error **err)
{
  double b,res;
  double M = exp(logM);
  
  cosmo_hm *model  = ((cosmo_hm_params *)params)->model; 
  double        a  = ((cosmo_hm_params *)params)->a; 
  
  b   = halo_bias(model, M, a, 1, err); 
  res = dn_dlnM(M,params,err)*b*Ngal(model, M, model->Mstellar_min, model->Mstellar_max);
  forwardError(*err, __LINE__, 0.0);
  
  return res;
}

double av_halo_mass(cosmo_hm *model, double a, error **err)
{
  /* returns the mean halo mass (depends on HOD). */
 
  double ng,res;
  cosmo_hm_params params;
  params.model = model;
  params.a    = a;
  params.asymptotic = 0;

  double min = logMmin;

  ng = ngal_den(model, a, logMmax, model->Mstellar_min, model->Mstellar_max, err);
  forwardError(*err, __LINE__, 0.0);
  res = int_gsl(int_for_av_halo_mass,(void*)&params,min,logMmax,1.e-3,err);
  forwardError(*err, __LINE__, 0.0);

  testErrorRet(ng == 0, ce_infnan, "Division by zero (ng)", *err, __LINE__, 0.0);
  
  return res/ng;
}

double int_for_av_halo_mass(double logM, void *params, error **err)
{
  double res;
  double M         = exp(logM);
  cosmo_hm *model  = ((cosmo_hm_params *)params)->model; 
  
  res = M*dn_dlnM(M,params,err)*Ngal(model, M, model->Mstellar_min, model->Mstellar_max);
  forwardError(*err, __LINE__, 0.0);

  return res;
}

double av_frsat(cosmo_hm *model, double a, error **err)
{
  /* Returns the fraction of satellite galaxies wrt redshift. */
  
  double ng,res;
  cosmo_hm_params params;
  params.model = model;
  params.a    = a;
  params.asymptotic = 0;

  double min = logMmin;

  ng = ngal_den(model, a, logMmax, model->Mstellar_min, model->Mstellar_max, err);
  forwardError(*err, __LINE__, 0.0);
  res = int_gsl(int_for_av_frsat,(void*)&params,min,logMmax,1.e-3,err);
  forwardError(*err, __LINE__, 0.0);

  testErrorRet(ng == 0, ce_infnan, "Division by zero (ng)", *err, __LINE__, 0.0);
  
  return 1.0 - res/ng;
}
double int_for_av_frsat(double logM, void *params, error **err)
{
  double res;
  double M         = exp(logM);
  cosmo_hm *model  = ((cosmo_hm_params *)params)->model; 
  
  res = dn_dlnM(M,params,err)*Ngal_c(model, M, model->Mstellar_min, model->Mstellar_max);
  forwardError(*err, __LINE__, 0.0);

  return res;
}

double mass_weighted_av_stellar_mass(cosmo_hm *model, double a, error **err){
  /* returns the mean stellar mass weighted by the stellar mass 
     (for Delta_sigma point source signal).*/

  // DEBUGGING
  // return av_stellar_mass(model, a, err);


  double num, den;
  cosmo_hm_params params;
  params.model = model;
  params.a     = a;
  
  num = int_gsl(int_for_mass_weighted_av_stellar_mass, (void*)&params, log(model->Mstellar_min), log(model->Mstellar_max), 1.e-3, err);
  forwardError(*err, __LINE__, 0.0);
  
  den = int_gsl(int_for_mass_denum_weighted_av_stellar_mass, (void*)&params, log(model->Mstellar_min), log(model->Mstellar_max), 1.e-3, err);
  forwardError(*err, __LINE__, 0.0);
  
  testErrorRet(den < 1.0e-15, ce_infnan, "Division by zero (ng)", *err, __LINE__, 0.0);
  
  return num/den;
}

double int_for_mass_weighted_av_stellar_mass(double logM, void *params, error **err)
{
  double res;
  double M         = exp(logM);
  cosmo_hm *model  = ((cosmo_hm_params *)params)->model; 
  double a         = ((cosmo_hm_params *)params)->a; 
  
  double dlog = 0.001;
  
  res  = M*M*ngal_den(model, a, logMmax, exp(logM - dlog/2.0), exp(logM + dlog/2.0), err)/dlog;
  forwardError(*err, __LINE__, 0.0);
  
  return res;
}

double int_for_mass_denum_weighted_av_stellar_mass(double logM, void *params, error **err){
  
  double res;
  double M         = exp(logM);
  cosmo_hm *model  = ((cosmo_hm_params *)params)->model; 
  double a         = ((cosmo_hm_params *)params)->a; 
  
  double dlog = 0.001;
  
  res  = M*ngal_den(model, a, logMmax, exp(logM - dlog/2.0), exp(logM + dlog/2.0), err)/dlog;
  forwardError(*err, __LINE__, 0.0);
  
  return res;
}


double av_stellar_mass(cosmo_hm *model, double a, error **err)
{
  /* returns the mean stellar mass (depends on HOD).*/
 
  double ng, res;
  cosmo_hm_params params;
  params.model = model;
  params.a     = a;
  
  // DEBUGGING
  //double max = MAX(model->Mstellar_max, 1.0e13);
  
  ng = ngal_den(model, a, logMmax, model->Mstellar_min, model->Mstellar_max, err);
  forwardError(*err, __LINE__, 0.0);

  res = int_gsl(int_for_av_stellar_mass, (void*)&params, log(model->Mstellar_min), log(model->Mstellar_max), 1.e-3, err);
  forwardError(*err, __LINE__, 0.0);

  testErrorRet(ng == 0, ce_infnan, "Division by zero (ng)", *err, __LINE__, 0.0);
    
  return res/ng;

}
double int_for_av_stellar_mass(double logM, void *params, error **err)
{
  double res;
  double M         = exp(logM);
  cosmo_hm *model  = ((cosmo_hm_params *)params)->model; 
  double a        = ((cosmo_hm_params *)params)->a; 
  
  double dlog = 0.01;
  
  res  = M*ngal_den(model, a, logMmax, exp(logM - dlog/2.0), exp(logM + dlog/2.0), err)/dlog;
  forwardError(*err, __LINE__, 0.0);
  
  return res;
}

double av_halo_bias(cosmo_hm *model, double a, error **err)
{
  /*returns the average halo bias wrt to redshift (independent of HOD).
   THIS IS NOT TESTED*/
  double res1,res2;
  cosmo_hm_params params;
  params.model = model;
  params.a    = a;
  params.asymptotic = 0;

  //Jean: to check
  //double min = log10(model->M_min);
  double min = logMmin;

  res1 = int_gsl(int_for_av_halo_bias,(void*)&params,min,logMmax,1.e-3,err);
  forwardError(*err, __LINE__, 0.0);
  res2 = int_gsl(dn_dlnM_integrand,(void*)&params,min,logMmax,1.e-3,err);
  forwardError(*err, __LINE__, 0.0);

  testErrorRet(res2 == 0, ce_infnan, "Division by zero (res2)", *err, __LINE__, 0.0);

  return res1/res2;
}

double int_for_av_halo_bias(double logM, void *params, error **err)
{
  double b,res;
  double M = exp(logM);

  cosmo_hm *model  = ((cosmo_hm_params *)params)->model; 
  double        a  = ((cosmo_hm_params *)params)->a; 
  
  b   = halo_bias(model, M, a, 1, err);
  res = dn_dlnM(M,params,err)*b;
  forwardError(*err, __LINE__, 0.0);

  return res;
}

//Used by everyone
double dn_dlnM_integrand(double logM, void *params, error **err)
{
  double res;
  double M = exp(logM);

  res = dn_dlnM(M,params,err);
  forwardError(*err, __LINE__, 0.0);
  
  return res;
}

/* ---------------------------------------------------------------- *
 * Physical quantities                                              *
 * ---------------------------------------------------------------- */

double vc(cosmo_hm *model, double zmin, double zmax, error **err)
{
  /* Comoving volume per unit solid angle between zmin and zmax */
  
  double res;
  cosmo_hm_params params;
  params.model = model;
  
  if (zmin > -0.5 && zmax > -0.5) {
    res = int_gsl(int_for_vc, (void*)&params, zmin, zmax, 1e-4, err);
    forwardError(*err, __LINE__, 0.0);
  } else {
     Nz_hist *nz = get_nz(model->redshift, 0, err);
     forwardError(*err, __LINE__, 0.0);
     int i;
     double sum = 0.0;
     for(i=0;i<nz->nbins;i++){ //loop over z
        sum += nz->n[i]*dvdz(model->cosmo,1.0/(1.0+nz->z[i]),err)/nz->max;
        //sum += nz->n[i]*dvdz(model->cosmo,1.0/(1.0+nz->z[i]),err);
        forwardError(*err, __LINE__, 0.0);
     }
     res = sum*nz->dz;
  }
  
  return res;
}

double int_for_vc(double z, void *params, error **err)
{
  double res;
  cosmo_hm *model  = ((cosmo_hm_params *)params)->model; 

  res = dvdz(model->cosmo,1.0/(1.0+z),err);
  forwardError(*err, __LINE__, 0.0);

  return res;
}


double logM_lim(cosmo_hm *model, double a, double r, error **err)
{
  double ng_triax, ng, logM, dlogM, logM_hi, logM_lo;
  ng_triax = ngal_triax(model, a, r, err);
  forwardError(*err, __LINE__, 0.0);
  
  dlogM   = 0.01;
  logM_lo = logMmin;
  logM_hi = logMmax;
  
  do {
    logM = (logM_hi+logM_lo)/2.0;
    ng = ngal_den(model, a, logM, model->Mstellar_min, model->Mstellar_max, err);
    forwardError(*err, __LINE__, 0.0);
    if (ng<ng_triax) logM_lo = logM;
    else logM_hi = logM;
  } while (logM_hi-logM_lo > dlogM);
     logM = (logM_hi+logM_lo)/2.0;
  return logM;
}

double ngal_triax(cosmo_hm *model, double a, double r,  error **err){
  int i,j,N = 40;
  double x,y,P,sum1,sum2,logM1,logM2,R1,R2;
  
  cosmo_hm_params params;
  params.model      = model;
  params.a          = a;
  params.asymptotic = 0;
  params.Mstellar_min   = model->Mstellar_min;
  params.Mstellar_max   = model->Mstellar_max;
  
  double dlogM = (logMmax - logMmin)/(double)N;
  
  sum1 = 0.0;
  for(i=0;i<N;i++){
    logM1 = logMmin + dlogM*(double)i;
    R1 = r_vir(model, exp(logM1), a, err);
    forwardError(*err, __LINE__, 0.0);
    sum2 = 0.0;
    for(j=0;j<N;j++){
      logM2 = logMmin + dlogM*(double)j;
      R2 = r_vir(model,exp(logM2), a, err);
      forwardError(*err, __LINE__, 0.0);
      x = r/(R1 + R2);
      // This matches Leauthaud et al. (2011): 
      //x = r/(R1);
      y = (x - 0.8)/0.29;
      if (y<0) {
	sum2 += 0.0;
      } else if (y>1) {
	sum2 +=  int_for_ngal_den(logM2,(void*)&params,err);
	forwardError(*err, __LINE__, 0.0);
      } else {
	P  = (3.0 - 2.0*y)*y*y;
	sum2 +=  int_for_ngal_den(logM2,(void*)&params,err)*P;
	forwardError(*err, __LINE__, 0.0);
      }
    }
    sum1 +=  int_for_ngal_den(logM1,(void*)&params,err)*sum2*dlogM;
    forwardError(*err, __LINE__, 0.0);
  }
  
  return sqrt(sum1*dlogM);
}


double ngal_den(cosmo_hm *model, double a, double logMlim, double Mstellar_min, double Mstellar_max, error **err)
/* ---------------------------------------------------------------- *
 * galaxy number density per unit volume for a given HOD. In h^3 Mpc^-3.
 * logMlim is the maximum halo mass to integrate. Mstellar_min/max 
 * are the stellar mass bin limit (set Mstellar_max = -1 for 
 * threshold sample)
 * ---------------------------------------------------------------- */
{


  
  double res;
  
  cosmo_hm_params params;
  params.model         = model;
  params.a             = a;
  params.asymptotic    = 0;
  params.err           = err;
  params.Mstellar_min  = Mstellar_min;
  params.Mstellar_max  = Mstellar_max;
  
  res = int_gsl(int_for_ngal_den, (void*)&params, logMmin, logMlim, 1.e-5, err);
  forwardError(*err, __LINE__, 0.0);
  
  return res;
}

double int_for_ngal_den(double logM, void *params,  error **err) {
  cosmo_hm *model;
  double res, M, Mstellar_min, Mstellar_max;
  
  model        = ((cosmo_hm_params *)params)->model;
  Mstellar_min = ((cosmo_hm_params *)params)->Mstellar_min;
  Mstellar_max = ((cosmo_hm_params *)params)->Mstellar_max;

  M     = exp(logM);
  
  res   = Ngal(model, M, Mstellar_min, Mstellar_max)*dn_dlnM(M, params, err);
  forwardError(*err, __LINE__, -1.0);
  
  return res;
}


double ngal_den_c(cosmo_hm *model, double a, double logMlim, double Mstellar_min, double Mstellar_max, error **err)
/* ---------------------------------------------------------------- *
 * central galaxy number density per unit volume for a given HOD. In h^3 Mpc^-3.
 * logMlim is the maximum halo mass to integrate. Mstellar_min/max 
 * are the stellar mass bin limit (set Mstellar_max = -1 for 
 * threshold sample)
 * ---------------------------------------------------------------- */
{
  
  double res;
  
  cosmo_hm_params params;
  params.model         = model;
  params.a             = a;
  params.asymptotic    = 0;
  params.err           = err;
  params.Mstellar_min  = Mstellar_min;
  params.Mstellar_max  = Mstellar_max;
  
  res = int_gsl(int_for_ngal_den_c, (void*)&params, logMmin, logMlim, 1.e-5, err);
  forwardError(*err, __LINE__, 0.0);
  
  return res;
}

double int_for_ngal_den_c(double logM, void *params,  error **err) {
  cosmo_hm *model;
  double res, M, Mstellar_min, Mstellar_max;
  
  model        = ((cosmo_hm_params *)params)->model;
  Mstellar_min = ((cosmo_hm_params *)params)->Mstellar_min;
  Mstellar_max = ((cosmo_hm_params *)params)->Mstellar_max;

  M     = exp(logM);
  
  res   = Ngal_c(model, M, Mstellar_min, Mstellar_max)*dn_dlnM(M, params, err);
  forwardError(*err, __LINE__, -1.0);
  
  return res;
}


double ngal_den_s(cosmo_hm *model, double a, double logMlim, double Mstellar_min, double Mstellar_max, error **err)
/* ---------------------------------------------------------------- *
 * central galaxy number density per unit volume for a given HOD. In h^3 Mpc^-3.
 * logMlim is the maximum halo mass to integrate. Mstellar_min/max 
 * are the stellar mass bin limit (set Mstellar_max = -1 for 
 * threshold sample)
 * ---------------------------------------------------------------- */
{
  
  double res;
  
  cosmo_hm_params params;
  params.model         = model;
  params.a             = a;
  params.asymptotic    = 0;
  params.err           = err;
  params.Mstellar_min  = Mstellar_min;
  params.Mstellar_max  = Mstellar_max;
  
  res = int_gsl(int_for_ngal_den_s, (void*)&params, logMmin, logMlim, 1.e-5, err);
  forwardError(*err, __LINE__, 0.0);
  
  return res;
}

double int_for_ngal_den_s(double logM, void *params,  error **err) {
  cosmo_hm *model;
  double res, M, Mstellar_min, Mstellar_max;
  
  model        = ((cosmo_hm_params *)params)->model;
  Mstellar_min = ((cosmo_hm_params *)params)->Mstellar_min;
  Mstellar_max = ((cosmo_hm_params *)params)->Mstellar_max;

  M     = exp(logM);
  
  res   = Ngal_s(model, M, Mstellar_min, Mstellar_max)*dn_dlnM(M, params, err);
  forwardError(*err, __LINE__, -1.0);
  
  return res;
}

double ngal_den_vol(cosmo_hm *model, error **err)
{
  /* 
     Compute the n(z) weighted number density of galaxies from the halo density 
     E.g. Hamana et al. (2004) eq. (11)                                         
     ATTENTION: NOT TESTED
  */
  
  double ans1, ans2;
  
  ans1 = ngal_weighted(model, err); forwardError(*err, __LINE__, 0.0);
  ans2 = vc(model, -1.0, -1.0, err);  forwardError(*err, __LINE__, 0.0);
  
  return ans1/ans2;
}

double ngal_weighted(cosmo_hm *model, error **err){
  /* Ngal times comoving volume per unit solid angle */
  
  Nz_hist *nz  = get_nz(model->redshift, 0, err);
  forwardError(*err, __LINE__, 0.0);
  int i;
  double sum = 0.0, tmp;
  
  for(i=0;i<nz->nbins;i++){ //loop over z
    tmp  = ngal_den(model, 1.0/(1.0+nz->z[i]), logMmax, model->Mstellar_min, model->Mstellar_max, err)*nz->n[i];
    forwardError(*err, __LINE__, 0.0);
    
    tmp *= dvdz(model->cosmo,1.0/(1.0+nz->z[i]),err)/nz->max;
    forwardError(*err, __LINE__, 0.0);
    sum += tmp;
    
    /* NOT TESTED
       sum += nz->n[i]*dvdz(model->cosmo,1.0/(1.0+nz->z[i]),err);
       forwardError(*err, __LINE__, 0.0);
    */
  }
  
  return sum*nz->dz;
}

/*----------------------------------------------------------------*
 *FFTLog                                                          *
 *----------------------------------------------------------------*/


double xi_from_Pkr(gsl_function *ak, double r_prime, FFTLog_config *fc, error **err){
  
  /*In order to get xi(r) from FFTLog we need to re-write
    the 2-D Hankel transform into a 3-D one. 
    We know that
    
    sqrt(2/pi) sin(x) = sqrt(x) J_(1/2)(x)
    
    so 

         infinity
           /
   xi(r) = | P(k,r) sin(kr)/(kr) k^2/(2 PI^2) dr
           /             
          0

    becomes

                               infinity
                                 /
   xi(r) r^(3/2) (2 PI)^(3/2) =  | P(k,r) k^(3/2) J_(1/2) r dr
                                 /             
                                 0
   HOWEVER, the fast FFT is no longer correct if P(k) = P(k,r), that's
   why the algorithm above is only valid at the point r = r' and requires
   an interpolation. It means that you have to run this routine as many times
   as you need to evaluate xi(r'). The number of point fc->N below is the 
   parameter to play with in order to get accurate and fast evaluation of
   xi(r').
  */

  int i;
  
  double *r    = malloc_err(fc->N*sizeof(double),err); forwardError(*err, __LINE__, 0.0);
  double *ar   = malloc_err(fc->N*sizeof(double),err); forwardError(*err, __LINE__, 0.0);
  double *logr = malloc_err(fc->N*sizeof(double),err); forwardError(*err, __LINE__, 0.0);
  

 
  FFTLog(fc,ak,r,ar,-1,err); forwardError(*err, __LINE__, 0.0);
  
  for(i=0;i<fc->N;i++) logr[i] = log(r[i]);

  /* interpolation */
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline    = gsl_spline_alloc (gsl_interp_cspline,fc->N);
  
  gsl_spline_init (spline, logr, ar, fc->N);
  double result = gsl_spline_eval (spline,log(r_prime), acc);
  
  free(r);
  free(ar);
  free(logr);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  
  return result*pow(2.0*pi*r_prime,-1.5);
}


void FFTLog(FFTLog_config *fc, const gsl_function *ar_in, double *k, double *ak_out, int dir, error **err){
  /* Hamilton 2000. http://casa.colorado.edu/~ajsh/FFTLog/

     The FFTLog algorithm for taking the discrete Hankel transform, equation (22), 
     of a sequence an of N logarithmically spaced points is:
     
     * FFT an to obtain the Fourier coefficients cm, equation (13);
     * multiply by um given by equations (18) and (19) to obtain cm um;
     * FFT cm um back to obtain the discrete Hankel transform ãn, equation (21). 

     A variant of the algorithm is to sandwich the above operations with power 
     law biasing and unbiasing operations. For example, one way to take the 
     unbiased continuous Hankel transform Ã(k) of a function A(r), equation (4), 
     is to bias A(r) and Ã(k) with power laws, equation (3), and take a biased Hankel transform, 
     equation (2). The discrete equivalent of this is:

     * Bias An with a power law to obtain an = An rn-q, equation (3);
     * FFT an to obtain the Fourier coefficients cm, equation (13);
     * multiply by um given by equations (18) and (19) to obtain cm um;
     * FFT cm um back to obtain the discrete Hankel transform ãn, equation (21);
     * Unbias ãn with a power law to obtain Ãn = ãnkn-q, equation (3). 
*/

  
  
  int i;
  
  double logrmin  = log(fc->min);
  double logrmax  = log(fc->max);
  double r, dlogr = (logrmax - logrmin)/(double)fc->N;
  double logrc    = (logrmax + logrmin)/2.0;
  double nc       = (double)(fc->N+1)/2.0-1;
  double logkc    = log(fc->kr)-logrc;
  
  /* write signal */
  for(i=0; i<fc->N; i++){
    k[i] = exp(logkc+((double)i-nc)*dlogr);
    r  = exp(logrc+((double)i-nc)*dlogr);
    fc->an[i][0] = ar_in->function(r,(void*)(ar_in->params))*pow(r,-(double)dir*fc->q);
    fc->an[i][1] = 0.0;
  }
  
 
  /* cm's: FFT forward */
  fftw_execute(fc->p_forward);

  /* um*cm */
  fc->cmum[0][0] = fc->cm[0][0]*fc->um[0][0] - fc->cm[0][1]*fc->um[0][1];
  fc->cmum[0][1] = fc->cm[0][0]*fc->um[0][1] + fc->cm[0][1]*fc->um[0][0];
  for(i=1;i<fc->N/2+1;i++){
    fc->cmum[i][0] = fc->cm[i][0]*fc->um[i][0] - fc->cm[i][1]*fc->um[i][1];
    fc->cmum[i][1] = fc->cm[i][0]*fc->um[i][1] + fc->cm[i][1]*fc->um[i][0];
    /* Hermitian symetry (i.e. to get a real signal after FFT back) */
    fc->cmum[fc->N-i][0] = fc->cmum[i][0];
    fc->cmum[fc->N-i][1] = -fc->cmum[i][1];
  }
  
  /* ak's: FFT backward */
   fftw_execute(fc->p_backward);
   
   /* reverse array... */
  for(i=0;i<fc->N/2;i++) FFTLog_SWAP(fc->ak[i][0],fc->ak[fc->N-i-1][0]);
  
  /* ...and write ak(k) */
  for(i=0;i<fc->N;i++) ak_out[i] = fc->ak[i][0]*pow(k[i],-(double)dir*fc->q)/(double)fc->N;

  


  
  return;
}

FFTLog_config *FFTLog_init(int N, double min, double max, double q, double mu){
  /* Initializes what FFTLog needs. */
  

  FFTLog_config *fc = (FFTLog_config*)malloc(sizeof(FFTLog_config));

  /* FFTW3 Initialization */
  fc->min        = min;
  fc->max        = max;
  fc->q          = q;
  fc->mu         = mu;
  fc->N          = N;
  fc->an         = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
  fc->ak         = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
  fc->cm         = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
  fc->um         = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
  fc->cmum       = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
  fc->p_forward  = fftw_plan_dft_1d(N,fc->an,fc->cm,FFTW_FORWARD,FFTW_ESTIMATE);
  fc->p_backward = fftw_plan_dft_1d(N,fc->cmum,fc->ak,FFTW_BACKWARD,FFTW_ESTIMATE);
  
  /* um's */
  FFTLog_complex z,result;
  double L = log(max)-log(min);
  fc->kr   = 1.0;
  int i;
  
  for(i=0;i<fc->N/2+1;i++){
    z.re   = fc->q;
    z.im   = 2.0*pi*(double)i/L;
    result = FFTLog_U_mu(mu,z);
    
    /* Multiply by (kr)^-2PIim/L */
    result.amp *= 1.0;
    result.arg += -2.0*pi*(double)i*log(fc->kr)/L;
    
    fc->um[i][0] = result.amp*cos(result.arg);
    fc->um[i][1] = result.amp*sin(result.arg);
  }
  
  /* If N even, mutiply by real part only */
  if(PARITY(fc->N) == EVEN) fc->um[fc->N/2][1] = 0.0;
  
  return fc;
}

void FFTLog_free(FFTLog_config *fc){
  
  fftw_destroy_plan(fc->p_forward);
  fftw_destroy_plan(fc->p_backward);
  fftw_free(fc->an);
  fftw_free(fc->ak);
  fftw_free(fc->cm);
  fftw_free(fc->um);
  fftw_free(fc->cmum);
  free(fc);

  return;
}

FFTLog_complex FFTLog_U_mu(double mu, FFTLog_complex z){
  /*Computes 2^z Gamma[(mu + 1 - z)/2]/Gamma[(mu + 1 - z)/2]
              1                2                 3
  */
  double amp1,arg1;
  gsl_sf_result lnamp2,arg2,lnamp3,arg3;
  
  FFTLog_complex result;
  
  /* 2^z */
  amp1 = exp(z.re*log(2.0));
  arg1 = z.im*log(2.0);
  
  /* Gamma 1 */
  FFTLog_complex zplus;
  zplus.re = (mu + 1.0 + z.re)/2.0;
  zplus.im = z.im/2.0;
  gsl_sf_lngamma_complex_e(zplus.re,zplus.im,&lnamp2,&arg2);
  
  /* Gamma 2 */
  FFTLog_complex zminus;
  zminus.re = (mu + 1.0 - z.re)/2.0;
  zminus.im = - z.im/2.0;
  gsl_sf_lngamma_complex_e(zminus.re,zminus.im,&lnamp3,&arg3);

  /* Result */
  result.amp = amp1*exp(lnamp2.val)*exp(-lnamp3.val);
  result.arg = arg1 + arg2.val - arg3.val;
  result.re = result.amp*cos(result.arg);
  result.im = result.amp*sin(result.arg);
  
  return result;
}

/*----------------------------------------------------------------*
 *Utils                                                           *
 *----------------------------------------------------------------*/

void updateFrom_hm(cosmo_hm* avant, cosmo_hm* apres, error **err)
{
   int change;

   if (change_massfct_params(avant, apres)) {
      apres->A     = -1.0;
   }
   if (change_Mstar(avant, apres)) {
      apres->Mstar = -1.0;
   }
   if (change_rhohat_halo(avant,apres)) {
      del_interTable2D(&(apres->rhohat));
   }
   if (change_Pthg(avant,apres)) {
     del_interTable(&(apres->xir));
     apres->a_xir   = -1.0;
   }
   if (change_Pth(avant,apres)) {
     del_interTable2D(&(apres->Pthdm));
   }
   if (change_massfct(avant, apres) || change_massfct_params(avant, apres)) {
     set_massfct(apres->massfct, &apres->nmz_a, &apres->nmz_p, err);
     forwardError(*err, __LINE__,);
   }
   if (change_sigma_R_sqr(avant, apres)) {
     del_splineTable(&(apres->sigRsqr));
   }
   change = change_zmean(avant->redshift, apres->redshift, err);
   forwardError(*err, __LINE__,);
   if (change) {
      del_interTable(&apres->xi_dm);
   }

   updateFrom(avant->cosmo, apres->cosmo, err);
   forwardError(*err, __LINE__,);

   updateFrom_redshift(avant->redshift, apres->redshift);
}

int change_vc(cosmo_hm *avant, cosmo_hm *apres)
{
   if (change_prob(avant->redshift, apres->redshift)) return 1;
   if (change_w(avant->cosmo, apres->cosmo)) return 1;

   if (NCOCLOSE(avant->cosmo,apres->cosmo,Omega_m) || NCOCLOSE(avant->cosmo,apres->cosmo,Omega_de) ||
       NCOCLOSE(avant->cosmo,apres->cosmo,w0_de)   || NCOCLOSE(avant->cosmo,apres->cosmo,w1_de) ||
       NCOEQ(avant->cosmo,apres->cosmo,de_param)   || NCOCLOSE(avant->cosmo,apres->cosmo,h_100) ||
       NCOCLOSE(avant->cosmo,apres->cosmo,Omega_b) || NCOCLOSE(avant->cosmo,apres->cosmo,Omega_nu_mass) ||
       NCOCLOSE(avant->cosmo,apres->cosmo,Neff_nu_mass))
     return 1;

   return 0;
}

int change_HOD(cosmo_hm *avant, cosmo_hm *apres)
{
   if (NCOCLOSE(avant, apres, M_min) || NCOCLOSE(avant, apres, M1) ||
       NCOCLOSE(avant, apres, M0) || NCOCLOSE(avant, apres, sigma_log_M) ||
       NCOCLOSE(avant, apres, alpha) || NCOEQ(avant, apres, hod)) return 1;
   return 0;
}

int change_Pthg(cosmo_hm* avant, cosmo_hm* apres)
{
   if (change_Pth(avant, apres)) return 1;
   if (change_HOD(avant, apres)) return 1;
   if (NCOEQ(avant, apres, hod)) return 1;

   return 0;
}

int change_w_of_theta(cosmo_hm *avant, cosmo_hm *apres)
{
   if (change_Pthg(avant,apres)) return 1;
   if (change_prob(avant->redshift,apres->redshift)) return 1;

   return 0;
}

/* ============================================================ *
 * Put the 'hist' format into 'Nz_hist' format.
 * The array par_nz for 'hist":
 * 0  1  | 2  3  ... n       | n+1 n+2 ... 2n
 * -----------------------------------------------
 * z0 zn | z1 z2 ... z_{n-1} | N0  N1  ... N_{n-1}
 *   nbins = n
 * ============================================================ */
Nz_hist *get_nz(redshift_t *nz_mk, int n_bin, error **err)
{
  testErrorRetVA(n_bin<0 || n_bin>=nz_mk->Nzbin, redshift_Nzbin,
		 "Requested redshift bin %d out of range [0; %d]",
		 *err, __LINE__, NULL, n_bin, nz_mk->Nzbin-1);

  testErrorRetVA(nz_mk->nofz[n_bin] != hist, redshift_type,
		 "Invalid redshift type nofz='%s'(%d), hase to be '(%s)'(%d)' for HOD",
		 *err, __LINE__, NULL, snofz_t(nz_mk->nofz[n_bin]), nz_mk->nofz[n_bin], snofz_t(hist), hist);

  int i;
  Nz_hist *result = malloc_err(sizeof(Nz_hist), err);                forwardError(*err, __LINE__, NULL);
 
  result->nbins   = (nz_mk->Nnz[n_bin]-1)/2;
  result->z       = malloc_err(result->nbins*sizeof(double), err);   forwardError(*err, __LINE__, NULL);
  result->n       = malloc_err(result->nbins*sizeof(double), err);   forwardError(*err, __LINE__, NULL);
  
  result->z[0]  = nz_mk->par_nz[n_bin*nz_mk->Nnz_max];
  result->n[0]  = nz_mk->par_nz[n_bin*nz_mk->Nnz_max + result->nbins+1];
  result->max   = result->n[0];

  for(i=1;i<result->nbins;i++){
    result->z[i] = nz_mk->par_nz[n_bin*nz_mk->Nnz_max + i+1];
    result->n[i] = nz_mk->par_nz[n_bin*nz_mk->Nnz_max + result->nbins+i+1];
    result->max = MAX(result->max,result->n[i]);
  }
  
  /* normalize AND return norm */
  result->dz = result->z[1]-result->z[0];
  result->Norm = normalize(result->n,result->nbins,result->dz,err);
  forwardError(*err, __LINE__, NULL);
  result->max /= result->Norm;

  return result;
}

void  free_nz( Nz_hist *nz){
  free(nz->z);
  free(nz->n);
  free(nz);
}

double normalize(double *data, size_t Nbins, double binsize,error **err){
  /*Normalize data with binsize AND return the norm*/
  size_t i;
  double Norm;
  Norm = 0.0;
  for(i=0;i<Nbins;i++) Norm += data[i]*binsize;
  for(i=0;i<Nbins;i++) data[i] /= Norm;

  return Norm;

}


wt_t *init_wt_t(double nlines, int wcov, error **err)
{
   wt_t *wtr;

   wtr        = malloc_err(sizeof(wt_t), err);                forwardError(*err, __LINE__, NULL);
   wtr->nbins = nlines;
   wtr->th    = malloc_err(sizeof(double)*nlines, err);       forwardError(*err, __LINE__, NULL);
   wtr->w     = malloc_err(sizeof(double)*nlines, err);       forwardError(*err, __LINE__, NULL);
   if (wcov==0) {
      wtr->werr = malloc_err(sizeof(double)*nlines, err);
      forwardError(*err, __LINE__, NULL);
      wtr->wcov = NULL;
   } else {
      wtr->werr = NULL;
      wtr->wcov = malloc_err(sizeof(double)*nlines*nlines, err);
      forwardError(*err, __LINE__, NULL);
   }
   //   wtr->RR = malloc_err(sizeof(double)*nlines, err);       forwardError(*err, __LINE__, NULL);

   return wtr;
}

void free_wt_t(wt_t **model)
{
   wt_t *s;

   s = *model;
   free(s->th);
   free(s->w);
   if (s->werr) free(s->werr);
   if (s->wcov) free(s->wcov);
   free(s->RR);
   free(s->th_RR);
   free(s);
   model = NULL;
}

wt_t *read_wtheta(const char *data_file_name, intconst_t intconst_type, double delta, 
		  double intconst, const char *ran_file_name, error **err)
{
  /* read a w(theta) file into an array and correct it from integral constraint */
  
  size_t i, Ncol, res;
  char line[NFIELD*NCHAR], item[NFIELD*NCHAR];
  
  FILE *fileIn = fopen_err(data_file_name, "r", err); forwardError(*err, __LINE__, NULL);
  
  /* Count the number of lines */
  size_t N;
  N = 0;
  rewind(fileIn);
  while(fgets(line,NFIELD*NCHAR,fileIn) != NULL)
    if(*(line) != '#') N++;
  
  wt_t *wtr;
  wtr = init_wt_t(N, 0, err);
  forwardError(*err, __LINE__, NULL);
  
  rewind(fileIn);
  i=0;
  while(fgets(line,NFIELD*NCHAR,fileIn) != NULL){
     res = getStrings(line,item," ",&Ncol,err);
     forwardError(*err, __LINE__, NULL);
     if(res) {
      wtr->th[i]   = getDoubleValue(item,1);
      wtr->w[i]    = getDoubleValue(item,2);
      wtr->werr[i] = getDoubleValue(item,3);
      if(intconst_type == constant){
	wtr->w[i]    = wtr->w[i]*pow(wtr->th[i],-delta)/
	  (pow(wtr->th[i],-delta)-intconst);
      }
      i++;
     }
  }
  
  fclose(fileIn);
  
  /* random points for model-dependent integral constraint calculation */
  if(intconst_type == random_file){ 
    
    fileIn = fopen_err(ran_file_name, "r", err); forwardError(*err, __LINE__, NULL);
    
    N = 0;
    rewind(fileIn);
    while(fgets(line,NFIELD*NCHAR,fileIn) != NULL)
      if(*(line) != '#') N++;
    
    wtr->nbins_RR = N; 
    wtr->RR       = malloc_err(sizeof(double)*N, err);       forwardError(*err, __LINE__, NULL);
    wtr->th_RR    = malloc_err(sizeof(double)*N, err);       forwardError(*err, __LINE__, NULL);

    rewind(fileIn);
    i=0;
    while(fgets(line,NFIELD*NCHAR,fileIn) != NULL){
      res = getStrings(line,item," ",&Ncol,err);
      forwardError(*err, __LINE__, NULL);
      if(res) {
	wtr->th_RR[i] = getDoubleValue(item,1);
	wtr->RR[i]    = getDoubleValue(item,2);
	i++;
      }
    }
    
    fclose(fileIn);
    
  }
    

  return wtr;
}

int getStrings(char *line, char *strings, char *delimit, size_t *N, error **err){
  /*Extract each word/number in line separated by delimit and returns 
    the array of items in strings.*/
  int i,j,begin,length;
  
  if(line == NULL || line[0] == '\n' || line[0] == '#' || line[0] == '\0') return 0;
  
  i = 0;
  j = 0;
  while(line[i] != '\0' && line[i] != '#' && line[i] != '\n'){
    begin = i;
    while(line[i] == *delimit || line[i] == '\t' && (line[i] != '\0' || line[i] != '#' || line[i] != '\n')) i++;
    begin = i;
    while(line[i] != *delimit && line[i] != '\t' && line[i] != '\0' && line[i] != '#' && line[i] != '\n') i++;
    length = i - begin;
    if(length > 0){
      strncpy(strings+NCHAR*j,&line[begin],length);
      strcpy(strings+NCHAR*j+length,"\0");
      j++;
    }
  }
  
  (*N) = j;
  
  if(*N > 0){
    return 1;
  }else{
    return 0;
  }
}


void read_config_hm_file(config_hm_info *config, char cname[], error **err)
{
   FILE *F;
   config_element c = {0, 0.0, ""};

   F = fopen_err(cname, "r", err);                   forwardError(*err, __LINE__,);

   CONFIG_READ_S(config, WTHETA, s, F, c, err);
   CONFIG_READ_S(config, COVNAME, s, F, c, err);
   CONFIG_READ_S(config, DNDZ, s, F, c, err);
   CONFIG_READ_S(config, OUTFILE, s, F, c, err);

   CONFIG_READ(config, ngal, d, F, c, err);
   CONFIG_READ(config, ngalerr, d, F, c, err);
   //   CONFIG_READ(config, area_deg, d, F, c, err);
   CONFIG_READ(config, intconst, d, F, c, err);
   CONFIG_READ(config, delta, d, F, c, err);
   CONFIG_READ(config, nbins, i, F, c, err);

   CONFIG_READ(config, alpha_min, d, F, c, err);
   CONFIG_READ(config, alpha_max, d, F, c, err);
   CONFIG_READ(config, dalpha, d, F, c, err);

   CONFIG_READ(config, M1_min, d, F, c, err);
   CONFIG_READ(config, M1_max, d, F, c, err);
   CONFIG_READ(config, dM1, d, F, c, err);

   CONFIG_READ(config, Mmin_min, d, F, c, err);
   CONFIG_READ(config, Mmin_max, d, F, c, err);
   CONFIG_READ(config, dMmin, d, F, c, err);

   fclose(F);
}


/*----------------------------------------------------------------*
 * Obsolete stuff                                                 *
 *----------------------------------------------------------------*/


nz_t *read_dndz_OBSOLETE(char *infile,error **err)
{
  FILE *inF; 
  char str[MAXCHAR];
  char vectok[] = {",() \t\n\r"};
  int i, nbin;
  double nzr1, nzr2;
  char *str2_new;
  double sum;
  double zm;
  nz_t *nzr;

  nzr      = malloc_err(sizeof(nz_t), err);          forwardError(*err, __LINE__, NULL);
  nzr->z   = malloc_err(sizeof(double)*NLINES, err); forwardError(*err, __LINE__, NULL);
  nzr->fac = malloc_err(sizeof(double)*NLINES, err); forwardError(*err, __LINE__, NULL);
  inF      = fopen_err(infile, "r", err);            forwardError(*err, __LINE__, NULL);
  
  i = 0;
  zm = sum = 0.0;
  while (fgets(str, MAXCHAR, inF)) {
    if (!(str2_new = strtok(str, vectok)))
      break;
    if (*str == (char)'#')
      continue;
    //theta = atof(str2_new); // not used
     nzr1        = (double)atof(str2_new);
     str2_new    = strtok(NULL, vectok);
     testErrorRet(!str2_new, hm_io, "Error reading the lines", *err, __LINE__, NULL);
     nzr2        = (double)atof(str2_new);
     nzr->z[i]   = fabs(nzr2-nzr1)/2.0 + nzr1;
     str2_new    = strtok(NULL, vectok);
     testErrorRet(!str2_new, hm_io, "Error reading the lines", *err, __LINE__, NULL);
     nzr->fac[i] = (double)atof(str2_new);
     sum        += (nzr2-nzr1)*nzr->fac[i];
     zm         += (nzr2-nzr1)*nzr->fac[i]*nzr->z[i];
     i++;
     testErrorRetVA(i>=NLINES, hm_overflow, "Too many lines (>=%d), increase NLINES in hod.h",
		    *err, __LINE__, NULL, NLINES);
  }
  
  fclose (inF);
  nbin = i;
  /* normalise */
  for (i=0;i<nbin;i++) {
    nzr->fac[i]=nzr->fac[i]/sum; 
    /*       printf ("%f %f\n", nzr->z[i], nzr->fac[i]); */
  }
  nzr->nbins=nbin;
  zm/=sum;
  nzr->zm=zm;
  printf ("zm=%f\n", zm);
  /* end */
  
  /* also initialize spline function here */ 
  nzr->ypn = sm2_vector(1,nzr->nbins,err); 
  forwardError(*err, __LINE__, NULL);
  
  sm2_spline(nzr->z-1,nzr->fac-1,nzr->nbins,1.0e30,1.0e30,nzr->ypn,err);
  forwardError(*err, __LINE__, NULL);
  
  return nzr;
}


double Ngal_mean(cosmo_hm *model, double a, error **err)
{
  /*returns the mean number of galaxies per halo. NOT TESTED.*/
  double res;
  cosmo_hm_params params;
  params.model = model;
  params.a    = a;
  params.asymptotic = 0;
  
  double ng = ngal_den(model, a, logMmax, model->Mstellar_min, model->Mstellar_max, err);
  forwardError(*err, __LINE__, 0.0);
  
  res = ng/int_gsl(dn_dlnM_integrand,(void*)&params,log(model->M_min),logMmax,1.e-3,err);
  forwardError(*err, __LINE__, 0.0);
  return res;
}



/*----------------------------------------------------------------*
 *Add by Melody Wolk 18/07/2011                                   *
 *----------------------------------------------------------------*/

/*----------------------------------------------------------------*
 *Compute total xiofr                                             *
 *----------------------------------------------------------------*/

#define EPS 1.0e-5
double xir_mwolk(cosmo_hm *model, double a, double r, error **err)
{
   double dr, logr, res;
   int j, nbins=60;
   
   testErrorRetVA(model->hod!=berwein02 && model->hod!=berwein02_hexcl, hm_hodtype,
		  "Invalid hod type (%d), only berwein02(%d) and berwein02_hexcl(%d) hod implemented for xi(r)",
		  *err, __LINE__, 0.0, model->hod, berwein02, berwein02_hexcl);
   
   /*
     testErrorRetVA(pofk_xir!=p1hgcs && pofk_xir!=p1hgss && pofk_xir!=p2hg && pofk_xir!=p1hgcs+p1hgss &&
     pofk_xir!=p1hgcs+p2hg && pofk_xir!=p1hgss+p2hg && pofk_xir!=p1hgcs+p1hgss+p2hg, hm_pofk,
     "Wrong pofk type (%d), has to be p1gcs(%d), p1hgss(%d), p2hg(%d) or a bit-wise combination thereof",
     *err, __LINE__, 0.0, pofk_xir, p1hgcs, p1hgss, p2hg);
     */

   /* Delete previous table if calculated for different scale factor a */
   // MKDEBUG TODO:
   // (1) maybe get rid of table xir
   // (2) a_xir only set at one place
   if (fabs(model->a_xir-a)>EPS) {
     del_interTable(&model->xir);
     //  printf("Check2 %lf %lf\n", model->a_xir, a);
   }

   double *ri =malloc(nbins*sizeof(double));
   double *xi =malloc(nbins*sizeof(double));
   
   if (model->xir==NULL) {
     dr       = (log10(RMAX)-log10(RMIN))/((double)nbins);
     model->xir   = init_interTable(nbins, log10(RMIN), log10(RMAX), dr, 0.0, 0.0, err);
     forwardError(*err, __LINE__, 0.0);
     model->a_xir = a;
     
     
     /* Earlier: Tabulation in a as well */
     
     for (j=0; j<nbins; j++) {
       ri[j] = pow(10., log10(RMIN)+(double)j*dr);
     }
     res=0.;
     double *xi1hgcs = xiofr(model, p1hgcs, ri, nbins, GG, err);     forwardError(*err, __LINE__, -1.0);
     double *xi1hgss = xiofr(model, p1hgss, ri, nbins, GG, err);     forwardError(*err, __LINE__, -1.0);
     double *xi2hg   = xiofr(model, p2hg, ri, nbins, GG, err);         forwardError(*err, __LINE__, -1.0);


     for (j=0; j<nbins; j++) {
       xi[j]= xi1hgcs[j]+xi1hgss[j]+xi2hg[j];
     }
     
     for (j=0; j<nbins; j++) {
       model->xir->table[j]=xi[j];
       // printf("xir: %lf %lf %d\n",model->xir->table[j], pow(10., log10(RMIN)+(double)j*dr), j);
     }


   }

   if (r<RMIN || r>RMAX) return 0.0;

   logr=log10(r) ;
  
   res = interpol_wr(model->xir, logr, err);
   forwardError(*err, __LINE__, 0.0);
   //printf("test: %lf %lf\n", model->xir->upper, model->xir->lower) ;
   // printf("r:%lf xir:%lf\n", r, res) ;
   return res;
}
#undef EPS


/*----------------------------------------------------------------*
 *Read mwolk format files                                         *
 *----------------------------------------------------------------*/

wt_t *read_wtheta_mwolk_OBSOLETE(const char *name, double delta, double intconst, error **err)
{
  /* read in a melodywolk-format wp-rp file into an array. */
  FILE *inF; 
  char str[MAXCHAR];
  double rp,wp, errwp;
  wt_t *wtr; 
  int i;

  wtr = init_wt_t(NLINES, 0, err);
  forwardError(*err, __LINE__, NULL);

  i = 0;
  inF = fopen_err(name, "r", err);                        forwardError(*err, __LINE__, NULL);

  while (fgets(str, MAXCHAR, inF)!=NULL) {
    sscanf(str,"%lf %lf %lf\n", &rp, &wp, &errwp) ;
    wtr->th[i]=rp ;
  
    wtr->werr[i] = errwp;
    
    /* Apply the correction for the integral constraint*/
    wtr->w[i]=wp ; /* I need to think about that */
    printf("coord: %lf %lf\n", wtr->th[i], wtr->w[i]) ;
    
     i++;
     
     testErrorRetVA(i>=NLINES, hm_overflow, "Too many lines (>=%d), increase NLINES in hod.h",
		    *err, __LINE__, NULL, NLINES);
  }
  /* Re-setting number of bines */
  wtr->nbins = i;
  printf("nb of bins: %d\n", wtr->nbins) ;
  realloc_err(wtr->th, wtr->nbins*sizeof(double), err);
  forwardError(*err, __LINE__, NULL);
  
  realloc_err(wtr->w, wtr->nbins*sizeof(double), err);
  forwardError(*err, __LINE__, NULL);
  
  realloc_err(wtr->werr, wtr->nbins*sizeof(double), err);
  forwardError(*err, __LINE__, NULL);
  
  
  fclose(inF);
  return(wtr);
}


/*----------------------------------------------------------------*
 *Likelihood for wp                                               *
 *----------------------------------------------------------------*/

double chi2_mwolk(cosmo_hm *model, const wt_t* data, ngal_fit_t ngal_fit_type,
	       double ngal, double ngalerr, double area, error **err)
{
   double chi2, ngd, ngd_obs, ngd_err, volume, rad2degsqr;

   // NOT USED ANYMORE see chi2_hm(...)

   /* TODO: combine several redshift bins! */

   //Observed 2d-number density of galaxies
   volume     = vc(model,model->zmin,model->zmax,err);  forwardError(*err, __LINE__, -1.0);
   rad2degsqr = (pi/180.0)*(pi/180.0);
   ngd_obs    = ngal/area/rad2degsqr/volume;
   ngd_err    = ngd_obs/ngal*ngalerr;
   
   //double ngd_obs_w;
   //ngd_obs_w = ngd_obs_weighted(model,ngal,area,err);
   //printf ("ngd=%5.2e ngd_err=%5.2e ngd_w=%5.2e\n",ngd_obs, ngd_err, ngd_obs_w);
   
   ngd  = 0.0;
   chi2 = compute_chisq_wp(model, data, ngd_obs, ngd_err, ngal_fit_type, &ngd, 0, err);
   forwardError(*err, __LINE__, 0.0);

   fprintf(stderr, "chi2_hm: ngden(mod,obs,err) = (%5.2e,%5.2e,%5.2e) ln(ngden)(mod,obs,err) = (%g,%g,%g, %g)\n",
   	   ngd, ngd_obs, ngd_err, log(ngd), log(ngd_obs), ngd_err/ngd_obs, log(ngd_err+ngd_obs) - log(ngd_obs-ngd_err));

   /* det C ... */
   return -0.5*chi2;
}



double test(double r, void *param, error **err)
{
   double answer, rp, a;
   cosmo_hm_params *cANDs;
   cosmo_hm *model;

   cANDs = (cosmo_hm_params*)param;
   model  = cANDs->model;
   a     = cANDs->a;
   rp    = cANDs->k;

   answer = r*xir_mwolk(model, a, r, err)*pow((r*r-rp*rp),-0.5);
   forwardError(*err, __LINE__, -1.0);

   return answer;
}

#define EPS 1.0e-4
double wp_mwolk(cosmo_hm *model, double rp, error **err)
{
   double answer2, a;
   cosmo_hm_params cANDs;

   a = 1.0/(1.0+zmean(model->redshift, 0, err));
   forwardError(*err, __LINE__, -1.0);

   cANDs.model = model;
   cANDs.a    = a;
   cANDs.k    = rp;
   //answer2 = 2*sm2_qromberg(test, 0, rp+EPS, 1000, EPS, err);
   answer2 = int_gsl(test, (void*)&cANDs, rp + EPS, 100.0, EPS, err);
   forwardError(*err, __LINE__, -1.0);

   return answer2;
}
#undef EPS

/* ============================================================= *
 * Returns chi^2 				                 *
 * ============================================================= */
double compute_chisq_wp(cosmo_hm *model, const wt_t *wth, double ngd_obs, double ngd_err,
		     ngal_fit_t ngal_fit_type, double *ngd, int dologw, error **err)
{
   int i, j;
   double chisq;
   double dchisq1,dchisq2;
   double lngd, lngd_obs, lngd_err;
   
   testErrorRet(model->redshift->Nzbin!=1, ce_overflow,
		"More than one redshift bin not yet supported in likelihood",
		*err, __LINE__, -1.0);
   
   if (ngal_fit_type!=ngal_lin_fit_only) {
    
     double *rp = malloc_err(sizeof(double)*wth->nbins, err);  forwardError(*err, __LINE__, -1.0);
     double *wh    = malloc_err(sizeof(double)*wth->nbins, err);  forwardError(*err, __LINE__, -1.0);
     
     /* the option model->FFTLog = 0 does not exist anymore. To be deprecaded... */
     model->FFTLog = 1;
     
     for (i=0; i<wth->nbins; i++) {
       rp[i] = wth->th[i];
       wh[i] = wp_mwolk(model, rp[i], err);
       forwardError(*err, __LINE__, -1.0);
     }   
     
     if (wth->wcov==NULL) {
       for (i=0,chisq=0.0; i<wth->nbins; i++) {
	 dchisq1  = dsqr((wth->w[i] - wh[i])/(wth->werr[i]));
	 chisq   += dchisq1;
	 fprintf (stderr, "%10d %5.2e %5.2e %5.2e %g\n", i, wh[i],wth->w[i],dchisq1, chisq);
       }
     } else {
       for (i=0,chisq=0.0; i<wth->nbins; i++) {
	 for (j=0; j<wth->nbins; j++) {
	   dchisq1  = (wth->w[i] - wh[i])*wth->wcov[i*wth->nbins+j]*(wth->w[j] - wh[j]);
	   testErrorRet(!finite(wth->wcov[i*wth->nbins+j]), ce_infnan, "inf or nan in logl", *err, __LINE__, -1.0);
	   chisq   += dchisq1;
	   if (i==j) fprintf (stderr, "%10d %10d %5.2e  %5.2e %5.2e   %5.2e %g\n",
			      i, j, wh[i], wth->w[i], wth->wcov[i*wth->nbins+j], dchisq1, chisq);
	}
       }
     }
     
     free(wh);
     free(rp);
     
   } else {
     
     /* Only fit number of galaxies (below) */
     chisq = 0.0;
     
   }
   
  
   if (ngal_fit_type != ngal_no_fit) {
     *ngd = ngal_den_vol(model, err);
  } else {
     *ngd = -1.0;
   }
   
   /* Number density of galaxies */
   switch (ngal_fit_type) {
   case ngal_log_fit :
     /* "Number density varied logarithmically with M1" [Hamana 2004] */
     lngd_obs = log(ngd_obs);
     
     /* Delta ln x = Delta x / x */
     lngd_err = ngd_err/ngd_obs;
     
     lngd = log(*ngd);
    
    dchisq2 = dsqr((lngd_obs - lngd)/lngd_err);
    break;

  case ngal_lin_fit : case ngal_lin_fit_only :
    dchisq2 = dsqr((ngd_obs - *ngd)/ngd_err);
    break;

  case ngal_no_fit :
    dchisq2 = 0.0;
    break;
    
    /* case ngal_match_M1 : */
    /* MKDEBUG: TODO... */
    /*--------------------------------------------------------------------//
      [jean]  We could use the formula given in Brown et al. 2008 (eq. 11)
      update (Feb 2013) useless - parameter removed
      //--------------------------------------------------------------------*/
  default :
    *err = addErrorVA(ce_unknown, "Wrong ngal_fit_type %d", *err, __LINE__, (int)ngal_fit_type);
    return -1.0;
  }
  
  chisq += dchisq2;
  
  fprintf(stderr, "\nCS2 chi^2: %g %g\n", chisq, dchisq2);
  //fprintf(stderr, "ngal_den(obs) ngal_den(pred): %5.2e %5.2e\n", exp(lngd_obs), exp(lngd));
  
  testErrorRetVA(chisq<0.0, math_negative, "Negative chi^2 %g. Maybe the covariance matrix is not positive",
		 *err, __LINE__, -1.0, chisq);
  testErrorRet(!finite(dchisq2), ce_infnan, "inf or nan in logl (dchisq2)", *err, __LINE__, -1.0);
  testErrorRet(!finite(chisq), ce_infnan, "inf or nan in logl", *err, __LINE__, -1.0);
  
  return chisq;
}


/*----------------------------------------------------------------*
 *End                                                             *
 *----------------------------------------------------------------*/
