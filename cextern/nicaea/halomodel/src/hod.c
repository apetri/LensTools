/* ============================================================ *
 * hod.c							*
 * Martin Kilbinger, Henry J. McCracken, Jean Coupon 2008-2012  *
 * Refs:							*
 *   - Coupon et al. 2011 (arXiv:1107.0616)                     *
 *   - Brown et al. 2008 (ApJ 682, 937)		        	*
 *   - Zheng et al. 2005 (ApJ, 633, 791)         		*
 *   - Tinker et al. 2005 (ApJ, 631, 41)	        	*
 *   - Kravstov et al. 2004 (ApJ, 609, 35)    	        	*
 *   - Hamana et al 2004 (MNRAS 347, 813)			*
 *   - Zheng et al. 2005 (ApJ, 633, 791)            		*
 *   - Berlind & Weinberg 2002 (ApJ 575, 587)			*
 *                                                              *
 * WARNINGS:                                                    *
 * - The Hamana model is currently NOT supported                *
 * - the FFTLog routines have not been tested in case           *
 *   of a change of cosmological parameters during a pmc run.   *
 *   Only tested for changes in HOD parameters.                 *
 * - sm2_romberg should never be used with eps > 1e-4           *
 *   (use int_gsl instead)                                      *
 * ============================================================ */

#include "hod.h"

double FFTLog_TMP;
#define FFTLog_SWAP(a,b) {FFTLog_TMP = (a); (a) = (b); (b) = FFTLog_TMP;}

/*----------------------------------------------------------------*
 *Add by Melody Wolk 18/07/2011                                   *
 *----------------------------------------------------------------*/

/*----------------------------------------------------------------*
 *Compute total xiofr                                             *
 *----------------------------------------------------------------*/

#define EPS 1.0e-5
double xir_mwolk(cosmo_hm *self, double a, double r, error **err)
{
   double dr, logr, res;
   int j, nbins=60;
   
   testErrorRetVA(self->hod!=berwein02 && self->hod!=berwein02_hexcl, hm_hodtype,
		  "Invalid hod type (%d), only berwein02(%d) and berwein02_hexcl(%d) hod implemented for xi(r)",
		  *err, __LINE__, 0.0, self->hod, berwein02, berwein02_hexcl);
   
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
   if (fabs(self->a_xir-a)>EPS) {
     del_interTable(&self->xir);
     //  printf("Check2 %lf %lf\n", self->a_xir, a);
   }

   double *ri =malloc(nbins*sizeof(double));
   double *xi =malloc(nbins*sizeof(double));
   
   if (self->xir==NULL) {
     dr       = (log10(RMAX)-log10(RMIN))/((double)nbins);
     self->xir   = init_interTable(nbins, log10(RMIN), log10(RMAX), dr, 0.0, 0.0, err);
     forwardError(*err, __LINE__, 0.0);
     self->a_xir = a;
     
     
     /* Earlier: Tabulation in a as well */
     
     for (j=0; j<nbins; j++) {
       ri[j] = pow(10., log10(RMIN)+(double)j*dr);
     }
     res=0.;
     double *xi1hgcs= xiofr(self, p1hgcs, ri, nbins, err);     forwardError(*err, __LINE__, -1.0);
     double *xi1hgss= xiofr(self, p1hgss, ri, nbins, err);     forwardError(*err, __LINE__, -1.0);
     double *xi2hg= xiofr(self, p2hg, ri, nbins, err);         forwardError(*err, __LINE__, -1.0);


     for (j=0; j<nbins; j++) {
       xi[j]= xi1hgcs[j]+xi1hgss[j]+xi2hg[j];
     }
     
     for (j=0; j<nbins; j++) {
       self->xir->table[j]=xi[j];
       // printf("xir: %lf %lf %d\n",self->xir->table[j], pow(10., log10(RMIN)+(double)j*dr), j);
     }


   }

   if (r<RMIN || r>RMAX) return 0.0;

   logr=log10(r) ;
  
   res = interpol_wr(self->xir, logr, err);
   forwardError(*err, __LINE__, 0.0);
   //printf("test: %lf %lf\n", self->xir->upper, self->xir->lower) ;
   // printf("r:%lf xir:%lf\n", r, res) ;
   return res;
}
#undef EPS


/*----------------------------------------------------------------*
 *Read mwolk format files                                         *
 *----------------------------------------------------------------*/

wt_t *read_wtheta_mwolk(const char *name, double delta, double intconst, error **err)
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

   /* TODO: combine several redshift bins! */

   /* Observed 2d-number density of galaxies */
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
   cosmo_hmANDstuff *cANDs;
   cosmo_hm *self;

   cANDs = (cosmo_hmANDstuff*)param;
   self  = cANDs->self;
   a     = cANDs->a;
   rp    = cANDs->k;

   answer = r*xir_mwolk(self, a, r, err)*pow((r*r-rp*rp),-0.5);
   forwardError(*err, __LINE__, -1.0);

   return answer;
}

#define EPS 1.0e-4
double wp(cosmo_hm *self, double rp, error **err)
{
   double answer2, a;
   cosmo_hmANDstuff cANDs;

   a = 1.0/(1.0+zmean(self->redshift, 0, err));
   forwardError(*err, __LINE__, -1.0);

   cANDs.self = self;
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
   
    
    model->FFTLog = 1;

 
    for (i=0; i<wth->nbins; i++) {
      rp[i] = wth->th[i];
      wh[i] = wp(model, rp[i], err);
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
  
  
  if (ngal_fit_type!=ngal_no_fit) {
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

  case ngal_match_M1 :
    /* MKDEBUG: TODO... */
    /*--------------------------------------------------------------------//
      [jean]  We could use the formula given in Brown et al. 2008 (eq. 11)
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


double chi2_DeltaSigma(cosmo_hm *model, halomode_t halomode, const wt_t* data, error **err){
  /* Returns chi^2 using Delta Sigma. */
  
  int i, j;
  double chi2;
    
  double z = (model->zmax + model->zmin)/2.0;
  
  /* Delta Sigma */
  double *deltaSigma = malloc_err(sizeof(double)*data->nbins, err);  forwardError(*err, __LINE__, -1.0);
  for (i=0; i< data->nbins; i++){
    deltaSigma[i] = DeltaSigma(model, data->th[i], pow(10.0, model->log10Mhalo), 1.0/(1.0+z), err);
  }
  
  /* chi2 */
  chi2 = 0.0;
  if (data->wcov == NULL){
    for (i=0; i < data->nbins; i++) chi2 += dsqr((data->w[i] - deltaSigma[i])/(data->werr[i]));
  }else{
    for (i=0; i < data->nbins; i++) {
      for (j=0; j < data->nbins; j++) {
	chi2 += (data->w[i] - deltaSigma[i])*data->wcov[i*data->nbins+j]*(data->w[j] - deltaSigma[j]);
	testErrorRet(!finite(data->wcov[i*data->nbins+j]), ce_infnan, "inf or nan in logl", *err, __LINE__, -1.0);
      }
    }
  }
  free(deltaSigma);
  
  testErrorRetVA(chi2<0.0, math_negative, "Negative chi^2 %g. Maybe the covariance matrix is not positive",*err, __LINE__, -1.0, chi2);
  testErrorRet(!finite(chi2), ce_infnan, "inf or nan chi2", *err, __LINE__, -1.0);
  
  /* det C ... */
  return -0.5*chi2;
}




double chi2_hm(cosmo_hm *model, halomode_t halomode, const wt_t* data, ngal_fit_t ngal_fit_type,
	       double ngal, double ngalerr, double area, error **err)
{
  /* 
     Returns chi^2 using the angular correlation function w(theta) and the galaxy number density.
  */
  int i, j;
  double chi2, ngd, ngd_obs, ngd_err, volume, rad2degsqr, zm;
  
  /* To do */
  testErrorRet(model->redshift->Nzbin!=1, ce_overflow,
	       "More than one redshift bin not yet supported in likelihood",
	       *err, __LINE__, -1.0);
  
  /* Observed 2d-number density of galaxies */
  volume     = vc(model,model->zmin,model->zmax,err);  forwardError(*err, __LINE__, -1.0);
  rad2degsqr = (pi/180.0)*(pi/180.0);
  ngd_obs    = ngal/area/rad2degsqr/volume;
  ngd_err    = ngd_obs/ngal*ngalerr;  
  zm         = zmean(model->redshift,0, err);          forwardError(*err, __LINE__, -1.0);
  ngd        = ngal_den(model,1.0/(1.0+zm),logMmax,err);   forwardError(*err, __LINE__, -1.0);
  
  /* Only if ngal used in chi2 */
  if(ngal_fit_type == ngal_lin_fit_only){
    chi2 = dsqr((ngd_obs - ngd)/ngd_err);
    testErrorRetVA(chi2<0.0, math_negative, "Negative chi^2 %g",*err, __LINE__, -1.0, chi2);
    testErrorRet(!finite(chi2), ce_infnan, "inf or nan chi^2", *err, __LINE__, -1.0);
    return -0.5*chi2;
  }
  
  /* w(theta) */
  double *wh     = malloc_err(sizeof(double)*data->nbins, err);  forwardError(*err, __LINE__, -1.0);
  double *w1hgcs = woftheta_FFTLog(model,p1hgcs,data->th,data->nbins,err); forwardError(*err, __LINE__, -1.0);
  double *w1hgss = woftheta_FFTLog(model,p1hgss,data->th,data->nbins,err); forwardError(*err, __LINE__, -1.0);
  double *w2hg   = woftheta_FFTLog(model,p2hg,data->th,data->nbins,err);   forwardError(*err, __LINE__, -1.0);
  for (i=0; i< data->nbins; i++) wh[i] = w1hgcs[i] + w1hgss[i] + w2hg[i];
  free(w1hgcs);free(w1hgss);free(w2hg);

  /* log(w(theta)) */
  if (halomode==galcorr_log) for (i=0; i< data->nbins; i++) wh[i] = log(wh[i]);
  
  /* chi2 */
  chi2 = 0.0;
  if (data->wcov==NULL){
    for (i=0; i < data->nbins; i++) chi2 += dsqr((data->w[i] - wh[i])/(data->werr[i]));
  }else{
    for (i=0; i < data->nbins; i++) {
      for (j=0; j < data->nbins; j++) {
	chi2 += (data->w[i] - wh[i])*data->wcov[i*data->nbins+j]*(data->w[j] - wh[j]);
	testErrorRet(!finite(data->wcov[i*data->nbins+j]), ce_infnan, "inf or nan in logl", *err, __LINE__, -1.0);
      }
    }
  }
  free(wh);
  
  /* Number density of galaxies */
  switch (ngal_fit_type) {
  case ngal_log_fit : /* "Number density varies logarithmically with M1" [Hamana 2004] */
    chi2     += dsqr((log(ngd_obs) - log(ngd))/(ngd_err/ngd_obs)); /* Delta ln x = Delta x / x */
    break;
  case ngal_lin_fit :
    chi2     += dsqr((ngd_obs - ngd)/ngd_err);
  break;
  case ngal_no_fit :
    break;
  case ngal_match_M1 :
    /* MKDEBUG: TODO */
  default :
    *err = addErrorVA(ce_unknown, "Wrong ngal_fit_type %d", *err, __LINE__, (int)ngal_fit_type);
    return -1.0;
  }
  
  if(ngal_fit_type != ngal_lin_fit_only){
    fprintf(stderr, "chi2_hm: ngden(mod,obs,err) = (%5.2e,%5.2e,%5.2e) ln(ngden)(mod,obs,err) = (%g,%g,%g, %g)\n",
	    ngd, ngd_obs, ngd_err, log(ngd), log(ngd_obs), ngd_err/ngd_obs, log(ngd_err+ngd_obs) - log(ngd_obs-ngd_err));
  }
  
  testErrorRetVA(chi2<0.0, math_negative, "Negative chi^2 %g. Maybe the covariance matrix is not positive",*err, __LINE__, -1.0, chi2);
  testErrorRet(!finite(chi2), ce_infnan, "inf or nan chi2", *err, __LINE__, -1.0);
  
  /* det C ... */
  return -0.5*chi2;
}

/*----------------------------------------------------------------*
 *w(theta)                                                        *
 *----------------------------------------------------------------*/

double *woftheta_FFTLog(cosmo_hm *model, pofk_t pofk, double *theta, int Ntheta, error **err)
{
  /*First compute xi, project using Limber equation at mean z and return w(theta). 
    See Bartelmann & Schneider (2001) eq. (2.79), Tinker et al. 2010 eq. (3), 
    Ross et al. 2009 eq (27), ...
    Theta in degree.
    N is the number of points that sample xi.The speed of the code
    is basically inversely proportional to this number...choose it carefuly!!
    You can also play with N in FFTLog routines.
  */
  
  double *result  = malloc_err(Ntheta*sizeof(double),err); 
  forwardError(*err, __LINE__, NULL);
   
  //tabulate xi(r)
  int i,j,k,N      = 40;
  double *u        = malloc_err(N*sizeof(double), err);  forwardError(*err, __LINE__, NULL);
  double *logu     = malloc_err(N*sizeof(double), err);  forwardError(*err, __LINE__, NULL);
  double umin      = 0.001, umax = 800;
  double dlogu     = log(umax/umin)/(double)N;

  for(i=0;i<N;i++){
    logu[i] = log(umin)+dlogu*(double)i;
    u[i]    = exp(logu[i]);
  }
  
  double *xi = xiofr(model,pofk,u,N,err);
  forwardError(*err, __LINE__, NULL);

  //interpolate xi(r)
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline    = gsl_spline_alloc (gsl_interp_cspline,N);
  gsl_spline_init(spline,logu,xi,N);
  
  Nz_hist *nz = get_nz(model->redshift,err);
  forwardError(*err, __LINE__, NULL);
  
  int wOmegar           = 0;
  double deg_to_rad_sqr = dsqr(pi/180.0);
  double ww,r,x,sum;
  
  //Limber equation - project xi to get w
  for(i=0;i<Ntheta;i++){ //loop over theta
    result[i] = 0.0;
    for(j=0;j<nz->nbins;j++){ //loop over z
      ww  = w(model->cosmo,1.0/(1.0+nz->z[j]), wOmegar, err); 
      forwardError(*err, __LINE__, NULL);
      x   = f_K(model->cosmo, ww, err);         
      forwardError(*err, __LINE__, NULL);
      sum = 0.0;
      for(k=0;k<N;k++){ //loop over u
	r    = sqrt(u[k]*u[k] + x*x*theta[i]*theta[i]*deg_to_rad_sqr);
	if(log(r) < logu[N-1]) {//to unsure log(r) lies within interpolation limits
	  sum += u[k]*gsl_spline_eval(spline,log(r),acc);
	}
	}
      result[i] += dsqr(nz->n[j])/drdz(model->cosmo,1.0/(1.0+nz->z[j]),err)*sum;
      forwardError(*err, __LINE__, NULL);
    }
    result[i] *= 2.0*nz->dz*dlogu;
    testErrorRetVA(!finite(result[i]), ce_infnan, "inf or nan in w(theta_%d=%g)",
		   *err, __LINE__, NULL, i, theta[i]/arcmin);
  }
  
  free(xi);
  free(u);
  free(logu);
  free_nz(nz);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  
  return result;
}

/* Returns an array of w(theta) for the sum of the three HOD terms (1hcs + 1hss + 2h) */
double* woftheta_FFTLog_total(cosmo_hm *model, double log_theta_min, double log_delta_theta, int Ntheta, error **err)
{
   int i;
   double *theta, *wh, *w1hgcs, *w1hgss, *w2hg;

   theta = malloc_err(sizeof(double)*Ntheta, err);  forwardError(*err, __LINE__, NULL);
   wh    = malloc_err(sizeof(double)*Ntheta, err);  forwardError(*err, __LINE__, NULL);

   for (i=0; i<Ntheta; i++) theta[i] = pow(10.0, log_theta_min + i*log_delta_theta);

   w1hgcs = woftheta_FFTLog(model, p1hgcs, theta, Ntheta, err); forwardError(*err, __LINE__, NULL);
   w1hgss = woftheta_FFTLog(model, p1hgss, theta, Ntheta, err); forwardError(*err, __LINE__, NULL);
   w2hg   = woftheta_FFTLog(model, p2hg,   theta, Ntheta, err); forwardError(*err, __LINE__, NULL);

  for (i=0; i<Ntheta; i++) wh[i] = w1hgcs[i] + w1hgss[i] + w2hg[i];

  free(w1hgcs); free(w1hgss); free(w2hg); free(theta);

  return wh;
}

/*----------------------------------------------------------------*
 *HOD P(k) and xi(r)                                              *
 *----------------------------------------------------------------*/

double *xiofr(cosmo_hm *model, pofk_t pofk, double *r, int N, error **err){
  /*Compute xi(r) at mean redshift*/
  

   double *res;
   double a = 1.0/(1.0 + zmean(model->redshift, 0, err));
   forwardError(*err, __LINE__, NULL);
  
   //q = 0.0, m = 1/2
   switch (pofk) {
      case p1hgcs :
	 res = xi_1hcs(model,a,r,N,err);
	 forwardError(*err, __LINE__, NULL);
	 break;
      case p1hgss :
	 res = xi_1hss(model,a,r,N,err);
	 forwardError(*err, __LINE__, NULL);
	 break;
      case p2hg   :
	 res = xi_2h(model,a,r,N,err);
	 forwardError(*err, __LINE__, NULL);
	 break;
      case pnl    :
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

double *xi_1hcs(cosmo_hm *model,double a,double *r, int N, error **err){
  /*1-halo central-satellite*/

  int i;
  double *result = malloc_err(N*sizeof(double), err); forwardError(*err, __LINE__, NULL);
  
  cosmo_hm_params params;
  params.cosmo  = model->cosmo;
  params.model  = model;
  params.a      = a;
  params.ng     = ngal_den(model,a,logMmax,err);
  forwardError(*err, __LINE__, NULL);
  params.err    = err;
  params.eps    = 1.0e-3;
  
  for(i=0;i<N;i++){
    if(r[i] > RMAX1){
      result[i] = 0.0;
    }else{
      params.r = r[i];
      result[i] = int_gsl(int_for_xi_1hcs, (void*)&params, log(M_vir(model,r[i],a)),
				 logMmax, params.eps, err);
      forwardError(*err, __LINE__, NULL);
    }
  }
  return result;
}

double int_for_xi_1hcs(double logM, void *params, error **err)
{
  double res, r_vir;
  double M         = exp(logM);
  cosmo_hm *model  =  ((cosmo_hm_params *)params)->model;
  double a         =  ((cosmo_hm_params *)params)->a; 
  double r         =  ((cosmo_hm_params *)params)->r;
  double ng        =  ((cosmo_hm_params *)params)->ng;
  
  //NFW profile, -1.0 is to re-compute r_vir
  r_vir = -1.0;
  double rho = rho_halo(model, r, M, a,&r_vir, err);
  
  //Hala mass function parameters
  cosmo_hmANDstuff2 intpar;
  intpar.self       = model;
  intpar.a          = a;
  intpar.asymptotic = 0;
  
  res  = dn_dlnM(M,(void*)&intpar,err);  forwardError(*err, __LINE__, 0.0);
  res *= Ngal_c(model,M)*Ngal_s(model,M)/(0.5*ng*ng)*rho/M;
  return res;
}

double *xi_1hss(cosmo_hm *model,double a,double *r, int N, error **err){
  /*1-halo satellite-satellite*/
  
  int i,j;
  double *result = malloc_err(N*sizeof(double), err);  forwardError(*err, __LINE__, NULL);
  
  cosmo_hm_params params;
  params.cosmo = model->cosmo;
  params.model = model;
  params.a     = a;
  params.ng    = ngal_den(model,a,logMmax,err);
  forwardError(*err, __LINE__, NULL);
  params.err   = err;
  params.eps   = 1.0e-5;
  
  FFTLog_config *fc = FFTLog_init(128,k_min,k_max_HOD,0.0,0.5);
  
  double *r_FFT    = malloc_err(fc->N*sizeof(double),err);  forwardError(*err, __LINE__, NULL);
  double *ar       = malloc_err(fc->N*sizeof(double),err);  forwardError(*err, __LINE__, NULL);
  double *logr_FFT = malloc_err(fc->N*sizeof(double),err);  forwardError(*err, __LINE__, NULL);
 
  gsl_function Pk;
  Pk.function = &P1hss;
  Pk.params   = &params;

  //FFTlog
  FFTLog(fc,&Pk,r_FFT,ar,-1,err); forwardError(*err, __LINE__, NULL);
  
  //Interpolation
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline    = gsl_spline_alloc (gsl_interp_cspline,fc->N);
  //Attention: N and fc->N are different
  for(j=0;j<fc->N;j++) logr_FFT[j] = log(r_FFT[j]);
  gsl_spline_init (spline,logr_FFT,ar,fc->N);
  
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
  double ng        = ((cosmo_hm_params *)params)->ng;
  double eps       = ((cosmo_hm_params *)params)->eps;
  error **err      = ((cosmo_hm_params *)params)->err;

  res = pow(k,1.5)*int_gsl(int_for_P1hss,params,logMmin,logMmax,eps,err)/(ng*ng);
  forwardError(*err, __LINE__, 0.0);
  
  return res;
}

double int_for_P1hss(double logM, void *params, error **err)
{
  double res;
  cosmo_hm *model  = ((cosmo_hm_params *)params)->model; 
  double a         = ((cosmo_hm_params *)params)->a; 
  double k         = ((cosmo_hm_params *)params)->k;
  double M         = exp(logM);
  
  //Fourrier transform of halo profile 
  double rhohat =  rhohat_halo(model, k, M, a, 1, err);
  forwardError(*err, __LINE__, 0.0);

  //Hala mass function parameters
  cosmo_hmANDstuff2 intpar;
  intpar.self       = model;
  intpar.a          = a;
  intpar.asymptotic = 0;

  res  = dn_dlnM(M,(void*)&intpar,err);
  forwardError(*err, __LINE__, 0.0);
  //res *= Ngal_c(model,M)*dsqr(Ngal_s(model,M))*rhohat*rhohat;
  res *= dsqr(Ngal_s(model,M))*rhohat*rhohat;
  
  return res;
}

#define EPS 1.0e-5
double *xi_2h(cosmo_hm *model,double a,double *r, int N, error **err)
{
  int i,j;
  double *result = malloc_err(N*sizeof(double), err);  forwardError(*err, __LINE__, NULL);
  
  double xi_dm, amean;
  cosmo_hm_params params;
  params.cosmo  = model->cosmo;
  params.model  = model;
  params.a      = a;
  params.ng     = ngal_den(model,a,logMmax,err);
  forwardError(*err, __LINE__, NULL);
  params.err    = err;
  params.eps    = 1.0e-3;
  
  FFTLog_config *fc;
  gsl_function Pk;

  /* New: consistency check */
  amean = 1.0 / (1.0 + zmean(model->redshift, 0, err));
  testErrorRetVA(fabs(amean - a) > EPS, hm_zmean_2h, "Precalculated xi_dm at a=%g not consistent with a=%g",
		 *err, __LINE__, NULL, amean, a);

  //Tabulate xi_dm
  if (model->xi_dm==NULL) {
    fc = FFTLog_init(128,k_min,k_max_HOD,0.0,0.5);
    double *r_FFT     = malloc_err(fc->N*sizeof(double),err);  forwardError(*err, __LINE__, NULL);
    double *ar        = malloc_err(fc->N*sizeof(double),err);  forwardError(*err, __LINE__, NULL);
    Pk.function       = &FFTLog_P_NL;
    Pk.params         = &params;
    FFTLog(fc,&Pk,r_FFT,ar,-1,err); forwardError(*err, __LINE__, NULL);
    double dlogr       = (log(r_FFT[fc->N-1]) - log(r_FFT[0]))/(double)fc->N;
    model->xi_dm       = init_interTable(fc->N,log(r_FFT[0]),log(r_FFT[fc->N-1]),dlogr, 0.0, 0.0, err);
    forwardError(*err, __LINE__, NULL);
    model->a_xir       = a;
    //fill in the table
    for(j=0;j<fc->N;j++) model->xi_dm->table[j] = ar[j]*pow(2.0*pi*r_FFT[j],-1.5);
    //free memory
    free(r_FFT);
    free(ar);
    FFTLog_free(fc);
  }
  
  //xi 2h
  Pk.function = &P2h;
  Pk.params = &params;
  fc = FFTLog_init(64,k_min,k_max_HOD,0.0,0.5);

  for(i=0;i<N;i++){
    if(r[i] < RMIN2){
      result[i] = 0.0;
    }else{
      if (model->hod == berwein02_hexcl || model->hod == leauthaud11){
	params.logMlim    = logM_lim(model,a,r[i],err);
	forwardError(*err, __LINE__, NULL);
	params.ngp        = ngal_den(model,a,params.logMlim,err);
	forwardError(*err, __LINE__, NULL);
	params.bias_func  = &bias_r;
	//xi_dm             = xi_dm_NL(model,a,r[i],err); // test(slower)
	xi_dm             = interpol_wr(model->xi_dm, log(r[i]), err);
	forwardError(*err, __LINE__, NULL);
	params.bias_fac   = pow(1.0+1.17*xi_dm,1.49)/pow(1.0+0.69*xi_dm,2.09);
      }else if (model->hod == berwein02){
	params.logMlim   = logMmax;
	params.ngp       = params.ng;
	params.bias_func = &bias_r;
	params.bias_fac  = 1;    /* No scale-dependent bias */
      }      
      
      if(params.ngp < 1.0e-8){
	result[i] = 0.0;
      }else{
	result[i] = dsqr(params.ngp/params.ng)*xi_from_Pkr(&Pk,r[i],fc,err);
	forwardError(*err, __LINE__, NULL);
      }
      
    }
  }
  FFTLog_free(fc);
  
  return result;
}
#undef EPS

/* Angular correlation function of dark matter */
#define EPS 1.0e-5
double *xi_P_NL(cosmo_hm *model, double a, double *r, int N, error **err)
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
  fc = FFTLog_init(64,k_min,k_max_HOD,0.0,0.5);
  
  for(i=0;i<N;i++){
    result[i] = xi_from_Pkr(&Pk, r[i], fc, err);
    forwardError(*err, __LINE__, NULL);
  }
  FFTLog_free(fc);  
  return result;
}

#undef EPS

double FFTLog_P_NL(double k, void *params)
{
  double res;
  ((cosmo_hm_params *)params)->k = k;
  cosmo_hm *model  = ((cosmo_hm_params *)params)->model; 
  double a         = ((cosmo_hm_params *)params)->a;
  error **err      = ((cosmo_hm_params *)params)->err;
  
  res = pow(k,1.5)*P_NL(model->cosmo,a,k,err);
  forwardError(*err, __LINE__, 0.0);
  return res;
}

double P2h(double k, void *params)
{
  double res;  
  ((cosmo_hm_params *)params)->k = k;
  cosmo_hm *model  = ((cosmo_hm_params *)params)->model; 
  double a         = ((cosmo_hm_params *)params)->a;
  double ngp       = ((cosmo_hm_params *)params)->ngp;
  double eps       = ((cosmo_hm_params *)params)->eps;
  double logMlim   = ((cosmo_hm_params *)params)->logMlim;
  error **err      = ((cosmo_hm_params *)params)->err;
  
  res  = pow(k,1.5)*P_NL(model->cosmo,a,k,err);
  forwardError(*err, __LINE__, 0.0);
  res *= dsqr(int_gsl(int_for_P2h,params,logMmin,logMlim,eps,err)/ngp);
  forwardError(*err, __LINE__, 0.0);
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
  
  //Fourrier transform of halo profile 
  double rhohat =  rhohat_halo(model, k, M, a, 1, err);
  forwardError(*err, __LINE__, 0.0);

  //Hala mass function parameters
  cosmo_hmANDstuff2 intpar;
  intpar.self       = model;
  intpar.a          = a;
  intpar.asymptotic = 0;
  
  double b = bias_func(M,params);
  
  res = dn_dlnM(M,(void*)&intpar,err)*b*Ngal(model,M)*rhohat;
  forwardError(*err, __LINE__, 0.0);
  return res;
}

/* ============================================================ *
 * Tinker 2005 eq B7. Depends on r via bias_fac.                *
 * Used to be called "bias_tinker_r". Now call to bias_sc is    *
 * allowed depending on self->halo_bias.			*
 * ============================================================ */
double bias_r(double M, void *params)
{
  double b;
  cosmo_hm *model  = ((cosmo_hm_params *)params)->model; 
  double a         = ((cosmo_hm_params *)params)->a;
  double bias_fac  = ((cosmo_hm_params *)params)->bias_fac;
  error **err      = ((cosmo_hm_params *)params)->err;
  
  // MKDEBUG: NEW
  //b  = bias_tinker(model,M,a,err);
  b = halo_bias(model, M, a, 1, err);
  forwardError(*err, __LINE__, 0.0);
  
  b *= sqrt(bias_fac);

  return b;
}

/*----------------------------------------------------------------*
 *HOD model functions                                             *
 *----------------------------------------------------------------*/

double log_fSHMR(double log10Mh, cosmo_hm *model){
  /* TO DO: interpolation + limits*/
  
  int status;
  int iter = 0, max_iter = 100;
  double x, x_lo = logMmin - log10Mh, x_hi = logMmax - log10Mh;
  
  gsl_function F;
  F.function = &log_fSHMR_inv_minus_Mh;
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
      status = gsl_root_test_interval (x_lo, x_hi, 0, 0.001);
    }
  while (status == GSL_CONTINUE && iter < max_iter);
  
  gsl_root_fsolver_free (s);
  
  return x;
}


double log_fSHMR_inv_minus_Mh(double log10Mstar, void *p){
  cosmo_hm *model = (cosmo_hm *)p;
  double result;
  
  
  double A = pow(10.0,log10Mstar)/model->Mstar0;
  result   = log10(model->M1) + model->beta*log10(A);
  result  += pow(A,model->delta)/(1.0 + pow(A,-model->gamma)) - 0.5;
  
  return result - model->x;
}


double Ngal_c(cosmo_hm *model, double M){
  /*Number of central galaxies per halo*/
  double arg, result;
  
  if(M < 10.0e10) return 0.0;

  if (model->hod == leauthaud11){
    arg    = (log10(model->Mstar_thres) - log_fSHMR(log10(M),model))/(sqrt(2.0)*model->sigma_log_M);
    result = 0.5*(1.0-gsl_sf_erf(arg));
  }else{
    arg    = (log10(M/model->M_min)/model->sigma_log_M);
    result =  0.5*(1.0+gsl_sf_erf(arg));
  }      
  return result;
}


double Ngal_s(cosmo_hm *model, double M){
  /*
    Number of satellite galaxies per halo. Two solutions here:
    1. N(M) = Ncen(1+Nsat) with Nsat \prop M^alpha and
    2. N(M) = Ncen + Nsat  with Nsat \prop Ncen * M^alpha 
    These are equivalent but <N(N-1)> must be changed accordinlgy. 
    Now: N(M) = Ncen(1+Nsat) because xi_cs and xi_ss are consistent 
    with this expression (see Blake et al. 2008).
    It seems to change quite a lot for unusual sigma_log_M
    
    As in Brown et al. (2008): Ngal_c(model, M)*pow((M-M0)/model->M1, model->alpha);
    As in Blake et al. (2008): pow(M/model->M1,model->alpha);
    As in Zheng et al. (2005): pow((M - M0)/model->M1,model->alpha);    
  */  
  double M0, result;
  
  if(M < 10.0e10) return 0.0;

  if (model->hod == leauthaud11){
    model->x = 0.0; /* to set Mh = 0 in log_fSHMR_inv_minus_Mh, so it computes log_fSHMR_inv */
    double log10Mh_thres = log_fSHMR_inv_minus_Mh(log10(model->Mstar_thres),model);
    double log10Msat     = log10(model->B_sat) + model->beta_sat*log10Mh_thres + 12.0*(1.0 - model->beta_sat);
    double log10Mcut     = log10(model->B_cut) + model->beta_cut*log10Mh_thres + 12.0*(1.0 - model->beta_cut);
    result = pow(M/pow(10.0,log10Msat),model->alpha_sat)*exp(-pow(10.0,log10Mcut)/M);
  }else{
    if (model->M0 < 0) {     /* Negative M0: set to M_min */
      M0 = model->M_min;
    } else {
      M0 = model->M0;
    }
    if (M - M0 < 0.0) return 0.0;
    result = pow((M - M0)/model->M1,model->alpha);    
  }      
  return result;  
}

double Ngal(cosmo_hm *model, double M){
  /*
    Total number of galaxies per halo.
    Zheng et al. (2005): Ngal_c(model, M) + Ngal_s(model, M);
  */
  return Ngal_c(model,M)*(1.0+Ngal_s(model,M));
  
}

/*----------------------------------------------------------------*
 *Deduced quantities                                              *
 *----------------------------------------------------------------*/

double av_gal_bias(cosmo_hm *model, double a, error **err)
{
  /*returns the average galaxy bias wrt to redshift (depends on HOD).*/
 
  double ng,res;
  cosmo_hmANDstuff2 intpar;
  intpar.self = model;
  intpar.a    = a;
  intpar.asymptotic = 0;

  double min = logMmin;

  ng = ngal_den(model,a,logMmax,err);
  forwardError(*err, __LINE__, 0.0);
  res = int_gsl(int_for_av_gal_bias,(void*)&intpar,min,logMmax,1.e-3,err);
  forwardError(*err, __LINE__, 0.0);
  
  return res/ng;
}
double int_for_av_gal_bias(double logM, void *intpar, error **err)
{
  double b,res;
  double M = exp(logM);
  
  cosmo_hm *model  = ((cosmo_hmANDstuff2 *)intpar)->self; 
  double        a  = ((cosmo_hmANDstuff2 *)intpar)->a; 
  
  b   = halo_bias(model, M, a, 1, err); 
  res = dn_dlnM(M,intpar,err)*b*Ngal(model,M);
  forwardError(*err, __LINE__, 0.0);
  
  return res;
}

double av_halo_mass(cosmo_hm *model, double a, error **err)
{
  /*returns the mean halo mass wrt to redshift (depends on HOD).*/
 
  double ng,res;
  cosmo_hmANDstuff2 intpar;
  intpar.self = model;
  intpar.a    = a;
  intpar.asymptotic = 0;

  double min = logMmin;

  ng = ngal_den(model,a,logMmax,err);
  forwardError(*err, __LINE__, 0.0);
  res = int_gsl(int_for_av_halo_mass,(void*)&intpar,min,logMmax,1.e-3,err);
  forwardError(*err, __LINE__, 0.0);

  testErrorRet(ng == 0, ce_infnan, "Division by zero (ng)", *err, __LINE__, 0.0);
  
  return res/ng;
}
double int_for_av_halo_mass(double logM, void *intpar, error **err)
{
  double res;
  double M         = exp(logM);
  cosmo_hm *model  = ((cosmo_hmANDstuff2 *)intpar)->self; 
  
  res = M*dn_dlnM(M,intpar,err)*Ngal(model,M);
  forwardError(*err, __LINE__, 0.0);

  return res;
}

double av_frsat(cosmo_hm *model, double a, error **err)
{
  /*returns the fraction of satellite galaxies wrt redshift (depends on HOD).*/
 
  double ng,res;
  cosmo_hmANDstuff2 intpar;
  intpar.self = model;
  intpar.a    = a;
  intpar.asymptotic = 0;

  double min = logMmin;

  ng = ngal_den(model,a,logMmax,err);
  forwardError(*err, __LINE__, 0.0);
  res = int_gsl(int_for_av_frsat,(void*)&intpar,min,logMmax,1.e-3,err);
  forwardError(*err, __LINE__, 0.0);

  testErrorRet(ng == 0, ce_infnan, "Division by zero (ng)", *err, __LINE__, 0.0);
  
  return 1.0 - res/ng;
}
double int_for_av_frsat(double logM, void *intpar, error **err)
{
  double res;
  double M         = exp(logM);
  cosmo_hm *model  = ((cosmo_hmANDstuff2 *)intpar)->self; 
  
  res = dn_dlnM(M,intpar,err)*Ngal_c(model,M);;
  forwardError(*err, __LINE__, 0.0);

  return res;
}


double av_halo_bias(cosmo_hm *model, double a, error **err)
{
  /*returns the average halo bias wrt to redshift (independent of HOD).
   THIS IS NOT TESTED*/
  double res1,res2;
  cosmo_hmANDstuff2 intpar;
  intpar.self = model;
  intpar.a    = a;
  intpar.asymptotic = 0;

  //Jean: to check
  //double min = log10(model->M_min);
  double min = logMmin;

  res1 = int_gsl(int_for_av_halo_bias,(void*)&intpar,min,logMmax,1.e-3,err);
  forwardError(*err, __LINE__, 0.0);
  res2 = int_gsl(dn_dlnM_integrand,(void*)&intpar,min,logMmax,1.e-3,err);
  forwardError(*err, __LINE__, 0.0);

  testErrorRet(res2 == 0, ce_infnan, "Division by zero (res2)", *err, __LINE__, 0.0);

  return res1/res2;
}

double int_for_av_halo_bias(double logM, void *intpar, error **err)
{
  double b,res;
  double M = exp(logM);

  cosmo_hm *model  = ((cosmo_hmANDstuff2 *)intpar)->self; 
  double        a  = ((cosmo_hmANDstuff2 *)intpar)->a; 
  
  b   = halo_bias(model, M, a, 1, err);
  res = dn_dlnM(M,intpar,err)*b;
  forwardError(*err, __LINE__, 0.0);

  return res;
}

double Ngal_mean(cosmo_hm *model, double a, error **err)
{
  /*returns the mean number of galaxies per halo. NOT TESTED.*/
  double res;
  cosmo_hmANDstuff2 intpar;
  intpar.self = model;
  intpar.a    = a;
  intpar.asymptotic = 0;
  
  double ng = ngal_den(model,a,logMmax,err);
  forwardError(*err, __LINE__, 0.0);
  
  res = ng/int_gsl(dn_dlnM_integrand,(void*)&intpar,log(model->M_min),logMmax,1.e-3,err);
  forwardError(*err, __LINE__, 0.0);
  return res;
}

//Used by everyone
double dn_dlnM_integrand(double logM, void *intpar, error **err)
{
  double res;
  double M = exp(logM);

  res = dn_dlnM(M,intpar,err);
  forwardError(*err, __LINE__, 0.0);
  
  return res;
}

/*----------------------------------------------------------------*
 *Physical quantities                                             *
 *----------------------------------------------------------------*/

double vc(cosmo_hm *model, double zmin, double zmax, error **err)
{
  /*Comoving volume per unit solid angle between zmin and zmax*/
  
  double res;
  cosmo_hm_params params;
  params.model = model;
  
  if (zmin > -0.5 && zmax > -0.5) {
     res = int_gsl(int_for_vc, (void*)&params, zmin, zmax, 1e-4, err);
     forwardError(*err, __LINE__, 0.0);
  } else {
     Nz_hist *nz = get_nz(model->redshift,err);
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

double M_vir(cosmo_hm *model, double r, double a)
{
  /*halo mass M corresponding to a virial radius r*/
  double Delta = Delta_vir(model,a);//critical overdensity for virialization
  return (4.0/3.0)*pi*r*r*r*rho_c0*model->cosmo->Omega_m*Delta;
}

double R_vir(cosmo_hm *model, double M, double a)
{
  /*virial radius R corresponding to a halo mass M*/
  double Delta = Delta_vir(model,a);//critical overdensity for virialization
  double result = (3.0/4.0)*M/(pi*rho_c0*model->cosmo->Omega_m*Delta);
  return  pow(result,1.0/3.0);
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
    ng = ngal_den(model,a,logM, err);
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
  
  cosmo_hmANDstuff2 intpar;
  intpar.self       = model;
  intpar.a          = a;
  intpar.asymptotic = 0;
  
  double dlogM = (logMmax - logMmin)/(double)N;
  
  sum1 = 0.0;
  for(i=0;i<N;i++){
    logM1 = logMmin + dlogM*(double)i;
    R1 = R_vir(model,exp(logM1),a);
    sum2 = 0.0;
    for(j=0;j<N;j++){
      logM2 = logMmin + dlogM*(double)j;
      R2 = R_vir(model,exp(logM2),a);
      x = r/(R1 + R2);
      y = (x - 0.8)/0.29;
      if (y<0) {
	sum2 += 0.0;
      } else if (y>1) {
	sum2 +=  int_for_ngal_den(logM2,(void*)&intpar,err);
	forwardError(*err, __LINE__, 0.0);
      } else {
	P  = (3.0 - 2.0*y)*y*y;
	sum2 +=  int_for_ngal_den(logM2,(void*)&intpar,err)*P;
	forwardError(*err, __LINE__, 0.0);
      }
    }
    sum1 +=  int_for_ngal_den(logM1,(void*)&intpar,err)*sum2*dlogM;
    forwardError(*err, __LINE__, 0.0);
  }
  
  return sqrt(sum1*dlogM);
}


double ngal_den(cosmo_hm *model, double a, double logMlim, error **err)
{
  /*galaxy number density*/
  
  double res;

  cosmo_hmANDstuff2 intpar;
  intpar.self       = model;
  intpar.a          = a;
  intpar.asymptotic = 0;
  intpar.err        = err;

  res = int_gsl(int_for_ngal_den, (void*)&intpar, logMmin, logMlim, 1.e-5, err);
  forwardError(*err, __LINE__, 0.0);

  return res;
}

double int_for_ngal_den(double logM, void *intpar,  error **err) {
  cosmo_hm *model;
  double res, M;
  
  model = ((cosmo_hmANDstuff2 *)intpar)->self;
  M     = exp(logM);

  res   = Ngal(model, M)*dn_dlnM(M, intpar, err);
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
  
  double ans1,ans2;
  
  ans1 = ngal_weighted(model,err); forwardError(*err, __LINE__, 0.0);
  ans2 = vc(model,-1.0,-1.0,err);  forwardError(*err, __LINE__, 0.0);
  
  return ans1/ans2;
}

double ngal_weighted(cosmo_hm *model, error **err){
  /*Ngal times comoving volume per unit solid angle*/
  
  Nz_hist *nz  = get_nz(model->redshift,err);
  forwardError(*err, __LINE__, 0.0);
  int i;
  double sum = 0.0, tmp;
  
  for(i=0;i<nz->nbins;i++){ //loop over z
     tmp  = ngal_den(model,1.0/(1.0+nz->z[i]),logMmax,err)*nz->n[i];
     forwardError(*err, __LINE__, 0.0);
     tmp *= dvdz(model->cosmo,1.0/(1.0+nz->z[i]),err)/nz->max;
     forwardError(*err, __LINE__, 0.0);
     sum += tmp;
     //sum += nz->n[i]*dvdz(model->cosmo,1.0/(1.0+nz->z[i]),err);
     //forwardError(*err, __LINE__, 0.0);
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

  //interpolation
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


void FFTLog(FFTLog_config *fc, const gsl_function *ar_in,double *k, double *ak_out, int dir, error **err){
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
  
  double logrmin = log(fc->min);
  double logrmax = log(fc->max);
  double r,dlogr = (logrmax - logrmin)/(double)fc->N;
  double logrc   = (logrmax+logrmin)/2.0;
  double nc      = (double)(fc->N+1)/2.0-1;
  double logkc   = log(fc->kr)-logrc;
  
  //write signal
  for(i=0; i<fc->N; i++){
    k[i] = exp(logkc+((double)i-nc)*dlogr);
    r  = exp(logrc+((double)i-nc)*dlogr);
    fc->an[i][0] = ar_in->function(r,(void*)(ar_in->params))*pow(r,-(double)dir*fc->q);
    fc->an[i][1] = 0.0;
  }
  
  //cm's: FFT forward
  fftw_execute(fc->p_forward);
  
  //um*cm
  fc->cmum[0][0] = fc->cm[0][0]*fc->um[0][0] - fc->cm[0][1]*fc->um[0][1];
  fc->cmum[0][1] = fc->cm[0][0]*fc->um[0][1] + fc->cm[0][1]*fc->um[0][0];
  for(i=1;i<fc->N/2+1;i++){
    fc->cmum[i][0] = fc->cm[i][0]*fc->um[i][0] - fc->cm[i][1]*fc->um[i][1];
    fc->cmum[i][1] = fc->cm[i][0]*fc->um[i][1] + fc->cm[i][1]*fc->um[i][0];
    //Hermitian symetry (i.e. to get a real signal after FFT back)
    fc->cmum[fc->N-i][0] = fc->cmum[i][0];
    fc->cmum[fc->N-i][1] = -fc->cmum[i][1];
  }
  
  //ak's: FFT backward
  fftw_execute(fc->p_backward);
  
  //reverse array
  for(i=0;i<fc->N/2;i++) FFTLog_SWAP(fc->ak[i][0],fc->ak[fc->N-i-1][0]);
  
  //and write ak(k)
  for(i=0;i<fc->N;i++) ak_out[i] = fc->ak[i][0]*pow(k[i],-(double)dir*fc->q)/(double)fc->N;
  
  return;
}

FFTLog_config *FFTLog_init(int N, double min, double max, double q, double mu){
  /*Initializes what FFTLog needs.*/
  
  FFTLog_config *fc = (FFTLog_config*)malloc(sizeof(FFTLog_config));

  //FFTW3 Initialization
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
  
  //um's
  FFTLog_complex z,result;
  double L = log(max)-log(min);
  fc->kr   = 1.0;
  int i;
  
  for(i=0;i<fc->N/2+1;i++){
    z.re   = fc->q;
    z.im   = 2.0*pi*(double)i/L;
    result = FFTLog_U_mu(mu,z);
    
    //Multiply by (kr)^-2PIim/L
    result.amp *= 1.0;
    result.arg += -2.0*pi*(double)i*log(fc->kr)/L;
    
    fc->um[i][0] = result.amp*cos(result.arg);
    fc->um[i][1] = result.amp*sin(result.arg);
  }
  
  //If N even, mutiply by real part only
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
  
  //2^z
  amp1 = exp(z.re*log(2.0));
  arg1 = z.im*log(2.0);
  
  //Gamma 1
  FFTLog_complex zplus;
  zplus.re = (mu + 1.0 + z.re)/2.0;
  zplus.im = z.im/2.0;
  gsl_sf_lngamma_complex_e(zplus.re,zplus.im,&lnamp2,&arg2);
  
  //Gamma 2
  FFTLog_complex zminus;
  zminus.re = (mu + 1.0 - z.re)/2.0;
  zminus.im = - z.im/2.0;
  gsl_sf_lngamma_complex_e(zminus.re,zminus.im,&lnamp3,&arg3);

  //Result
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

Nz_hist *get_nz(redshift_t *nz_mk,error **err){
  /*
    Put Martin Kilbinger's n(z) format into Nz_hist
    format.

    0  1  | 2  3  ... n       | n+1 n+2 ... 2n
    -----------------------------------------------
    z0 zn | z1 z2 ... z_{n-1} | N0  N1  ... N_{n-1}
    nbins = n
  */
  
  int i;
  Nz_hist *result = malloc_err(sizeof(Nz_hist), err);                forwardError(*err, __LINE__, NULL);
 
  result->nbins   = (nz_mk->Nnz[0]-1)/2;
  result->z       = malloc_err(result->nbins*sizeof(double), err);   forwardError(*err, __LINE__, NULL);
  result->n       = malloc_err(result->nbins*sizeof(double), err);   forwardError(*err, __LINE__, NULL);
  
  result->z[0]  = nz_mk->par_nz[0];
  result->n[0]  = nz_mk->par_nz[result->nbins+1];
  result->max   = result->n[0];

  for(i=1;i<result->nbins;i++){
    result->z[i] = nz_mk->par_nz[i+1];
    result->n[i] = nz_mk->par_nz[result->nbins+i+1];
    result->max = MAX(result->max,result->n[i]);
  }
  
  //normalize AND return norm
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

double int_gsl(funcwithpars func, void *params, double a, double b, double eps, error **err)
{
  int n = 1000;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (n);
  double result, result_err;
  
  gsl_function F;
  F.function = &integrand_gsl;
  
  gsl_int_params p;
  p.func   = func;
  p.err    = err;
  p.params = params;
  F.params = &p;
  
  gsl_integration_qag (&F,a,b,eps,eps,n,GSL_INTEG_GAUSS51,w,&result,&result_err);
  forwardError(*err, __LINE__, 0.0);
  gsl_integration_workspace_free (w);
  
  return result;
}

double integrand_gsl(double x,void *p)
{
  double res;
  error **err  =  ((gsl_int_params *)p)->err;
  void *params =  ((gsl_int_params *)p)->params;
  
  res = ((gsl_int_params *)p)->func(x,params,err);
  forwardError(*err, __LINE__, 0.0);
  return res;
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
   wtr->RR = malloc_err(sizeof(double)*nlines, err);       forwardError(*err, __LINE__, NULL);

   return wtr;
}

void free_wt_t(wt_t **self)
{
   wt_t *s;

   s = *self;
   free(s->th);
   free(s->w);
   if (s->werr) free(s->werr);
   if (s->wcov) free(s->wcov);
   free(s->RR);
   free(s);
   self = NULL;
}

wt_t *read_wtheta(const char *name, double delta, double intconst, error **err)
{
  /* read in a hjmcc-format w-theta file into an array. */
  
  size_t i, Ncol, res;
  char line[NFIELD*NCHAR], item[NFIELD*NCHAR];
  
  FILE *fileIn = fopen_err(name, "r", err); forwardError(*err, __LINE__, NULL);
  
  //Count the number of lines
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
      /* OLD FORMAT
      wtr->th[i]   = pow(10.0,getDoubleValue(item,1));
      wtr->w[i]    = getDoubleValue(item,8);
      wtr->werr[i] = getDoubleValue(item,10);
      */
      wtr->th[i]   = getDoubleValue(item,1);
      wtr->w[i]    = getDoubleValue(item,2);
      wtr->werr[i] = getDoubleValue(item,3);
      wtr->w[i]    = wtr->w[i]*pow(wtr->th[i],-delta)/
      	(pow(wtr->th[i],-delta)-intconst);
      i++;
    }
  }

  fclose(fileIn);
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


nz_t *read_dndz(char *infile,error **err)
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
   CONFIG_READ(config, area_deg, d, F, c, err);
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

