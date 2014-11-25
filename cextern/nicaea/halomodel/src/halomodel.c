/* ============================================================ *
 * halomodel.c							*
 * Martin Kilbinger, Jean Coupon 2006-2012                      *
 * Refs:							*
 *   - Takada & Jain 2003					*
 *   - Eisenstein & Hu 1998					*
 *   - Percival 2005						*
 *   - Weinberg & Kamionkowski 2003				*
 *   - Kitayama & Suto 1996					*
 *   - Cooray & Hu 2001						*
 *   - Problem set from W.Hu's homepage				*
 *   - Bartelmann & Schneider 2001				*
 * ============================================================ */

#include "halomodel.h"

cosmo_hm* init_parameters_hm(double OMEGAM, double OMEGADE, double W0_DE, double W1_DE,
			     double *W_POLY_DE, int N_POLY_DE,
			     double H100, double OMEGAB, double OMEGANUMASS, 
			     double NEFFNUMASS, double NORM, double NSPEC,
			     int Nzbin, const int *Nnz, const nofz_t *nofz, double *par_nz,
			     double zmin, double zmax,
			     nonlinear_t NONLINEAR, transfer_t TRANSFER, growth_t GROWTH,
			     de_param_t DEPARAM, norm_t normmode,
			     double C0, double ALPHANFW, double BETANFW, massfct_t MASSFCT, halo_bias_t HALO_BIAS,
			     double M_min, double M1, double M0, double sigma_log_M, double alpha,
			     double Mstar0, double beta, double delta, double gamma, double B_cut, double B_sat, 
			     double beta_cut, double beta_sat, double Mstellar_min, double Mstellar_max, double eta,
			     double fcen1, double fcen2,
			     hod_t HOD, double pi_max, error **err)
{
   cosmo_hm *res;
   double amin;

   res = (cosmo_hm*)malloc_err(sizeof(cosmo_hm), err);
   forwardError(*err, __LINE__, NULL);

   res->redshift = init_redshift(Nzbin, Nnz, nofz, par_nz, NULL, err);
   forwardError(*err, __LINE__, NULL);

   res->zmin = zmin;
   res->zmax = zmax;

   amin = get_amin(res->redshift, err);
   forwardError(*err, __LINE__, NULL);

   res->cosmo = init_parameters(OMEGAM, OMEGADE, W0_DE, W1_DE, W_POLY_DE, N_POLY_DE,
				H100, OMEGAB, OMEGANUMASS, NEFFNUMASS, NORM, NSPEC,
				NONLINEAR, TRANSFER, GROWTH,
				DEPARAM, normmode, amin, err);
   forwardError(*err, __LINE__, NULL);

   /* Halomodel parameters */
   res->c0        = C0;
   res->alpha_NFW = ALPHANFW;
   res->beta_NFW  = BETANFW;
   res->massfct   = MASSFCT;
   set_massfct(res->massfct, &(res->nmz_a), &(res->nmz_p), err);
   res->halo_bias = HALO_BIAS;

   testErrorRetVA(M_min>pow(logMmax,10.0), hm_Mmin, "M_min=%d larger than maximum mass (%d)",
		  *err, __LINE__, NULL, M_min, pow(logMmax,10.0));

   res->hod         = HOD;

   /* All models */
   res->M1          = M1;
   res->sigma_log_M = sigma_log_M;
   res->alpha       = alpha;
   
   /* hamana, berwein models */
   res->M_min       = M_min;
   res->M0          = M0;
   res->M_min       = M_min;
   res->eta         = eta;

   /* leauthaud11 model */
   res->Mstar0      = Mstar0;
   res->beta        = beta;  
   res->delta       = delta;
   res->gamma       = gamma;
   res->B_cut       = B_cut; 
   res->B_sat       = B_sat; 
   res->beta_cut    = beta_cut; 
   res->beta_sat    = beta_sat;
   res->fcen1        = fcen1;
   res->fcen2        = fcen2;
   res->Mstellar_min    = Mstellar_min;
   res->Mstellar_max    = Mstellar_max;

   testErrorRet(res->hod == leauthaud11 && res->Mstellar_min < 0, hm_undef,
                "Minimum stellar mass cannot be set to a negative value (undef) for hod=leauthaud",
                *err, __LINE__, NULL);
   testErrorRet(res->hod == leauthaud11 && res->Mstellar_max < 0, hm_undef,
                "Maximum stellar mass cannot be set to a negative value (undef) for hod=leauthaud",
                *err, __LINE__, NULL);
   testErrorRet(res->Mstellar_min > res->Mstellar_max && res->Mstellar_max > 0, hm_undef,
                "Maximum has to be larger than minimum stellar mass",
                *err, __LINE__, NULL);
   
   res->pi_max      = pi_max;

   

   /* Reset precomputed values */
   res->A          = -1.0;
   res->Mstar      = -1.0;
   res->rhohat     = NULL;
   res->sigRsqr    = NULL;
   res->Pthdm      = NULL;
   res->xir        = NULL;
   res->xi_dm      = NULL;


   res->a_xir      = -1.0;

   return res;
}

cosmo_hm* copy_parameters_hm_only(cosmo_hm* source, error **err)
{
   cosmo_hm* res;

   res = init_parameters_hm(source->cosmo->Omega_m, source->cosmo->Omega_de, source->cosmo->w0_de,
			    source->cosmo->w1_de, source->cosmo->w_poly_de, source->cosmo->N_poly_de,
			    source->cosmo->h_100, source->cosmo->Omega_b,
			    source->cosmo->Omega_nu_mass, source->cosmo->Neff_nu_mass,
			    source->cosmo->normalization, source->cosmo->n_spec,
			    source->redshift->Nzbin, source->redshift->Nnz, source->redshift->nofz,
			    source->redshift->par_nz,
			    source->zmin,source->zmax,
			    source->cosmo->nonlinear, source->cosmo->transfer,
			    source->cosmo->growth, source->cosmo->de_param, 
			    source->cosmo->normmode, source->c0, source->alpha_NFW, source->beta_NFW,
			    source->massfct, source->halo_bias, source->M_min, source->M1, 
			    source->M0, source->sigma_log_M, source->alpha, 
			    source->Mstar0, source->beta, source->delta, source->gamma, source->B_cut, source->B_sat, 
			    source->beta_cut, source->beta_sat, source->Mstellar_min, source->Mstellar_max, source->eta,
			    source->fcen1, source->fcen2,
			    source->hod, source->pi_max, err);
   forwardError(*err, __LINE__, NULL);

   return res;
}

cosmo_hm *copy_parameters_hm(cosmo_hm *source, error **err)

/* ********************************
 * ********************************
 * THIS FUNCTION IS NOT UP TO DATE
 * ********************************
 * *********************************/
  
{
   cosmo_hm *res;
   int Nzcorr;

   Nzcorr = source->redshift->Nzbin*(source->redshift->Nzbin+1)/2;

   

   res = copy_parameters_hm_only(source, err);
   forwardError(*err, __LINE__, NULL);



   /* Reset cosmo and redshift before copying parameters and tables */
   free_parameters(&res->cosmo);
   res->cosmo = copy_parameters(source->cosmo, err);         forwardError(*err, __LINE__, NULL);


   res->zmin      = source->zmin;
   res->zmax      = source->zmax;

   free_redshift(&res->redshift);
   res->redshift  = copy_redshift(source->redshift, err);    forwardError(*err, __LINE__, NULL);

   res->Mstar      = source->Mstar;
   res->A          = source->A;
   res->Pthdm      = copy_interTable2D(source->Pthdm, err);     forwardError(*err, __LINE__, NULL);
   res->xir        = copy_interTable(source->xir, err);         forwardError(*err, __LINE__, NULL);
   res->xi_dm      = copy_interTable(source->xi_dm, err);       forwardError(*err, __LINE__, NULL);
   res->rhohat     = copy_interTable2D(source->rhohat, err);    forwardError(*err, __LINE__, NULL);
   res->sigRsqr    = copy_splineTable(source->sigRsqr, err);    forwardError(*err, __LINE__, NULL);
   res->a_xir      = source->a_xir;

   return res;
}

void read_cosmological_parameters_hm(cosmo_hm **model, FILE *F, error **err)
{
   cosmo_hm *tmp;
   config_element c = {0, 0.0, ""};
   struct { char cosmo_file[CSLENS], nofz_file[CSLENS], shod[CSLENS], smassfct[CSLENS], shalo_bias[CSLENS]; } tmp2;
   int j;
   FILE *FD;

   tmp = set_cosmological_parameters_to_default_hm(err);
   forwardError(*err, __LINE__,);

   /* Cosmology */
   CONFIG_READ_S(&tmp2, cosmo_file, s, F, c, err);
   if (strcmp(tmp2.cosmo_file, "-")!=0) {
      FD = fopen_err(tmp2.cosmo_file, "r", err);
      forwardError(*err, __LINE__,);
   } else {
      FD = F;
   }
   read_cosmological_parameters(&(tmp->cosmo), FD, err);
   forwardError(*err, __LINE__,);
   if (strcmp(tmp2.cosmo_file, "-")!=0) fclose(FD);

   /* Redshift distribution */
   CONFIG_READ_S(&tmp2, nofz_file, s, F, c, err);
   if (strcmp(tmp2.nofz_file, "-")!=0) {
      FD = fopen_err(tmp2.nofz_file, "r", err);
      forwardError(*err, __LINE__,);
   } else {
      FD = F;
   }
   read_redshift_info(&(tmp->redshift), FD, err);
   forwardError(*err, __LINE__,);
   
   /* not used anymore 
   CONFIG_READ(tmp,zmin, d, F, c, err);
   CONFIG_READ(tmp,zmax, d, F, c, err);
   */

   /* Halomodel parameters (dark matter) */
   CONFIG_READ(tmp, alpha_NFW, d, F, c, err);
   CONFIG_READ(tmp, c0, d, F, c, err);
   CONFIG_READ(tmp, beta_NFW, d, F, c, err);
   CONFIG_READ_S(&tmp2, smassfct, s, F, c, err);
   STRING2ENUM(tmp->massfct, tmp2.smassfct, massfct_t, smassfct_t, j, Nmassfct_t, err);

   CONFIG_READ_S(&tmp2, shalo_bias, s, F, c, err);
   STRING2ENUM(tmp->halo_bias, tmp2.shalo_bias, halo_bias_t, shalo_bias_t, j, Nhalo_bias_t, err);

   /* for wp(rp) */
   CONFIG_READ(tmp, pi_max, d, F, c, err);
   
   /* HOD model */
   CONFIG_READ_S(&tmp2, shod, s, F, c, err);
   STRING2ENUM(tmp->hod, tmp2.shod, hod_t, shod_t, j, Nhod_t, err);
   testErrorRetVA(tmp->hod!=hod_none && tmp->hod!=hamana04 && tmp->hod!=berwein02 && tmp->hod!=berwein02_hexcl && tmp->hod!=leauthaud11,
		  hm_hodtype, "HOD type (%d) unknown", *err, __LINE__,, tmp->hod, hamana04);
   /* sample properties */
   CONFIG_READ(tmp, Mstellar_min, d, F, c, err);
   CONFIG_READ(tmp, Mstellar_max, d, F, c, err);

   /* HOD parameters */
   switch (tmp->hod) {
      case berwein02 : case berwein02_hexcl: case hamana04 :

         CONFIG_READ(tmp, M_min, d, F, c, err);
         CONFIG_READ(tmp, M1, d, F, c, err);
         CONFIG_READ(tmp, M0, d, F, c, err);
         CONFIG_READ(tmp, sigma_log_M, d, F, c, err);
         CONFIG_READ(tmp, alpha, d, F, c, err);
         CONFIG_READ(tmp, eta, d, F, c, err);
         break;

      case leauthaud11:

         //testErrorRet(Mstellar_min < 0 || Mstellar_max < 0, hm_undef, "Mstellar_min and Mstellar_max have to be defined (>0)",
         //*err, __LINE__,);

         CONFIG_READ(tmp, M1, d, F, c, err);
         CONFIG_READ(tmp, Mstar0, d, F, c, err);
         CONFIG_READ(tmp, beta, d, F, c, err);
         CONFIG_READ(tmp, delta, d, F, c, err);
         CONFIG_READ(tmp, gamma,d, F, c, err);
         CONFIG_READ(tmp, sigma_log_M, d, F, c, err);
         CONFIG_READ(tmp, B_cut, d, F, c, err);
         CONFIG_READ(tmp, B_sat, d, F, c, err);
         CONFIG_READ(tmp, beta_cut, d, F, c, err);
         CONFIG_READ(tmp, beta_sat, d, F, c, err);
         CONFIG_READ(tmp, alpha, d, F, c, err);
         CONFIG_READ(tmp, fcen1, d, F, c, err);
         CONFIG_READ(tmp, fcen2, d, F, c, err);
         break;

      default:
         break;

   }

   *model = copy_parameters_hm_only(tmp, err);
   forwardError(*err, __LINE__,);
}

#define NZBIN 1
#define NNZ 5
cosmo_hm *set_cosmological_parameters_to_default_hm(error **err)
{
   int    Nnz[NZBIN]        = {NNZ};
   double par_nz[NZBIN*NNZ] = {0.0, 6.0, 0.612, 8.125, 0.62};
   nofz_t nofz[NZBIN]       = {ymmk};

   return init_parameters_hm(0.25, 0.75, -1.0, 0.0, NULL, 0, 0.70, 0.044, 0.0, 0.0, 0.80, 1.0,
         NZBIN, Nnz, nofz, par_nz,
         0.0, 10.0,			     
         smith03, eisenhu, growth_de, linder, norm_s8,
         9.0, 1.5, 0.13, st2, halo_bias_sc, 1.6e11, 8.0e12, 1.6e11, 0.3, 0.75, 
         1.0e11, 0.5, 0.6, 1.5, 1.5, 10.62, -0.13, 0.9, 1.0e10, -1, 1.0, 0.15, 0.5, hamana04, 60.0, err);
}
#undef NZBIN
#undef NNZ

void free_parameters_hm(cosmo_hm** model)
{
   cosmo_hm *s;
   int Nzcorr;

   s = *model;

   Nzcorr = s->redshift->Nzbin*(s->redshift->Nzbin+1)/2;

   del_interTable2D(&s->Pthdm);
   del_interTable(&s->xir);
   del_interTable(&s->xi_dm);
   del_interTable2D(&s->rhohat);
   del_splineTable(&s->sigRsqr);

   free_parameters(&s->cosmo);
   free_redshift(&s->redshift);

   free(s);
   s = NULL;
}

/* Sets (p,a=q), the parameters of the mass function */
void set_massfct(massfct_t massfct, double *nmz_a, double *nmz_p, error **err)
{
   switch (massfct) {
      case ps  : *nmz_p = 0.0; *nmz_a = 1.0;           break;  /* Press-Schechter    */
      case st  : *nmz_p = 0.3; *nmz_a = 0.75;          break;  /* Sheth-Tormen       */
      case st2 : *nmz_p = 0.3; *nmz_a = 1.0/sqrt(2.0); break;  /* aussi Sheth-Tormen */
      case j01 : *nmz_p = 0.3; *nmz_a = 1.0/sqrt(2.0); break;  /* Jenkins            */
      default  : *err = addError(ce_unknown, "Wrong or not supported mass function", *err, __LINE__);
	         return;
   }
}

void dump_param_only_hm(cosmo_hm* model, FILE *F)
{

  if (!F) F = stderr;
  fprintf(F, "#     c0  a_NFW  b_NFW  mnz_a  mnz_p mass_func   bias             model\n");
  fprintf(F, "# %6.2f % .3f % .3f % .3f % .3f %s(%d) %s(%d) %s(%d)\n",
	  model->c0, model->alpha_NFW, model->beta_NFW,
	  model->nmz_a, model->nmz_p,
	  smassfct_t(model->massfct), model->massfct, shalo_bias_t(model->halo_bias), model->halo_bias,
	  shod_t(model->hod), model->hod);
  if (model->hod == berwein02_hexcl || model->hod == berwein02){
    fprintf(F, "#   M_min  M1    M0    slogM alpha m b h\n");
    fprintf(F, "# %6.2e %6.2e %6.2e %.3f %.3f\n",
	    model->M_min, model->M1, model->M0, model->sigma_log_M, model->alpha);
  }else if (model->hod == leauthaud11){
    fprintf(F, "# M1       Mstar0   beta  delta gamma sigma_log_M B_cut B_sat  beta_cut beta_sat \n");
    fprintf(F, "# %6.2e %6.2e %.3f %.3f %.3f %.3f        %.3f %.3f %.3f   %.3f\n",
	    model->M1, model->Mstar0, model->beta, model->delta, model->gamma, model->sigma_log_M, 
	    model->B_cut, model->B_sat, model->beta_cut, model->beta_sat); 
    fprintf(F, "# alpha fcen1  fcen2   Mstellar_min Mstellar_max \n");
    fprintf(F, "# %.3f %.3f %.3f   %6.2e     %6.2e\n",
	    model->alpha, model->fcen1, model->fcen2, model->Mstellar_min, model->Mstellar_max); 
  }
  
  
}

void dump_param_hm(cosmo_hm* model, FILE *F, error **err)
{
   if (!F) F = stderr;
   dump_param(model->cosmo, F);
   dump_redshift(model->redshift, F, err);
   forwardError(*err, __LINE__,);
   dump_param_only_hm(model, F);
}

/* Could go into smith2 */
#define JMAX 40
double sm2_rtbis(double (*func)(double, void *, error **), double x1, double x2,
		 double xacc, void *param, error **err)
{
   int j;
   double dx,f,fmid,xmid,rtb;

   f=(*func)(x1,param,err);         forwardError(*err, __LINE__, 0);
   fmid=(*func)(x2,param,err);      forwardError(*err, __LINE__, 0);
   testErrorRet(f*fmid>=0.0, ce_overflow, "Root must be bracketed for bisection", *err, __LINE__, -1);

   rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
   for (j=1;j<=JMAX;j++) {
      fmid=(*func)(xmid=rtb+(dx *= 0.5),param,err);
      if (fmid <= 0.0) rtb=xmid;
      if (fabs(dx) < xacc || fmid == 0.0) return rtb;
   }

   *err = addError(ce_tooManySteps, "Too many bisections steps", *err, __LINE);
   return -1;

}

/* Cosine and Sine integral */

dcomplex Complex(double re, double im)
{
	dcomplex c;
	c.r=re;
	c.i=im;
	return c;
}

dcomplex Cadd(dcomplex a, dcomplex b)
{
	dcomplex c;
	c.r=a.r+b.r;
	c.i=a.i+b.i;
	return c;
}

dcomplex Cmul(dcomplex a, dcomplex b)
{
	dcomplex c;
	c.r=a.r*b.r-a.i*b.i;
	c.i=a.i*b.r+a.r*b.i;
	return c;
}

dcomplex Cdiv(dcomplex a, dcomplex b)
{
	dcomplex c;
	double r,den;
	if (fabs(b.r) >= fabs(b.i)) {
		r=b.i/b.r;
		den=b.r+r*b.i;
		c.r=(a.r+r*a.i)/den;
		c.i=(a.i-r*a.r)/den;
	} else {
		r=b.r/b.i;
		den=b.i+r*b.r;
		c.r=(a.r*r+a.i)/den;
		c.i=(a.i*r-a.r)/den;
	}
	return c;
}

dcomplex RCmul(double x, dcomplex a)
{
	dcomplex c;
	c.r=x*a.r;
	c.i=x*a.i;
	return c;
}

#undef JMAX
#define EPS 6.0e-8
#define EULER 0.57721566
#define MAXIT 100
#define PIBY2 1.5707963
#define FPMIN 1.0e-30
#define TMIN 2.0
#define TRUE 1
#define ONE Complex(1.0,0.0)
void sm2_cisi(double x, double *ci, double *si, error **err)
{
	int i,k,odd;
	double a,fact,sign,sum,sumc,sums,t,term, error;
	dcomplex h,b,c,d,del;

	t=fabs(x);
	if (t == 0.0) {
		*si=0.0;
		*ci = -1.0/FPMIN;
		return;
	}
	if (t > TMIN) {
		b=Complex(1.0,t);
		c=Complex(1.0/FPMIN,0.0);
		d=h=Cdiv(ONE,b);
		for (i=2;i<=MAXIT;i++) {
			a = -(i-1)*(i-1);
			b=Cadd(b,Complex(2.0,0.0));
			d=Cdiv(ONE,Cadd(RCmul(a,d),b));
			c=Cadd(b,Cdiv(Complex(a,0.0),c));
			del=Cmul(c,d);
			h=Cmul(h,del);
			if (fabs(del.r-1.0)+fabs(del.i) < EPS) break;
		}
		testErrorRet(i>MAXIT, ce_tooManySteps, "cf failed in cisi", *err, __LINE,);
		h=Cmul(Complex(cos(t),-sin(t)),h);
		*ci = -h.r;
		*si=PIBY2+h.i;
	} else {
		if (t < sqrt(FPMIN)) {
			sumc=0.0;
			sums=t;
		} else {
			sum=sums=sumc=0.0;
			sign=fact=1.0;
			odd=TRUE;
			for (k=1;k<=MAXIT;k++) {
				fact *= t/k;
				term=fact/k;
				sum += sign*term;
				error=term/fabs(sum);
				if (odd) {
					sign = -sign;
					sums=sum;
					sum=sumc;
				} else {
					sumc=sum;
					sum=sums;
				}
				if (error<EPS) break;
				odd=!odd;
			}
			testErrorRet(k>MAXIT, ce_tooManySteps, "cf failed in cisi", *err, __LINE,);
		}
		*si=sums;
		*ci=sumc+log(t)+EULER;
	}
	if (x < 0.0) *si = -(*si);
}
#undef EPS
#undef EULER
#undef MAXIT
#undef PIBY2
#undef FPMIN
#undef TMIN
#undef TRUE
#undef ONE

int change_delta_c(cosmo_hm *avant, cosmo_hm *apres)
{
   if (change_w_de(avant->cosmo, apres->cosmo)) return 1;

   return 0;
}

/* Critical collapse overdensity */
double delta_c(cosmo *model, double a, error **err)
{
   double delta_EdS, alpha, deltac;

   delta_EdS = 1.68647;

   /* Per05, D_+ ? */
   /* return delta_EdS; */

   /* WK03 (18) */
   alpha = 0.131 + model->w0_de*(0.555 + model->w0_de*(1.128 + 
           model->w0_de*(1.044 + model->w0_de*0.353)));

   /* KS96 (A.6). Note: typo (0.123) in astro-ph version. Correct in ApJ. */
   //alpha = 0.0123;

   deltac =  delta_EdS*(1. + alpha*log10(Omega_m_a(model, a, -1.0, err)));
   forwardError(*err, __LINE__, 0.0);

   return deltac;
}

/* Bisection for Mstar. */
double bis_Mstar(double logM, void *param, error **err)
{
   cosmo_hm *model;
   double diff, M;
 
   model = (cosmo_hm*)param;
   M     = exp(logM);

   diff  = delta_c(model->cosmo, 1.0, err);     forwardError(*err, __LINE__, 0);
   /* New: omitting D+(a=1, normalised=1) = 1 */
   diff /= sqrt(sigmasqr_M(model, M, err));     forwardError(*err, __LINE__, 0);
   
   return diff-1.0;
}

/* Bisection for Mstar_a. This redshift-dependent routine is not used in general. */
double bis_Mstar_a(double logM, void *param, error **err)
{
   cosmo_hm *model;
   cosmo_hm_params *cANDd;
   double diff, M, a;
 
   cANDd = (cosmo_hm_params*)param;
   model = cANDd->model;
   a    = cANDd->a;
   M     = exp(logM);

   diff  = delta_c(model->cosmo, a, err);     forwardError(*err, __LINE__, 0);
   diff /= D_plus(model->cosmo, a, 1, err);   forwardError(*err, __LINE__, 0);
   diff /= sqrt(sigmasqr_M(model, M, err));   forwardError(*err, __LINE__, 0);

   return diff-1.0;
}

/* dsigma^2(R)/dR, R in Mpc/h */
double dsigma_R_sqr_dR(cosmo_hm *model, double R, error **err)
{
   /* Numerical derivative */
   double res, h;
   h = R/20.0;
   /* Problematic ... ? */
   res = (sigma_R_sqr(model, R+h, err)-sigma_R_sqr(model, R-h, err))/(2.0*h);
   forwardError(*err, __LINE__, 0);

   if (-res < 0.0) {
      printf("dsigma_R_sqr_dR %g %g %g  %g %g %g\n", R, R+h, R-h, res, sigma_R_sqr(model, R+h, err), sigma_R_sqr(model, R-h, err));
   }

   testErrorRetVA(-res<0.0, ce_negative,  "-dsigma_R_sqr_dR=-%g negative", *err, __LINE__, 0.0, res);

   return res;
}

int change_massfct_params(cosmo_hm *avant, cosmo_hm *apres)
{
   if (NCOEQ(avant, apres, nmz_p) || NCOEQ(avant, apres, nmz_p)) return 1;
   return 0;
}

int change_massfct(cosmo_hm *avant, cosmo_hm *apres)
{
   if (NCOEQ(avant,apres,massfct) || NCOCLOSE(avant,apres,c0) || NCOCLOSE(avant,apres,alpha_NFW)
       || NCOCLOSE(avant,apres,beta_NFW)) return 1;
   return 0;
}

int change_halo_bias(cosmo_hm *avant, cosmo_hm *apres)
{
   if (NCOEQ(avant, apres, halo_bias)) return 1;
   return 0;
}

/* Mass function part */
double nufnu(cosmo_hm *model, double nu, int asymptotic, error **err)
{
   double res, qnusqr;
   
   if (model->A<0) {

      /* Sets int[f(nu) nu] = 1 */
      model->A = 1.0/(1.0 + pow(2.0, -model->nmz_p)*exp(gammln(0.5-model->nmz_p))/sqrt(pi));
      /* = 0.322 */

      /* matches linear power spectrum on large scales, CH01: with st2. *
       * int[f(nu) nu] about 19% off.					*/
      //model->A = 0.383;

   }

   testErrorRet(nu<0, ce_negative, "nu<0", *err, __LINE__, -1);

   if (asymptotic==0) {
      qnusqr = model->nmz_a*nu*nu;
      res = sqrt(2.0/pi*qnusqr)*(1.0 + pow(qnusqr, -model->nmz_p))*exp(-qnusqr/2.0);
   } else if (asymptotic==1) {
      res = pow(nu, 0.5-model->nmz_p);
   } else {
      *err = addError(ce_wrongValue, "The flag 'asymptotic' has to be 0 or 1", *err, __LINE__);
      return 0.0;
   }

   return res*model->A;
}

/* Mass function fo j01. *
 * x = D+ * sigma_m      */
double nufnu_j01(double x)
{
  return 0.315 * exp(-pow(fabs(log(1.0/(x)) + 0.61), 3.8));
  
  /* Tinker et al. (2008)
     TO DO: implement it as a user's choice
     double log_Delta = log10(250.0);
     double A = 0.1*log_Delta - 0.05;
     double a = 1.43 + pow(log_Delta - 2.30,  1.5);
     double b = 1.00 + pow(log_Delta - 1.60, -1.5);
     double c = 1.20 + pow(log_Delta - 2.35,  1.6);
     
     return A * (pow(x/b, -a)+1.0)*exp(-c/(x*x));
  */
  
}

int change_sigma_R_sqr(cosmo_hm *avant, cosmo_hm *apres)
{
   return change_Tsqr(avant->cosmo, apres->cosmo);
}

/* sigma(R)^2, power spectrum variance smoothed on a sphere with radius R, R in Mpc/h */
#define kmax  3.0e5
#define Rmin 1.0e-5
#define Rmax  700.0
#define Ns     2000
#define eps  1.0e-8
double sigma_R_sqr(cosmo_hm *model, double R, error **err)
{
   double h, a, b, RR, dlogR, res, norm;
   cosmoANDdouble cANDd;
   int i;

   if (model->sigRsqr==NULL) {

      model->sigRsqr = init_splineTable(Ns, err);
      forwardError(*err, __LINE__, 0);

      cANDd.self = model->cosmo;
      dlogR = (log(Rmax)-log(Rmin))/(Ns-1.0);
      norm  = 9.0*dsqr(model->cosmo->sigma_8)/sigma_8_sqr(model->cosmo, err);
      forwardError(*err, __LINE__, 0);
      for (i=1,RR=Rmin; i<=Ns; i++,RR*=exp(dlogR)) {
         cANDd.r = RR;
         model->sigRsqr->x[i] = log(RR);
         model->sigRsqr->y[i] = norm/(2.0*pi*pi)*sm2_qromberg(int_for_sigma_R, (void*)&cANDd, 
               log(k_min), log(kmax), eps, err);
         forwardError(*err, __LINE__, 0);
         model->sigRsqr->y[i] = log(model->sigRsqr->y[i]);
      }

      h = Rmin/20.0;
      cANDd.r = Rmin+h;
      a = norm/(2.0*pi*pi)*sm2_qromberg(int_for_sigma_R, (void*)&cANDd, log(k_min),
            log(kmax), eps, err);
      forwardError(*err, __LINE__, 0);
      cANDd.r = Rmin-h;
      b = norm/(2.0*pi*pi)*sm2_qromberg(int_for_sigma_R, (void*)&cANDd, log(k_min),
            log(kmax), eps, err);
      forwardError(*err, __LINE__, 0);
      model->sigRsqr->yp1 = (a-b)/(2.0*h)*Rmin/exp(model->sigRsqr->y[1]);
      model->sigRsqr->yp1 = 0.0;

      h = Rmax/20.0;
      cANDd.r = Rmax+h;
      a = norm/(2.0*pi*pi)*sm2_qromberg(int_for_sigma_R, (void*)&cANDd, log(k_min),
            log(kmax), eps, err);
      forwardError(*err, __LINE__, 0);
      cANDd.r = Rmax-h;
      b = norm/(2.0*pi*pi)*sm2_qromberg(int_for_sigma_R, (void*)&cANDd, log(k_min), 
            log(kmax), eps, err);
      forwardError(*err, __LINE__, 0);
      model->sigRsqr->ypn = (a-b)/(2.0*h)*Rmax/exp(model->sigRsqr->y[Ns]);

      sm2_spline(model->sigRsqr->y, model->sigRsqr->x, model->sigRsqr->n, model->sigRsqr->yp1, 
            model->sigRsqr->ypn, model->sigRsqr->y2, err);
      forwardError(*err, __LINE__, 0);
   }

   testErrorRetVA(R<Rmin||R>Rmax, ce_interpoloutofrange, "R = %g out of range [%g;%g]",
		  *err, __LINE__, 0.0, R, Rmin, Rmax);

   sm2_splint(model->sigRsqr->x, model->sigRsqr->y, model->sigRsqr->y2, model->sigRsqr->n, log(R), &res, err);
   forwardError(*err, __LINE__, 0.0);

   return exp(res);
}
#undef Rmin
#undef Rmax
#undef Ns
#undef eps

/* Calls sigma(R)^2, with M (in M_sol/h) mass in sphere of radius R */
double sigmasqr_M(cosmo_hm *model, double M, error **err)
{
   double R, res;


   // DEBUGGING
   double a = 1.0/(1.0+0.633);
   R = cbrt(3.0*M/(4.0*pi*Omega_m_halo(model, a, err)*rho_crit_halo(model, a, err)));
   res = sigma_R_sqr(model, R, err);
   forwardError(*err, __LINE__, 0);

   return res;
}

/* Returns dsigma^-1/dlnM. M in M_sol/h */
double dsigma_m1_dlnM(cosmo_hm *model, double M, error **err)
{
   double res, sigma, R, rhobar;
 // DEBUGGING
   double a = 1.0/(1.0+0.633);
   rhobar = Omega_m_halo(model, a, err)*rho_crit_halo(model, a, err);
   R      = cbrt(3.0*M/(4.0*pi*rhobar));
   sigma  = sqrt(sigma_R_sqr(model, R, err));             forwardError(*err, __LINE__, 0);
   res    = M/(2.0*sigma*sigma*sigma);
   res   *= dsigma_R_sqr_dR(model, R, err) * cbrt(1.0/(36.0*pi*rhobar*M*M));
   forwardError(*err, __LINE__, 0);

   testErrorRetVA(-res<0.0, ce_negative,  "-dsigma_m1_dlnM=-%g negative", *err, __LINE__, 0.0, res);

   return -res;
}

/* M in M_sol/h */
double dnu_dlnM(cosmo_hm *model, double M, double a, error **err)
{
  double res;
  
  res  = delta_c(model->cosmo, a, err);   forwardError(*err, __LINE__, 0.0);
  res /= D_plus(model->cosmo, a, 1, err); forwardError(*err, __LINE__, 0.0);
  res *= dsigma_m1_dlnM(model, M, err);   forwardError(*err, __LINE__, 0.0);
 
  testErrorRetVA(res<0.0, ce_negative,  "dnu_dlnM=%g negative", *err, __LINE__, 0.0, res);
  
  return res;
}

double dn_dlnM_lnM(double logM, void *params, error **err)
{
  double res, M;
  M = exp(logM);
  res = M*dn_dlnM(M, params, err);
  forwardError(*err, __LINE__, 0);
  return res;
}

/* Mass function, M in M_sol/h, user-friendly version */
double dn_dlnM_uf(double M, cosmo_hm *model, double a, error **err)
{
   cosmo_hm_params params;
   double res;

   params.model       = model;
   params.a          = a;
   params.asymptotic = 0;

   res = dn_dlnM(M, (void*)&params, err);
   forwardError(*err, __LINE__, -1.0);

   return res;
}

/* ========================================================================= *
 * Mass function: halo number density per log unit mass.                     *
 * dn / dlnM = dn / d ln sigma^{-1} * d ln sigma^{-1} / dlnM                 *
 *           = rho_0 / M * f(sigma) * sigma * dsigma^{-1} / dln M            *
 *           = rho_0 / M * f(sigma) * sigma * dsigma{^-1} / dnu * dnu / dlnM *
 *	     = rho_0 / M * f(sigma) * sigma * sigma^{-1} / nu * dnu / dlnM   *
 *           = rho_0 / M * f(sigma) / nu * dnu / dlnM                        *
 * The function f(sigma) is defined as                                       *
 *   f(sigma) = M / rho_0 dn / d ln sigma^{-1} (Jenkins et al. 2001, eq. 4)  *
 *            = nu f(nu)                       (Cooray & Sheth 2001)	     *
 * Units: [dn/dlnM] = h^3 Mpc^-3  ln(M_sol^-1 h).                            *
 * [M] = M_sol/h.                                                            *
 * Non-user-friendly version, called by user-friendly version dn_dlnM_uf.    *
 * ========================================================================= */
double dn_dlnM(double M, void *params, error **err)
{
  cosmo_hm_params *cANDs;
  cosmo_hm *model;
  double res = 0.0, nu, a, dp, sM, nfn, dnudlnM;
  int asym;   

  testErrorRetVA(M<0.0, ce_negative, "Mass M = %g, has to be positive", *err, __LINE__, 0.0, M);

  cANDs = (cosmo_hm_params*)params;
  model = cANDs->model;
  a     = cANDs->a;
  asym  = cANDs->asymptotic;

  if (model->massfct != ps && model->massfct != st && model->massfct != st2 && model->massfct != j01) {
     *err = addErrorVA(ce_wrongValue, "Unknown mass function type %d", *err, __LINE__, model->massfct);
     return 0.0;
  }
  
  dp = D_plus(model->cosmo, a, 1, err);          forwardError(*err, __LINE__, 0);
  sM = sqrt(sigmasqr_M(model, M, err));          forwardError(*err, __LINE__, 0);
  nu = delta_c(model->cosmo, a, err)/(dp*sM);    forwardError(*err, __LINE__, 0);
    
  testErrorRet(nu<0.0, ce_negative, "nu not positive",      *err, __LINE__, 0.0);
  testErrorRet(dp<0.0, ce_negative, "D+ not positive",      *err, __LINE__, 0.0);
  testErrorRet(sM<0.0, ce_negative, "sigma_M not positive", *err, __LINE__, 0.0);

  if (model->massfct == j01) {
     /* Jenkins 01: nu f(nu) = f(sigma) */
     /* MKDEBUG TODO: sigma or D+ * sigma ??? */
    nfn = nufnu_j01(dp * sM);
  } else {
    nfn = nufnu(model, nu, asym, err);  forwardError(*err, __LINE__, 0);
  }

  dnudlnM = dnu_dlnM(model, M, a, err);   forwardError(*err, __LINE__, 0);
  res     = Omega_m_halo(model, a, err)*rho_crit_halo(model, a, err)/M * nfn/nu * dnudlnM;

  testErrorRetVA(res<0.0, ce_negative,  "dn_dlnM=%g negative", *err, __LINE__, 0.0, res);
  testErrorRet(!finite(res), ce_infnan, "dn_dlnM inf or nan",  *err, __LINE__, 0.0);

  return res;
}



int change_Mstar(cosmo_hm *avant, cosmo_hm *apres)
{
   if (change_D_plus(avant->cosmo, apres->cosmo)) return 1;
   if (change_delta_c(avant, apres)) return 1;
   if (change_sigma_R_sqr(avant, apres)) return 1;
   return 0;
}

#define logMsmin 10.0
#define logMsmax 35.0
#define xacc 0.001

/* Returns M_* in M_sol/h */
double Mstar(cosmo_hm *model, error **err)
{
  if (model->Mstar<0) {
    model->Mstar = sm2_rtbis(bis_Mstar, logMsmin, logMsmax, xacc, (void*)model, err);
    forwardError(*err, __LINE__, -1);
    model->Mstar = exp(model->Mstar);
    /* model->Mstar = 1.0e14; */ 
  }
  
  return model->Mstar;
}

/* Returns M_* in M_sol/h, redshift-dependent! In general, M_star is used (z=0) */
double Mstar_a(cosmo_hm *model, double a, error **err)
{
   cosmo_hm_params cANDd;
   double Mstar;

   cANDd.model = model;
   cANDd.a     = a;

   Mstar = sm2_rtbis(bis_Mstar_a, logMsmin, logMsmax, xacc, (void*)(&cANDd), err);
   forwardError(*err, __LINE__, -1);

   return exp(Mstar);
}
#undef logMsmin
#undef logMsmax
#undef xacc

/* Concentration parameter, M in M_sol/h				 *
 * So far no de dependance!						 *
 * see Dolag et al 2004: alpha=0.1, c to be corrected by D+(de)/D+(LCDM) */
double concentration(cosmo_hm *model, double Mh, double a, error **err)
{
   double c, Ms;
   
   Ms = Mstar(model, err);
   forwardError(*err, __LINE__, 0);
   c = model->c0*a*pow(Mh/Ms, -model->beta_NFW);
   
   return c;
   
   /*
     This is Munoz-Cuartas et al. (2011),
     used in Leauthaud et al. (2011-12). 
     [Jean]: very small impact on w(theta)
   
     double z = 1.0/a - 1.0;
     double aa = 0.029*z - 0.097;
     double bb = -110.001/(z+16.885) + 2469.720/dsqr(z+16.885);
     double log10_c = aa*log10(Mh) + bb;
     return pow(10.0, log10_c);
   */
   
}


/* Virial overdensity, Delta_vir ~ 200 at z = 1 */
double Delta_vir(cosmo_hm *model, double a)
{
   double w_vir, D, om, ov;

   Omega_a(model->cosmo, a, &om, &ov);
   w_vir = 1.0/om - 1.0; // *pow(1+z, model->cosmo->w0_de * 3.0) ??

   /* KS96 (A.5) */
   //D = 18.0*pi*pi*(1.0 + 0.4093*pow(w_vir, 0.9052));

   /* WK03 Eqs 16 & 17 */
   double aa, b;
   aa = 0.399 - 1.309*(pow(fabs(model->cosmo->w0_de), 0.426) - 1.0);
   b = 0.941 - 0.205*(pow(fabs(model->cosmo->w0_de), 0.938) - 1.0);
   D = 18.0*pi*pi*(1.0 + aa * pow(w_vir, b));
   
   return D;
}


double M_vir(cosmo_hm *model, double r_vir, double a, error **err)
{
  /* Mass of a halo with radius r_vir. This is NOT Mh(r) */
  double Delta, rhocrit, Omega_m;

  Delta   = Delta_h(model, a, err);       forwardError(*err, __LINE__, 0.0);
  rhocrit = rho_crit_halo(model, a, err); forwardError(*err, __LINE__, 0.0);
  Omega_m = Omega_m_halo(model, a, err);  forwardError(*err, __LINE__, 0.0);
  
  return (4.0/3.0)*pi*r_vir*r_vir*r_vir*Delta*rhocrit*Omega_m;
}

double r_vir(cosmo_hm *model, double M, double a, error **err)
{
  /* Virial radius of a halo with mass M */
  double Delta, rhocrit, Omega_m;
  
  Delta   = Delta_h(model, a, err);       forwardError(*err, __LINE__, 0.0);
  rhocrit = rho_crit_halo(model, a, err); forwardError(*err, __LINE__, 0.0);
  Omega_m = Omega_m_halo(model, a, err);  forwardError(*err, __LINE__, 0.0);
  
  return cbrt((3.0/4.0)*M/(pi*Delta*rhocrit*Omega_m));  
}

double Delta_h(cosmo_hm *model, double a, error **err){
  /* Overdensity of dark matter haloes. Used to define r_s, 
   * the radius within which the mean matter density is 
   * rho_h = Delta_h X rho_bar
   * TO DO: allow user to set Delta_h and also in bias_tinker10 */
  
  return Delta_vir(model, a);
  //return 200.0;
  //return 200.0/Omega_m_a(model->cosmo, a, -1.0, err);
  //forwardError(*err, __LINE__, 0.0);
}

double rho_crit(cosmo_hm *model, double a, error **err){
  /* returns rho_crit at redshift z in [M_sol h^2 / Mpc^3]*/
  /* present critical density is  rho_c0 = 2.7754e11 */
  
  double G = 4.302e-9;                                                                 /* in  [km^2 s^-2 M_sun^-1 Mpc^1]  */
  double H = 100.0*sqrt(model->cosmo->Omega_m*pow(a, -3.0)+model->cosmo->Omega_de); /* H(z) in h km s^-1 Mpc^-1 */
  
  return 3.0*H*H/(8.0*pi*G);
}


/* ----------------------------------------------------------- *
 * WARNING: do not change the two functions below. 
 *  Those have only been use to check consistency of
 *  the 1-h term in DeltaSigma (where como_to_phys 
 *  should be changed to como_to_phys = 1.0).
 *  Everywhere in this code rho_crit_halo should be 
 *  =rho_c0 and Omega_m_halo=model->cosmo->Omega_m as
 *  all quantities are in comoving units.
 ----------------------------------------------------------- */
double rho_crit_halo(cosmo_hm *model, double a, error **err){
  /* Debugging stuff, do not touch */
  double res;
  
  res = rho_c0;
  /*
  res = rho_crit(model,  a,  err);
  forwardError(*err, __LINE__, 0.0);
  */
  return res;
}
double Omega_m_halo(cosmo_hm *model, double a, error **err){
  /* Debugging stuff, do not touch */
  double res;
  
  res = model->cosmo->Omega_m;
  /*
    res = Omega_m_a(model->cosmo, a, -1.0, err);
    forwardError(*err, __LINE__, 0.0);
  */
  
  return res;
}

#define xmin 1.0e-4
double rho_halo(cosmo_hm *model, double r, double a, double Mh, double c, error **err){
  /* Returns NFW profile in h^2 M_sol/Mpc^3. r in Mpc/h, M in M_sol/h.
   * Set c = -1 to have it computed in rho_halo. 
   */
  double result, rvir, Delta,  rhocrit, Omega_m;
  
  if(c < 0.0){
    /* if computed externally */
    c = concentration(model, Mh, a, err);
    forwardError(*err, __LINE__, 0.0);
  }
  
  rvir    = r_vir(model, Mh, a, err);     forwardError(*err, __LINE__, 0.0);
  Delta   = Delta_h(model, a, err);       forwardError(*err, __LINE__, 0.0);
  rhocrit = rho_crit_halo(model, a, err); forwardError(*err, __LINE__, 0.0);
  Omega_m = Omega_m_halo(model, a, err);  forwardError(*err, __LINE__, 0.0);
  
  /* If truncated halo. This matches Leauthaud et al. (2011) */
  /*
  if(r > rvir){
    return 0.0;
  }
  */
  
  double r_s   = rvir/c;
  double rho_s = rhocrit*Delta*Omega_m/3.0*c*c*c/(log(1.0+c)-c/(1.0+c)); /* rho_s = rho_crit * delta_char */
  forwardError(*err, __LINE__, 0.0);
  double x     = r/r_s;  
  
  if (x > xmin) {
    result = pow(x, -model->alpha_NFW)*pow(1.0+x, model->alpha_NFW-3.0);
  } else { 
    /* Approximation finite for r=0 for alpha_NFW<2 */
    result = pow(x, 2.0-model->alpha_NFW);
  }
  
  testErrorRet(!finite(result), ce_infnan, "inf or nan encountered", *err, __LINE__, 0);
  
  return rho_s*result;
}
#undef xmin


/* Returns Delta Sigma(r) in h M_sol/pc^2 (r in Mpc/h, M in M_sol/h).
 * See Wright & Brainerd (2000), Eqs 11-16.
 */

#define eps 1.0e-10
double DeltaSigma_WB2000(cosmo_hm *model, double r, const double a, const double Mh,  double c, double Delta, error **err){
  
  double result, rvir, Omega_m, rhocrit;
  
  if(c < 0.0){
    c = concentration(model, Mh, a, err);
    forwardError(*err, __LINE__, 0.0);
  }
  if(Delta < 0.0){
    Delta   = Delta_h(model, a, err);      forwardError(*err, __LINE__, 0.0);
    forwardError(*err, __LINE__, 0.0);
  }
  rvir    = r_vir(model, Mh, a, err);        forwardError(*err, __LINE__, 0.0);
  Omega_m = Omega_m_halo(model, a, err);  forwardError(*err, __LINE__, 0.0);
  rhocrit = rho_crit_halo(model, a, err); forwardError(*err, __LINE__, 0.0);
   
  double r_s   = rvir/c;
  double rho_s = rhocrit*Delta*Omega_m/3.0*c*c*c/(log(1.0+c)-c/(1.0+c)); /* rho_s = rho_crit * delta_char */
  forwardError(*err, __LINE__, 0.0);
  double x     = r/r_s;
  
  if(1.0 - eps < x && x < 1.0 + eps){ /* x = 1 */
    result = r_s*rho_s*(10.0/3.0 + 4.0*log(0.5));
  }else if(x < 1.0){                  /* x < 1 */
    result = r_s*rho_s*g_inf(x,err);
  }else{                              /* x > 1 */
    result = r_s*rho_s*g_sup(x,err);
  }                                  
  
  return result/1.0e12;
}

double g_inf(double x, error **err){
  double result;
  
  result  = 8.0*atanh(sqrt((1.0 - x)/(1.0 + x)))/(x*x*sqrt(1.0 - x*x));
  result += 4.0/(x*x) * log(x/2.0);
  result += -2.0/(x*x - 1.0);
  result += 4.0*atanh(sqrt((1.0 - x)/(1.0 + x)))/((x*x - 1.0)*pow(1.0 - x*x,0.5));
  
  return result;
}

double g_sup(double x, error **err){
  double result;
  
  result  = 8.0*atan(sqrt((x - 1.0)/(1.0+x)))/(x*x*sqrt(x*x - 1.0));
  result += 4.0/(x*x) * log(x/2.0);
  result += -2.0/(x*x - 1.0);
  result += 4.0*atan(sqrt((x - 1.0)/(1.0 + x)))/(pow(x*x - 1.0,1.5));
  
  return result;
}
#undef eps

int change_rhohat_halo(cosmo_hm *avant, cosmo_hm *apres)
{
   if (change_massfct(avant, apres)) return 1;
   if (change_Mstar(avant, apres)) return 1;
   if (change_Esqr(avant->cosmo, apres->cosmo)) return 1;
   if (change_w_de(avant->cosmo, apres->cosmo)) return 1;

   return 0;
}


/* ============================================================ *
 * Fourier Transform of halo profile.				*
 * k in h/Mpc, M in M_sol/h.					*
 * Integration of rho up to xvir*r_vir. Should be unity to be   *
 * consistent with closed NFW formula.				*
 * ============================================================ */


#define xvir       1.0
#define EPS        1.0e-6

double rhohat_halo(cosmo_hm *model, double k, double M, double a, double c, error **err)
{

  double res, rvir;
  double  f, eta, cieta, sieta, cieta1pc, sieta1pc;
  cosmo_hm_params params;
  
  rvir = r_vir(model, M, a, err);
  forwardError(*err, __LINE__, 0.0);
  if(c < 0.0){
    c = concentration(model, M, a, err);
    forwardError(*err, __LINE__, 0.0);
  }

  /* Closed formula for NFW profile (alpha=1) */
  if (fabs(model->alpha_NFW-1.0)<EPS) {
  
    f = 1.0/(log(1.0+c) - c/(1.0+c));
    eta = k*rvir/c;
    sm2_cisi(eta, &cieta, &sieta, err);                 forwardError(*err, __LINE__, 0);
    sm2_cisi(eta*(1.0+c), &cieta1pc, &sieta1pc, err);   forwardError(*err, __LINE__, 0);
    
    /* TJ03 (17) */
    res = f*(sin(eta)*(sieta1pc - sieta) + cos(eta)*(cieta1pc - cieta)
	     - sin(eta*c)/(eta*(1.0+c)));
    
  } else {
       
    double norm;
       
#define logrmin -6.0
       
    params.model  = model;
    params.k     = k;
    params.M     = M;
    params.a     = a;
    params.c     = c;
    
    res          = 0.0;
    
    params.logintegrate = +1;
    res += 4.0*pi/M*int_gsl(int_for_rhohat, (void*)&params, logrmin, log(xvir*rvir), EPS, err);
    forwardError(*err, __LINE__, 0);
       
    params.logintegrate = -1;
    res += 4.0*pi/M*int_gsl(int_for_rhohat, (void*)&params, 0.0, exp(logrmin), EPS, err);
    forwardError(*err, __LINE__, 0);
       
    /* Normalization -> rhohat(k=0) = 1 */
    params.k = 0.0;
       
    params.logintegrate = +1;
    norm = 4.0*pi/M*int_gsl(int_for_rhohat, (void*)&params, logrmin, log(xvir*rvir), EPS, err);
    forwardError(*err, __LINE__, 0);
    params.logintegrate = -1;
    norm += 4.0*pi/M*int_gsl(int_for_rhohat, (void*)&params, 0.0, exp(logrmin), EPS, err);
    forwardError(*err, __LINE__, 0);

    res = res/norm;
       
#undef logrmin
    
  }
     

   return res;
}
#undef xvir  
#undef EPS



/* DEPRECATED, BECAUSE: when interpolated, rhohat gives unstable 
   integration results. And the interpolation is slower...
*/

#define Nc          600
#define Neta        600
#define logcmin    -5.0
#define logcmax    10.5
#define logetamin -24.5
#define logetamax  14.5
#define EPS        1.0e-6
#define xvir       1.0
double rhohat_haloOLD(cosmo_hm *model, double k, double M, double a, int interp, error **err)
{

  
   double res, rvir;
   double c, f, eta, cieta, sieta, cieta1pc, sieta1pc;
   cosmo_hm_params params;
   
   rvir = r_vir(model, M, a, err);
   forwardError(*err, __LINE__, 0.0);
   c     = concentration(model, M, a, err);
   forwardError(*err, __LINE__, 0.0);

   /* Closed form for NFW cannot be used */
   if (fabs(model->alpha_NFW-1.0)>EPS) goto ninterp;

   if (interp==1) {

      double dlogc, dlogeta, logc, logeta;
      interTable2D *tab;
      int i=0, j=0;

      if (model->rhohat==NULL) {

	 dlogc   = (logcmax - logcmin)/(Nc-1.0);
	 dlogeta = (logetamax - logetamin)/(Neta-1.0);
	 tab     = init_interTable2D(Nc, logcmin, logcmax, dlogc, Neta, logetamin, logetamax,
				     dlogeta, 0.0, 0.0, err);
	 forwardError(*err, __LINE__, 0.0);

	 for (i=0,logc=logcmin; i<Nc; i++,logc+=dlogc) {
	    c = exp(logc);

	    for (j=0,logeta=logetamin; j<Neta; j++,logeta+=dlogeta) {
	       eta = exp(logeta);
	       f = 1/(log(1.0+c) - c/(1.0+c));
	       sm2_cisi(eta, &cieta, &sieta, err);
	       forwardError(*err, __LINE__, 0);
	       sm2_cisi(eta*(1.0+c), &cieta1pc, &sieta1pc, err);
	       forwardError(*err, __LINE__, 0);
	       /* TJ03 (17) */
	       tab->table[i][j] = f*(sin(eta)*(sieta1pc - sieta) + cos(eta)*(cieta1pc - cieta)
				     - sin(eta*c)/(eta*(1.0+c)));
	    }

	 }

	 model->rhohat = tab;

      }

      logc   = log(c);
      logeta = log(k*rvir) - logc;

      if (logc<logcmin || logc>logcmax || logeta<logetamin || logeta>logetamax)
	goto ninterp;
      


      res = interpol2D(model->rhohat, logc, logeta, err);
      forwardError(*err, __LINE__, 0);

   } else {
     
   ninterp:
     /* Closed formula for NFW profile (alpha=1) */
     if (fabs(model->alpha_NFW-1.0)<EPS) {
       
       f = 1.0/(log(1.0+c) - c/(1.0+c));
       eta = k*rvir/c;
       sm2_cisi(eta, &cieta, &sieta, err);                 forwardError(*err, __LINE__, 0);
       sm2_cisi(eta*(1.0+c), &cieta1pc, &sieta1pc, err);   forwardError(*err, __LINE__, 0);
       
       /* TJ03 (17) */
       res = f*(sin(eta)*(sieta1pc - sieta) + cos(eta)*(cieta1pc - cieta)
		- sin(eta*c)/(eta*(1.0+c)));
       
     } else {
       
       double norm;
       
#define logrmin -6.0
       
       params.model  = model;
       params.k     = k;
       params.M     = M;
       params.a     = a;
       params.c     = c;
       
       res          = 0.0;
       
       params.logintegrate = +1;
       res += 4.0*pi/M*int_gsl(int_for_rhohat, (void*)&params, logrmin, log(xvir*rvir), EPS, err);
       forwardError(*err, __LINE__, 0);
       
       params.logintegrate = -1;
       res += 4.0*pi/M*int_gsl(int_for_rhohat, (void*)&params, 0.0, exp(logrmin), EPS, err);
       forwardError(*err, __LINE__, 0);
       
       /* Normalization -> rhohat(k=0) = 1 */
       params.k = 0.0;
       
       params.logintegrate = +1;
       norm = 4.0*pi/M*int_gsl(int_for_rhohat, (void*)&params, logrmin, log(xvir*rvir), EPS, err);
       forwardError(*err, __LINE__, 0);
       params.logintegrate = -1;
       norm += 4.0*pi/M*int_gsl(int_for_rhohat, (void*)&params, 0.0, exp(logrmin), EPS, err);
       forwardError(*err, __LINE__, 0);
       
       res = res/norm;
       
#undef logrmin
       
     }
     
   }
   
   

   return res;
}
#undef Nc
#undef Neta
#undef logcmin
#undef logcmax
#undef logetamin
#undef logetamax
#undef EPS


double int_for_rhohat(double logr, void *params, error **err)
{
  cosmo_hm *model;
  double k, M, a, res, r, c;
  int logintegrate;
  
  model  = ((cosmo_hm_params *)params)->model;
  k      = ((cosmo_hm_params *)params)->k;
  M      = ((cosmo_hm_params *)params)->M;
  a      = ((cosmo_hm_params *)params)->a;
  c      = ((cosmo_hm_params *)params)->c;
  
  logintegrate = ((cosmo_hm_params *)params)->logintegrate;
  
  if (logintegrate==1) r = exp(logr);
  else r = logr;
  res    = r*r*rho_halo(model, r, a, M, c, err)*sinc(r*k);
  forwardError(*err, __LINE__, 0);
  if (logintegrate==1) res *= r;
  
  return res;
}


/* ============================================================== *
 * CS02 (68), M in M_sol/h. k is the order of the bias expansion, *
 * not the scale.						  *
 * ============================================================== */
double bias(cosmo_hm *model, double M, double a, int k, error **err)
{
   double b, eps[3], E[3], a2, deltac, qnusqr;

   testErrorRetVA(model->halo_bias != halo_bias_sc, hm_halo_bias, "Invalid halo bias type %d, has to be %d",
		  *err, __LINE__, 0.0, model->halo_bias, halo_bias_sc);

   if (k==0) return 1.0;

   deltac  = delta_c(model->cosmo, a, err);                              forwardError(*err, __LINE__, 0.0);
   /*  deltac = 1.686;  */
   qnusqr  = model->nmz_a*dsqr(deltac/D_plus(model->cosmo, a, 1, err));   forwardError(*err, __LINE__, 0.0);
   qnusqr /= sigmasqr_M(model, M, err);                                  forwardError(*err, __LINE__, 0.0);
   eps[1]  = (qnusqr - 1.0)/deltac;
   E[1]    = 2.0*model->nmz_p/(deltac*(1.0 + pow(qnusqr, model->nmz_p)));

   switch (k) {

      case 1  : b = 1.0 + eps[1] + E[1];
	        break;

      case 2  : eps[2] = qnusqr/deltac*(qnusqr-3.0)/deltac;
	        E[2]   = ((1.0 + 20.0*model->nmz_p)/deltac + 2.0*eps[1])*E[1];
		a2     = -17.0/21.0;
		b      = 2.0*(1.0 + a2)*(eps[1] + E[1]) + eps[2] + E[2];
		break;

      default : *err = addError(ce_unknown, "bias order too large", *err, __LINE__);
	        b = -1.0;

   }

   return b;
}

/* Tinker et al. 2010 (6)/Table 2 */
double bias_tinker10(cosmo_hm *model, double M, double a, error **err)
{
   double nu, deltac, bb;

   testErrorRetVA(model->halo_bias != halo_bias_tinker10, hm_halo_bias, "Invalid halo bias type %d, has to be %d",
		  *err, __LINE__, 0.0, model->halo_bias, halo_bias_tinker10);

   deltac  = delta_c(model->cosmo, a, err);            forwardError(*err, __LINE__, 0.0);
   nu      = deltac/D_plus(model->cosmo, a, 1, err);   forwardError(*err, __LINE__, 0.0);
   nu     /= sqrt(sigmasqr_M(model, M, err));          forwardError(*err, __LINE__, 0.0);
   

   /* TO DO give choice of Delta_h (defined as overdensity to the mean density) */
   double y = log10(Delta_vir(model, model->nmz_a));
    
   double A = 1.0+0.24*y*exp(-pow(4.0/y,4.0));
   double aa = 0.44*y-0.88;
   double B = 0.183;
   double b = 1.5;
   double C = 0.019+0.107*y+0.19*exp(-pow(4.0/y,4.0));
   double c = 2.4;
   
   bb = 1.0-A*pow(nu,aa)/(pow(nu,aa)+pow(deltac,aa))+B*pow(nu,b)+C*pow(nu,c);

   return bb;
}


/* Tinker et al. 2005 (A1) */
double bias_tinker(cosmo_hm *model, double M, double a, error **err)
{
   double qnusqr, deltac, bb, b, c;

   testErrorRetVA(model->halo_bias != halo_bias_tinker05, hm_halo_bias, "Invalid halo bias type %d, has to be %d",
		  *err, __LINE__, 0.0, model->halo_bias, halo_bias_tinker05);

   deltac  = delta_c(model->cosmo, a, err);                              forwardError(*err, __LINE__, 0.0);
   qnusqr  = model->nmz_a*dsqr(deltac/D_plus(model->cosmo, a, 1, err));   forwardError(*err, __LINE__, 0.0);
   qnusqr /= sigmasqr_M(model, M, err);                                  forwardError(*err, __LINE__, 0.0);
   
   bb  = 0.35;
   c   = 0.8;
   
   b   = sqrt(model->nmz_a)*qnusqr;
   b  += sqrt(model->nmz_a)*bb*pow(qnusqr, 1.0-c);
   b  -= pow(qnusqr, c)/(pow(qnusqr, c) + bb*(1.0-c)*(1.0-c/2.0));
   b  *= 1.0/sqrt(model->nmz_a)/deltac;
   b  += 1.0;

   return b;
}

/* ============================================================ *
 * Returns the bias of halos wrt to the smooth dark matter bg.  *
 * ============================================================ */
double halo_bias(cosmo_hm *model, double M, double a, int k, error **err)
{
   double res;
   
   switch (model->halo_bias) {
      case halo_bias_sc :
	res = bias(model, M, a, k, err);
	forwardError(*err, __LINE__, 0.0);
	break;
   case halo_bias_tinker05 :
     res = bias_tinker(model, M, a, err);
     forwardError(*err, __LINE__, 0.0);
     break;
   case halo_bias_tinker10 :
     res = bias_tinker10(model, M, a, err);
     forwardError(*err, __LINE__, 0.0);
     break;
   default :
     res = 0.0;
     *err = addErrorVA(hm_halo_bias, "Invalid bias type %d", *err, __LINE__, model->halo_bias);
   }
   
  return res;
}

double int_for_bias_norm(double logM, void *params, error **err)
{
   double res, M, a, dp, sM, nu;
   int k;
   cosmo_hm_params *cANDs;
   cosmo_hm *model;

   M     = exp(logM);

   cANDs = (cosmo_hm_params*)params;
   model  = cANDs->model;
   a     = cANDs->a;
   k     = cANDs->i;

   dp    = D_plus(model->cosmo, a, 1, err);             forwardError(*err, __LINE__, 0.0);
   sM    = sqrt(sigmasqr_M(model, M, err));             forwardError(*err, __LINE__, 0.0);
   nu    = delta_c(model->cosmo, a, err)/(dp*sM);       forwardError(*err, __LINE__, 0.0);

   if (model->massfct == j01) {
      res = nufnu_j01(dp * sM)/nu;
   } else {
      res = nufnu(model, nu, 0, err)/nu;                forwardError(*err, __LINE__, 0.0);
   }
   res  *= dnu_dlnM(model, M, a, err);                  forwardError(*err, __LINE__, 0.0);
   res  *= halo_bias(model, M, a, k, err);              forwardError(*err, __LINE__, 0.0);

   return res;
}

/* Returns int(dlogM M^2/rhobar n(M) b(M)) = int(dlogM nu f(nu)/nu dnu/dlogM b(nu). *
 * Used in 2h-term to normalize P2h to P_lin on large scales.			    */
#define EPS 1.0e-5
double bias_norm(cosmo_hm *model, double a, error **err)
{
   double norm;
   cosmo_hm_params cANDs;

   cANDs.model = model;
   cANDs.i    = 1;
   cANDs.a    = a;
   /* MKDEBUG: New, lower limit was logMmin-2 */
   norm = sm2_qromberg(int_for_bias_norm, (void*)&cANDs, logMmin, logMmax, EPS, err);
   forwardError(*err, __LINE__, 0);

   return norm;
}
#undef EPS

/* General mass-integrand for second-order halo terms (1h, 2h) */
#define EPS 1.0e-12
double int_for_M_ij(double logM, void *params, error **err)
{
   int i, j, n;
   double a, b, dndlnM, rhohat, M, Moverrho, res, rhohatfirst=-1.0;
   cosmo_hm_params *cANDs;
   cosmo_hm *model;

   cANDs = (cosmo_hm_params*)params;
   model  = cANDs->model;
   i     = cANDs->i;
   j     = cANDs->j;
   a     = cANDs->a;
   M     = exp(logM);

   /* MKDEBUG New: Tinker bias also here for dm-only */
   b     = halo_bias(model, M, a, i, err);
   forwardError(*err, __LINE__, 0);

   dndlnM = dn_dlnM_uf(M, model, a, err);      forwardError(*err, __LINE__, 0);

   for (n=0,rhohat=1.0; n<j; n++) {
      if (n>0 && fabs(cANDs->kk[n]-cANDs->kk[0])<EPS) {   /* same k as first k          */
         rhohat *= rhohatfirst;
      } else {			    	                /* different k -> recalculate */
         rhohat *= rhohat_halo(model, cANDs->kk[n], M, a, -1, err);
         forwardError(*err, __LINE__, 0);
         rhohatfirst = rhohat;
      }
   }

   for (n=0,Moverrho=1.0; n<j; n++) {
      Moverrho *= M/(Omega_m_halo(model, a, err)*rho_crit_halo(model, a, err));
   }

   res = Moverrho * dndlnM * b * rhohat;
   return res;
}
#undef EPS

/* ============================================================ *
 * CS02 (98). k is a j-dim. vector				*
 * ============================================================ */
double M_ij(cosmo_hm *model, int i, int j, double a, const double *k, error **err)
{
   double Mij;
   cosmo_hm_params params;
   int n;

   testErrorRet(i<0 || i>2 || j<=0 || j>3, ce_unknown, "indices out of range", *err, __LINE__, 0);

   params.model = model;
   params.i    = i;
   params.j    = j;
   params.a    = a;
   params.kk   = malloc(sizeof(double)*j);
   for (n=0; n<j; n++) params.kk[n] = k[n];
   
   Mij = sm2_qromberg(int_for_M_ij, (void*)&params, logMmin, logMmax, 1.e-4, err);
   forwardError(*err, __LINE__, 0);
   free(params.kk);

   return Mij;
}

/* 1-halo term of dark matter power spectrum, k in h/Mpc */
double P1h_dm(cosmo_hm *model, double a, double k, error **err)
{
   double K[2], res;

   K[0] = K[1] = k;
   res = M_ij(model, 0, 2, a, K, err);           forwardError(*err, __LINE__, 0);
   return res;
}

/* 2-halo term of the power spectrum, k in h/Mpc */
double P2h_dm(cosmo_hm *model, double a, double k, error **err)
{
   double p2h;

   p2h = dsqr(M_ij(model, 1, 1, a, &k, err));    forwardError(*err, __LINE__, 0);
   p2h /= dsqr(bias_norm(model, a, err));        forwardError(*err, __LINE__, 0);
   p2h *= P_L(model->cosmo, a, k, err);          forwardError(*err, __LINE__, 0);

   return p2h;
}

int change_Pth(cosmo_hm* avant, cosmo_hm* apres)
{
   if (change_rhohat_halo(avant,apres)) return 1;
   if (change_P_NL(avant->cosmo, apres->cosmo)) return 1;

   if (NCOEQ(avant->cosmo, apres->cosmo, growth) ||
       NCOEQ(avant->cosmo, apres->cosmo, transfer) ||
       NCOEQ(avant->cosmo, apres->cosmo, de_param))
     return 1;

   if (change_Esqr(avant->cosmo, apres->cosmo)) return 1;

   if (NCOCLOSE(avant->cosmo, apres->cosmo, h_100) ||
       NCOCLOSE(avant->cosmo, apres->cosmo, Omega_b) ||
       NCOCLOSE(avant->cosmo, apres->cosmo, n_spec)||
       NCOCLOSE(avant->cosmo, apres->cosmo, normalization) ||
       NCOCLOSE(avant->cosmo, apres->cosmo, Neff_nu_mass))
     return 1;

   if (change_w_de(avant->cosmo, apres->cosmo)) return 1;
   if (change_massfct(avant, apres)) return 1;
   if (change_halo_bias(avant, apres)) return 1;
   if (change_redshift(avant->redshift, apres->redshift)) return 1;

   return 0;
}



double int_gsl(funcwithpars func, void *params, double a, double b, double eps, error **err)
{
  int n = 1000, status;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (n);
  double result, result_err;
  
  gsl_function F;
  F.function = &integrand_gsl;
  
  gsl_int_params p;
  p.func   = func;
  p.err    = err;
  p.params = params;
  F.params = &p;
  
  gsl_set_error_handler_off();
  status = gsl_integration_qag (&F, a, b, eps, eps, n, GSL_INTEG_GAUSS51, w, &result, &result_err);
  forwardError(*err, __LINE__, 0.0);
  
  // THIS IS UNSTABLE, WHY ??
  //testErrorRetVA(status != 0, hm_gsl_int, "gsl integration error, gsl returned with status=%s", *err, __LINE__, 0.0, gsl_strerror (status));
  
  testErrorExitVA(status != 0, hm_gsl_int, "gsl integration error, gsl returned with status=%s", *err, __LINE__, 0.0, gsl_strerror (status));
  
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


/* ==================================================================== *
 * OBSOLETE STUFF
 * ==================================================================== */  










#define EPS 1.0e-8
double xi_dm_NL_OBSOLETE(cosmo_hm *model, double a, double r, error **err)
{

/* ==================================================================== *
 * Real-space correlation function for the non-linear power spectrum.	*
 * ==================================================================== */  

  double val, dk, k;
  
  cosmo_hm_params params; 

  // variables for integration 
  params.r    = r; 
  params.a    = a; 
  params.model = model;

  testErrorRetVA(r<EPS, math_infnan, "Division by zero (r=%g)", *err, __LINE__, 0.0, r);

  k   = k_min;
  val = 0;
  dk  = pi/r/2.0/50;

  while (k+dk<=k_max_HOD) {
    val += int_for_xi_dm_NL_OBSOLETE(k, (void*)&params, err);
    forwardError(*err, __LINE__, 0.0);
    k   += dk;
  }
  val = val*dk/(2.0*pi*pi);

  return val;
}
#undef EPS

double int_for_xi_dm_NL_OBSOLETE(double k, void *params, error **err)
{
  double val;
  cosmo_hm_params *cANDs;
  cosmo_hm *model;
  double a, r;

  cANDs = (cosmo_hm_params *)params;
  a     = cANDs->a;
  r     = cANDs->r;
  model  = cANDs->model;

  val   = k*k*sin(k*r)/(k*r);
  //  fprintf (stderr,"%f\n",a);fflush(stderr);
  val   *= P_NL(model->cosmo, a, k, err);

  forwardError(*err, __LINE__, 0.0);

  return val;
}
