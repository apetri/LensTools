/* ============================================================ *
 * halomodel.c							*
 * Martin Kilbinger, Jean Coupon 2006-2012                      *
 * Refs:							*
 *   - Takada & Jain 2003					*
 *   - Eisenstein & Hu 1998					*
 *   - Percival 2005						*
 *   - Weinberg & Kamionkowski 2002				*
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
			     hod_t HOD, error **err)
{
   cosmo_hm *res;
   double amin;

   res = (cosmo_hm*)malloc_err(sizeof(cosmo_hm), err);
   forwardError(*err, __LINE__, NULL);

   res->redshift = init_redshift(Nzbin, Nnz, nofz, par_nz, NULL, err);
   forwardError(*err, __LINE__, NULL);

   //[jean]------------------------------------------------
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

   res->M_min       = M_min;
   res->M1          = M1;
   res->M0          = M0;
   res->sigma_log_M = sigma_log_M;
   res->alpha       = alpha;
   res->hod         = HOD;

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
			    source->redshift->par_nz, source->zmin,source->zmax,
			    source->cosmo->nonlinear, source->cosmo->transfer,
			    source->cosmo->growth, source->cosmo->de_param, 
			    source->cosmo->normmode, source->c0, source->alpha_NFW, source->beta_NFW,
			    source->massfct, source->halo_bias, source->M_min, source->M1, 
			    source->M0, source->sigma_log_M, source->alpha, source->hod,
			    err);
   forwardError(*err, __LINE__, NULL);

   return res;
}

cosmo_hm *copy_parameters_hm(cosmo_hm *source, error **err)
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

void read_cosmological_parameters_hm(cosmo_hm **self, FILE *F, error **err)
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

   CONFIG_READ(tmp,zmin, d, F, c, err);
   CONFIG_READ(tmp,zmax, d, F, c, err);

   /* Halomodel parameters (dark matter) */
   CONFIG_READ(tmp, alpha_NFW, d, F, c, err);
   CONFIG_READ(tmp, c0, d, F, c, err);
   CONFIG_READ(tmp, beta_NFW, d, F, c, err);
   CONFIG_READ_S(&tmp2, smassfct, s, F, c, err);
   STRING2ENUM(tmp->massfct, tmp2.smassfct, massfct_t, smassfct_t, j, Nmassfct_t, err);

   CONFIG_READ_S(&tmp2, shalo_bias, s, F, c, err);
   STRING2ENUM(tmp->halo_bias, tmp2.shalo_bias, halo_bias_t, shalo_bias_t, j, Nhalo_bias_t, err);

   /* HOD parameters */
   CONFIG_READ(tmp, M_min, d, F, c, err);
   CONFIG_READ(tmp, M1, d, F, c, err);
   CONFIG_READ(tmp, M0, d, F, c, err);
   CONFIG_READ(tmp, sigma_log_M, d, F, c, err);
   CONFIG_READ(tmp, alpha, d, F, c, err);
   CONFIG_READ_S(&tmp2, shod, s, F, c, err);

   STRING2ENUM(tmp->hod, tmp2.shod, hod_t, shod_t, j, Nhod_t, err);

   testErrorRetVA(tmp->hod!=hod_none && tmp->hod!=hamana04 && tmp->hod!=berwein02 && tmp->hod!=berwein02_hexcl,
		  hm_hodtype, "HOD type (%d) unknown", *err, __LINE__,, tmp->hod, hamana04);

   *self = copy_parameters_hm_only(tmp, err);
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
			     9.0, 1.5, 0.13, st2, halo_bias_sc, 1.6e11, 8.0e12, 1.6e11, 0.3, 0.75, hamana04, err);
}
#undef NZBIN
#undef NNZ

void free_parameters_hm(cosmo_hm** self)
{
   cosmo_hm *s;
   int Nzcorr;

   s = *self;

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
	 //case j01 : *nmz_p = *nmz_a = -1.0; break;                /* Jenkins            */
      case j01 : *nmz_p = 0.3; *nmz_a = 1.0/sqrt(2.0); break;                /* Jenkins            */
      default  : *err = addError(ce_unknown, "Wrong or not supported mass function", *err, __LINE__);
	         return;
   }
}

void dump_param_only_hm(cosmo_hm* self, FILE *F)
{
   if (!F) F = stderr;
   fprintf(F, "#     c0  a_NFW  b_NFW  mnz_a  mnz_p       M_min  M1    M0    slogM alpha m b h\n");
   fprintf(F, "# %6.2f % .3f % .3f % .3f % .3f %6.2e %6.2e %6.2e %.3f % .3f %s(%d) %s(%d) %s(%d)\n",
	   self->c0, self->alpha_NFW, self->beta_NFW,
	   self->nmz_a, self->nmz_p, self->M_min, self->M1, self->M0, self->sigma_log_M, self->alpha,
	   smassfct_t(self->massfct), self->massfct, shalo_bias_t(self->halo_bias), self->halo_bias,
	   shod_t(self->hod), self->hod);
}

void dump_param_hm(cosmo_hm* self, FILE *F, error **err)
{
   if (!F) F = stderr;
   dump_param(self->cosmo, F);
   dump_redshift(self->redshift, F, err);
   forwardError(*err, __LINE__,);
   dump_param_only_hm(self, F);
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
double delta_c(cosmo *self, double a, error **err)
{
   double delta_EdS, alpha, deltac;

   delta_EdS = 1.68647;

   /* Per05, D_+ ? */
   /* return delta_EdS; */

   /* WK02 (18) */
   alpha = 0.131 + self->w0_de*(0.555 + self->w0_de*(1.128 + 
           self->w0_de*(1.044 + self->w0_de*0.353)));

   /* KS96 (A.6). Note: typo (0.123) in astro-ph version. Correct in ApJ. */
   //alpha = 0.0123;

   deltac =  delta_EdS*(1. + alpha*log10(Omega_m_a(self, a, -1.0, err)));
   forwardError(*err, __LINE__, 0.0);

   return deltac;
}

/* Bisection for Mstar. */
double bis_Mstar(double logM, void *param, error **err)
{
   cosmo_hm *self;
   double diff, M;
 
   self = (cosmo_hm*)param;
   M     = exp(logM);

   diff  = delta_c(self->cosmo, 1.0, err);     forwardError(*err, __LINE__, 0);
   /* New: omitting D+(a=1, normalised=1) = 1 */
   diff /= sqrt(sigmasqr_M(self, M, err));     forwardError(*err, __LINE__, 0);
   
   return diff-1.0;
}

/* Bisection for Mstar_a. This redshift-dependent routine is not used in general. */
double bis_Mstar_a(double logM, void *param, error **err)
{
   cosmo_hm *self;
   cosmo_hmANDstuff *cANDd;
   double diff, M, a;
 
   cANDd = (cosmo_hmANDstuff*)param;
   self = cANDd->self;
   a    = cANDd->a;
   M     = exp(logM);

   diff  = delta_c(self->cosmo, a, err);     forwardError(*err, __LINE__, 0);
   diff /= D_plus(self->cosmo, a, 1, err);   forwardError(*err, __LINE__, 0);
   diff /= sqrt(sigmasqr_M(self, M, err));   forwardError(*err, __LINE__, 0);

   return diff-1.0;
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
double Mstar(cosmo_hm *self, error **err)
{
   if (self->Mstar<0) {
      self->Mstar = sm2_rtbis(bis_Mstar, logMsmin, logMsmax, xacc, (void*)self, err);
      forwardError(*err, __LINE__,-1);
      self->Mstar = exp(self->Mstar);
      /* self->Mstar = 1.0e14; */
   }

   return self->Mstar;
}

/* Returns M_* in M_sol/h, redshift-dependent! In general, M_star is used (z=0) */
double Mstar_a(cosmo_hm *self, double a, error **err)
{
   cosmo_hmANDstuff cANDd;
   double Mstar;

   cANDd.self = self;
   cANDd.a     = a;

   Mstar = sm2_rtbis(bis_Mstar_a, logMsmin, logMsmax, xacc, (void*)(&cANDd), err);
   forwardError(*err, __LINE__,-1);

   return exp(Mstar);
}
#undef logMsmin
#undef logMsmax
#undef xacc

/* Concentration parameter, M in M_sol/h				 *
 * So far no de dependance!						 *
 * see Dolag et al 2004: alpha=0.1, c to be corrected by D+(de)/D+(LCDM) */
double concentration(cosmo_hm *self, double M, double a, error **err)
{
   double c, Ms;

   Ms = Mstar(self, err);                   forwardError(*err, __LINE__, 0);
   //Ms = 1e14;
   
   c = self->c0*a*pow(M/Ms, -self->beta_NFW);
  
   return c;
}

/* Virial overdensity */
double Delta_vir(cosmo_hm *self, double a)
{
   double w_vir, D, om, ov;

   Omega_a(self->cosmo, a, &om, &ov);
   w_vir = 1.0/om - 1.0;

   /* KS96 (A.5) */
   //D = 18.0*pi*pi*(1.0 + 0.4093*pow(w_vir, 0.9052));

   /* WK02 (16) */
   double aa, b;
   aa = 0.399 - 1.309*(pow(fabs(self->cosmo->w0_de), 0.426) - 1.0);
   b = 0.941 - 0.205*(pow(fabs(self->cosmo->w0_de), 0.938) - 1.0);
   D = 18.0*pi*pi*(1.0 + aa * pow(w_vir, b));
   
   //printf("Delta_vir(z=%g) = %g\n", 1.0/a - 1.0, D);

   return D;
}

/* dsigma^2(R)/dR, R in Mpc/h */
double dsigma_R_sqr_dR(cosmo_hm *self, double R, error **err)
{
   /* Numerical derivative */
   double res, h;
   h = R/20.0;
   /* Problematic ... ? */
   res = (sigma_R_sqr(self, R+h, err)-sigma_R_sqr(self, R-h, err))/(2.0*h);
   forwardError(*err, __LINE__, 0);
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
double nufnu(cosmo_hm *self, double nu, int asymptotic, error **err)
{
   double res, qnusqr;

   if (self->A<0) {

      /* Sets int[f(nu) nu] = 1 */
      self->A = 1.0/(1.0 + pow(2.0, -self->nmz_p)*exp(gammln(0.5-self->nmz_p))/sqrt(pi));
      /* = 0.322 */

      /* matches linear power spectrum on large scales, CH01: with st2. *
       * int[f(nu) nu] about 19% off.					*/
      //self->A = 0.383;

   }

   testErrorRet(nu<0, ce_negative, "nu<0", *err, __LINE__, -1);

   if (asymptotic==0) {
      qnusqr = self->nmz_a*nu*nu;
      res = sqrt(2.0/pi*qnusqr)*(1.0 + pow(qnusqr, -self->nmz_p))*exp(-qnusqr/2.0);
   } else if (asymptotic==1) {
      res = pow(nu, 0.5-self->nmz_p);
   } else {
      *err = addError(ce_wrongValue, "The flag 'asymptotic' has to be 0 or 1", *err, __LINE__);
      return 0.0;
   }

   return res*self->A;
}

/* Mass function fo j01. *
 * x = D+ * sigma_m      */
double nufnu_j01(double x)
{
   return 0.315 * exp(-pow(fabs(log(1.0/(x)) + 0.61), 3.8));
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
double sigma_R_sqr(cosmo_hm *self, double R, error **err)
{
   double h, a, b, RR, dlogR, res, norm;
   cosmoANDdouble cANDd;
   int i;

   if (self->sigRsqr==NULL) {

      self->sigRsqr = init_splineTable(Ns, err);
      forwardError(*err, __LINE__, 0);

      cANDd.self = self->cosmo;
      dlogR = (log(Rmax)-log(Rmin))/(Ns-1.0);
      norm  = 9.0*dsqr(self->cosmo->sigma_8)/sigma_8_sqr(self->cosmo, err);
      forwardError(*err, __LINE__, 0);
      for (i=1,RR=Rmin; i<=Ns; i++,RR*=exp(dlogR)) {
	 cANDd.r = RR;
	 self->sigRsqr->x[i] = log(RR);
	 self->sigRsqr->y[i] = norm/(2.0*pi*pi)*sm2_qromberg(int_for_sigma_R, (void*)&cANDd, 
							     log(k_min), log(kmax), eps, err);
	 forwardError(*err, __LINE__, 0);
	 self->sigRsqr->y[i] = log(self->sigRsqr->y[i]);
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
      self->sigRsqr->yp1 = (a-b)/(2.0*h)*Rmin/exp(self->sigRsqr->y[1]);
      self->sigRsqr->yp1 = 0.0;

      h = Rmax/20.0;
      cANDd.r = Rmax+h;
      a = norm/(2.0*pi*pi)*sm2_qromberg(int_for_sigma_R, (void*)&cANDd, log(k_min),
					log(kmax), eps, err);
      forwardError(*err, __LINE__, 0);
      cANDd.r = Rmax-h;
      b = norm/(2.0*pi*pi)*sm2_qromberg(int_for_sigma_R, (void*)&cANDd, log(k_min), 
					log(kmax), eps, err);
      forwardError(*err, __LINE__, 0);
      self->sigRsqr->ypn = (a-b)/(2.0*h)*Rmax/exp(self->sigRsqr->y[Ns]);

      sm2_spline(self->sigRsqr->y, self->sigRsqr->x, self->sigRsqr->n, self->sigRsqr->yp1, 
		 self->sigRsqr->ypn, self->sigRsqr->y2, err);
      forwardError(*err, __LINE__, 0);
   }

   testErrorRetVA(R<Rmin||R>Rmax, ce_interpoloutofrange, "R = %g out of range [%g;%g]",
		  *err, __LINE__, 0.0, R, Rmin, Rmax);

   sm2_splint(self->sigRsqr->x, self->sigRsqr->y, self->sigRsqr->y2, self->sigRsqr->n, log(R), &res, err);
   forwardError(*err, __LINE__, 0.0);
   return exp(res);
}
#undef Rmin
#undef Rmax
#undef Ns
#undef eps

/* Calls sigma(R)^2, with M (in M_sol/h) mass in sphere of radius R */
double sigmasqr_M(cosmo_hm *self, double M, error **err)
{
   double R, res;

   R = cbrt(3.0*M/(4.0*pi*self->cosmo->Omega_m*rho_c0));
   res = sigma_R_sqr(self, R, err);
   forwardError(*err, __LINE__, 0);

   return res;
}

/* Returns dsigma^-1/dlnM. M in M_sol/h */
double dsigma_m1_dlnM(cosmo_hm *self, double M, error **err)
{
   double res, sigma, R, rhobar;

   rhobar = self->cosmo->Omega_m*rho_c0;
   R      = cbrt(3.0*M/(4.0*pi*rhobar));
   sigma  = sqrt(sigma_R_sqr(self, R, err));             forwardError(*err, __LINE__, 0);
   res    = M/(2.0*sigma*sigma*sigma);
   res   *= dsigma_R_sqr_dR(self, R, err) * cbrt(1.0/(36.0*pi*rhobar*M*M));
   forwardError(*err, __LINE__, 0);

   return -res;
}

/* M in M_sol/h */
double dnu_dlnM(cosmo_hm *self, double M, double a, error **err)
{
  double res;
  
  res  = delta_c(self->cosmo, a, err);   forwardError(*err, __LINE__, 0.0);
  res /= D_plus(self->cosmo, a, 1, err); forwardError(*err, __LINE__, 0.0);
  res *= dsigma_m1_dlnM(self, M, err);   forwardError(*err, __LINE__, 0.0);
  
  return res;
}

double dn_dlnM_lnM(double logM, void *intpar, error **err)
{
  double res, M;
  M = exp(logM);
  res = M*dn_dlnM(M, intpar, err);
  forwardError(*err, __LINE__, 0);
  return res;
}

/* Mass function, M in M_sol/h, user-friendly version */
double dn_dlnM_uf(double M, cosmo_hm *self, double a, error **err)
{
   cosmo_hmANDstuff2 cANDs2;
   double res;

   cANDs2.self       = self;
   cANDs2.a          = a;
   cANDs2.asymptotic = 0;

   res = dn_dlnM(M, (void*)&cANDs2, err);
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
double dn_dlnM(double M, void *intpar, error **err)
{
  cosmo_hmANDstuff2 *cANDs;
  cosmo_hm *self;
  double res = 0.0, nu, a, dp, sM, nfn, dnudlnM;
  int asym;   

  testErrorRetVA(M<0.0, ce_negative, "Mass M = %g, has to be positive", *err, __LINE__, 0.0, M);

  cANDs = (cosmo_hmANDstuff2*)intpar;
  self  = cANDs->self;
  a     = cANDs->a;
  asym  = cANDs->asymptotic;

  if (self->massfct != ps && self->massfct != st && self->massfct != st2 && self->massfct != j01) {
     *err = addErrorVA(ce_wrongValue, "Unknown mass function type %d", *err, __LINE__, self->massfct);
     return 0.0;
  }
  
  dp = D_plus(self->cosmo, a, 1, err);          forwardError(*err, __LINE__, 0);
  sM = sqrt(sigmasqr_M(self, M, err));          forwardError(*err, __LINE__, 0);
  nu = delta_c(self->cosmo, a, err)/(dp*sM);    forwardError(*err, __LINE__, 0);
    
  testErrorRet(nu<0.0, ce_negative, "nu not positive",      *err, __LINE__, 0.0);
  testErrorRet(dp<0.0, ce_negative, "D+ not positive",      *err, __LINE__, 0.0);
  testErrorRet(sM<0.0, ce_negative, "sigma_M not positive", *err, __LINE__, 0.0);

  if (self->massfct == j01) {
     /* Jenkins 01: nu f(nu) = f(sigma) */
     /* MKDEBUG TODO: sigma or D+ * sigma ??? */
     nfn = nufnu_j01(dp * sM);
  } else {
     nfn = nufnu(self, nu, asym, err);  forwardError(*err, __LINE__, 0);
  }

  dnudlnM = dnu_dlnM(self, M, a, err);   forwardError(*err, __LINE__, 0);
  res     = self->cosmo->Omega_m*rho_c0/M * nfn/nu * dnudlnM;
    
  /*
    // Old version: Bug for j01
    dp   = D_plus(self->cosmo, a, 1, err);            forwardError(*err, __LINE__, 0);
    sM   = sqrt(sigmasqr_M(self, M, err));            forwardError(*err, __LINE__, 0.0);
    sM  *= dp;
    res  = self->cosmo->Omega_m*rho_c0/M;
    res *= 0.315 * exp(-pow(fabs(log(1.0/sM) + 0.61), 3.8));
    res *= dsigma_m1_dlnM(self, M, err);              forwardError(*err, __LINE__, 0.0);
    res /= dp;

    res *= sM;  // New: this was omitted earlier
  */

  testErrorRetVA(res<0.0, ce_negative,  "dn_dlnM=%g negative", *err, __LINE__, 0.0, res);
  testErrorRet(!finite(res), ce_infnan, "dn_dlnM inf or nan",  *err, __LINE__, 0.0);

  return res;
}



#define xmin 1.0e-4
double rho_halo(cosmo_hm *model, double r, double M, double a, double *r_vir, error **err){
  /* Returns NFW profile in h^2 M_sol/Mpc^3. r in Mpc/h, M in M_sol/h. Normalized to rho_c0.
   * Set r_vir = -1.0 to have it computed in rho_halo(..., &r_vir,...). 
   * TO DO: compute rho_s and r_vir outside: 2-elements array or NULL.
   * Merge rho_halo and DeltaSigma
   */
  
  double c, rhobar, rsqrrho, x, Dvir, rho_s;
  
  rhobar = rho_c0*model->cosmo->Omega_m;
  Dvir   = Delta_vir(model, a);
  c      = concentration(model, M, a, err);
  forwardError(*err, __LINE__, 0.0);
  
  /* Normalisation */
  rho_s  = Dvir/3.0*rhobar*c*c*c/(log(1.0+c)-c/(1.0+c));
  
  if (*r_vir < 0.0) {
    *r_vir  = 3.0*M/(4.0*pi*rhobar*Dvir);
    *r_vir  = cbrt(*r_vir);
    testErrorRet(!finite(*r_vir), ce_infnan, "inf or nan encountered (r_vir)", *err, __LINE__, 0);
    testErrorRet(*r_vir<=0, ce_infnan, "Division by zero (r_vir)", *err, __LINE__, 0);
  }
  
  x = r*c / *r_vir;
  if (x > xmin) {
    rsqrrho = pow(x, -model->alpha_NFW)*pow(1.0+x, model->alpha_NFW-3.0);
  } else { /* Approximation finite for r=0 for alpha_NFW<2 */
    rsqrrho = pow(x, 2.0-model->alpha_NFW);
  }
  
  testErrorRet(!finite(rsqrrho), ce_infnan, "inf or nan encountered", *err, __LINE__, 0);
  
  return rho_s*rsqrrho;
}
#undef xmin


#define eps 1.0e-10
double DeltaSigma(cosmo_hm *model, double r, const double M, const double a, error **err){
  /* Returns Delta Sigma(r) in h M_sol/pc^2 (gg lensing convention, r in Mpc/h, M in M_sol/h).
   * See Wright & Brainerd (2000), Eqs 11-16.
   */
  double result;
  
  double rhobar = rho_c0*model->cosmo->Omega_m;
  double Dvir   = Delta_vir(model, a);
  double c      = concentration(model, M, a, err);
  double rho_s  = Dvir/3.0*rhobar*c*c*c/(log(1.0+c)-c/(1.0+c)); 
  double r_vir  = cbrt(3.0*M/(4.0*pi*rhobar*Dvir));
  double r_s    = r_vir/c;
  
  double x      = r/r_s;
  
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

double int_for_rhohat(double logr, void *intpar, error **err)
{
   cosmo_hmANDstuff3 *cANDs;
   cosmo_hm *self;
   double k, M, a, r_vir, res, r;
   int logintegrate;

   cANDs = (cosmo_hmANDstuff3*)intpar;
   self  = cANDs->self;
   k     = cANDs->k;
   M     = cANDs->M;
   a     = cANDs->a;
   r_vir = cANDs->r_vir;
   logintegrate = cANDs->logintegrate;

   if (logintegrate==1) r = exp(logr);
   else r = logr;
   res    = r*r*rho_halo(self, r, M, a, &r_vir, err)*sinc(r*k);
   forwardError(*err, __LINE__, 0);
   if (logintegrate==1) res *= r;

   return res;
}

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
 * ============================================================ */
#define Nc          200
#define Neta        200
#define logcmin    -5.0
#define logcmax    10.5
#define logetamin -24.5
#define logetamax  14.5
#define EPS        1.0e-6
/* Integration of rho up to xvir*r_vir. Should be unity to be consistent with closed NFW formula. *
 * Usually called with interp=1.								  */
#define xvir        1.0
double rhohat_halo(cosmo_hm *self, double k, double M, double a, int interp, error **err)
{
   double res, r_vir=0.0, rhobar;
   double c, f, eta, cieta, sieta, cieta1pc, sieta1pc;
   cosmo_hmANDstuff3 intpar;

   rhobar = rho_c0*self->cosmo->Omega_m;
   r_vir  = cbrt(3.0*M/(4.0*pi*rhobar*Delta_vir(self, a)));
   c = concentration(self, M, a, err);
   forwardError(*err, __LINE__, 0);

   /* Closed form for NFW cannot be used */
   if (fabs(self->alpha_NFW-1.0)>EPS) goto ninterp;

   if (interp==1) {

      double dlogc, dlogeta, logc, logeta;
      interTable2D *tab;
      int i=0, j=0;

      if (self->rhohat==NULL) {

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

	 self->rhohat = tab;

      }

      logc   = log(c);
      logeta = log(k*r_vir) - logc;

      if (logc<logcmin || logc>logcmax || logeta<logetamin || logeta>logetamax)
	goto ninterp;

      res = interpol2D(self->rhohat, logc, logeta, err);
      forwardError(*err, __LINE__, 0);

   } else {

   ninterp:
      /* Closed formula for NFW profile (alpha=1) */
      if (fabs(self->alpha_NFW-1.0)<EPS) {

	 f = 1.0/(log(1.0+c) - c/(1.0+c));
	 eta = k*r_vir/c;
	 sm2_cisi(eta, &cieta, &sieta, err);                 forwardError(*err, __LINE__, 0);
	 sm2_cisi(eta*(1.0+c), &cieta1pc, &sieta1pc, err);   forwardError(*err, __LINE__, 0);

	 /* TJ03 (17) */
	 res = f*(sin(eta)*(sieta1pc - sieta) + cos(eta)*(cieta1pc - cieta)
		  - sin(eta*c)/(eta*(1.0+c)));

      } else {

	 double norm;

#define logrmin -6.0

	 intpar.self  = self;
	 intpar.k     = k;
	 intpar.M     = M;
	 intpar.a     = a;
	 //printf("3 r_vir = %g\n", r_vir);
	 intpar.r_vir = r_vir;
	 res          = 0.0;

	 intpar.logintegrate = +1;
	 res += 4.0*pi/M*sm2_qromberg(int_for_rhohat, (void*)&intpar, logrmin, log(xvir*r_vir), EPS, err);
	 forwardError(*err, __LINE__, 0);
	 intpar.logintegrate = -1;
	 res += 4.0*pi/M*sm2_qromberg(int_for_rhohat, (void*)&intpar, 0.0, exp(logrmin), EPS, err);
	 forwardError(*err, __LINE__, 0);

	 /* Normalization -> rhohat(k=0) = 1 */
	 intpar.k = 0.0;

	 intpar.logintegrate = +1;
	 norm = 4.0*pi/M*sm2_qromberg(int_for_rhohat, (void*)&intpar, logrmin, log(xvir*r_vir), EPS, err);
	 forwardError(*err, __LINE__, 0);
	 intpar.logintegrate = -1;
	 norm += 4.0*pi/M*sm2_qromberg(int_for_rhohat, (void*)&intpar, 0.0, exp(logrmin), EPS, err);
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

/* ============================================================== *
 * CS02 (68), M in M_sol/h. k is the order of the bias expansion, *
 * not the scale.						  *
 * ============================================================== */
double bias(cosmo_hm *self, double M, double a, int k, error **err)
{
   double b, eps[3], E[3], a2, deltac, qnusqr;

   testErrorRetVA(self->halo_bias != halo_bias_sc, hm_halo_bias, "Invalid halo bias type %d, has to be %d",
		  *err, __LINE__, 0.0, self->halo_bias, halo_bias_sc);

   if (k==0) return 1.0;

   deltac  = delta_c(self->cosmo, a, err);                              forwardError(*err, __LINE__, 0.0);
   /*  deltac = 1.686;  */
   qnusqr  = self->nmz_a*dsqr(deltac/D_plus(self->cosmo, a, 1, err));   forwardError(*err, __LINE__, 0.0);
   qnusqr /= sigmasqr_M(self, M, err);                                  forwardError(*err, __LINE__, 0.0);
   eps[1]  = (qnusqr - 1.0)/deltac;
   E[1]    = 2.0*self->nmz_p/(deltac*(1.0 + pow(qnusqr, self->nmz_p)));

   switch (k) {

      case 1  : b = 1.0 + eps[1] + E[1];
	        break;

      case 2  : eps[2] = qnusqr/deltac*(qnusqr-3.0)/deltac;
	        E[2]   = ((1.0 + 20.0*self->nmz_p)/deltac + 2.0*eps[1])*E[1];
		a2     = -17.0/21.0;
		b      = 2.0*(1.0 + a2)*(eps[1] + E[1]) + eps[2] + E[2];
		break;

      default : *err = addError(ce_unknown, "bias order too large", *err, __LINE__);
	        b = -1.0;

   }

   return b;
}

/* Tinker et al. 2005 (A1) */
double bias_tinker(cosmo_hm *self, double M, double a, error **err)
{
   double qnusqr, deltac, bb, b, c;

   testErrorRetVA(self->halo_bias != halo_bias_tinker05, hm_halo_bias, "Invalid halo bias type %d, has to be %d",
		  *err, __LINE__, 0.0, self->halo_bias, halo_bias_tinker05);

   deltac  = delta_c(self->cosmo, a, err);                              forwardError(*err, __LINE__, 0.0);
   qnusqr  = self->nmz_a*dsqr(deltac/D_plus(self->cosmo, a, 1, err));   forwardError(*err, __LINE__, 0.0);
   qnusqr /= sigmasqr_M(self, M, err);                                  forwardError(*err, __LINE__, 0.0);
   
   bb  = 0.35;
   c   = 0.8;
   
   b   = sqrt(self->nmz_a)*qnusqr;
   b  += sqrt(self->nmz_a)*bb*pow(qnusqr, 1.0-c);
   b  -= pow(qnusqr, c)/(pow(qnusqr, c) + bb*(1.0-c)*(1.0-c/2.0));
   b  *= 1.0/sqrt(self->nmz_a)/deltac;
   b  += 1.0;

   return b;
}

/* ============================================================ *
 * Returns the bias of halos wrt to the smooth dark matter bg.  *
 * ============================================================ */
double halo_bias(cosmo_hm *self, double M, double a, int k, error **err)
{
   double res;

   switch (self->halo_bias) {
      case halo_bias_sc :
	 res = bias(self, M, a, k, err);
	 forwardError(*err, __LINE__, 0.0);
	 break;
      case halo_bias_tinker05 :
	 res = bias_tinker(self, M, a, err);
	 forwardError(*err, __LINE__, 0.0);
	 break;
      default :
	 res = 0.0;
	 *err = addErrorVA(hm_halo_bias, "Invalid bias type %d", *err, __LINE__, self->halo_bias);
   }
 
  return res;
}

double int_for_bias_norm(double logM, void *intpar, error **err)
{
   double res, M, a, dp, sM, nu;
   int k;
   cosmo_hmANDstuff_dm *cANDs;
   cosmo_hm *self;

   M     = exp(logM);

   cANDs = (cosmo_hmANDstuff_dm*)intpar;
   self  = cANDs->self;
   a     = cANDs->a;
   k     = cANDs->i;

   dp    = D_plus(self->cosmo, a, 1, err);             forwardError(*err, __LINE__, 0.0);
   sM    = sqrt(sigmasqr_M(self, M, err));             forwardError(*err, __LINE__, 0.0);
   nu    = delta_c(self->cosmo, a, err)/(dp*sM);       forwardError(*err, __LINE__, 0.0);

   if (self->massfct == j01) {
      res = nufnu_j01(dp * sM)/nu;
   } else {
      res = nufnu(self, nu, 0, err)/nu;                forwardError(*err, __LINE__, 0.0);
   }
   res  *= dnu_dlnM(self, M, a, err);                  forwardError(*err, __LINE__, 0.0);
   res  *= halo_bias(self, M, a, k, err);              forwardError(*err, __LINE__, 0.0);

   return res;
}

/* Returns int(dlogM M^2/rhobar n(M) b(M)) = int(dlogM nu f(nu)/nu dnu/dlogM b(nu). *
 * Used in 2h-term to normalize P2h to P_lin on large scales.			    */
#define EPS 1.0e-5
double bias_norm(cosmo_hm *self, double a, error **err)
{
   double norm;
   cosmo_hmANDstuff_dm cANDs;

   cANDs.self = self;
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
double int_for_M_ij(double logM, void *intpar, error **err)
{
   int i, j, n;
   double a, b, dndlnM, rhohat, M, Moverrho, res, rhohatfirst=-1.0;
   cosmo_hmANDstuff_dm *cANDs;
   cosmo_hm *self;

   cANDs = (cosmo_hmANDstuff_dm*)intpar;
   self  = cANDs->self;
   i     = cANDs->i;
   j     = cANDs->j;
   a     = cANDs->a;
   M     = exp(logM);

   /* MKDEBUG New: Tinker bias also here for dm-only */
   b     = halo_bias(self, M, a, i, err);
   forwardError(*err, __LINE__, 0);
   
   dndlnM = dn_dlnM_uf(M, self, a, err);      forwardError(*err, __LINE__, 0);

   for (n=0,rhohat=1.0; n<j; n++) {
      if (n>0 && fabs(cANDs->k[n]-cANDs->k[0])<EPS) {   /* same k as first k          */
	 rhohat *= rhohatfirst;
      } else {			    	                /* different k -> recalculate */
	 rhohat *= rhohat_halo(self, cANDs->k[n], M, a, 1, err);
	 forwardError(*err, __LINE__, 0);
	 rhohatfirst = rhohat;
      }
   }

   for (n=0,Moverrho=1.0; n<j; n++) {
      Moverrho *= M/(self->cosmo->Omega_m*rho_c0);
   }

   res = Moverrho * dndlnM * b * rhohat;
   return res;
}
#undef EPS

/* ============================================================ *
 * CS02 (98). k is a j-dim. vector				*
 * ============================================================ */
double M_ij(cosmo_hm *self, int i, int j, double a, const double *k, error **err)
{
   double Mij;
   cosmo_hmANDstuff_dm intpar;
   int n;

   testErrorRet(i<0 || i>2 || j<=0 || j>3, ce_unknown, "indices out of range", *err, __LINE__, 0);

   intpar.self = self;
   intpar.i    = i;
   intpar.j    = j;
   intpar.a    = a;
   intpar.k    = malloc(sizeof(double)*j);
   for (n=0; n<j; n++) intpar.k[n] = k[n];
   
   Mij = sm2_qromberg(int_for_M_ij, (void*)&intpar, logMmin, logMmax, 1.e-4, err);
   forwardError(*err, __LINE__, 0);
   free(intpar.k);

   return Mij;
}

/* 1-halo term of dark matter power spectrum, k in h/Mpc */
double P1h_dm(cosmo_hm *self, double a, double k, error **err)
{
   double K[2], res;

   K[0] = K[1] = k;
   res = M_ij(self, 0, 2, a, K, err);           forwardError(*err, __LINE__, 0);
   return res;
}

/* 2-halo term of the power spectrum, k in h/Mpc */
double P2h_dm(cosmo_hm *self, double a, double k, error **err)
{
   double p2h;

   p2h = dsqr(M_ij(self, 1, 1, a, &k, err));    forwardError(*err, __LINE__, 0);
   p2h /= dsqr(bias_norm(self, a, err));        forwardError(*err, __LINE__, 0);
   p2h *= P_L(self->cosmo, a, k, err);          forwardError(*err, __LINE__, 0);

   //printf("MKDEBUG %g %g %g\n", k, dsqr(M_ij(self, 1, 1, a, &k, err)), dsqr(bias_norm(self, a, err)));

   return p2h;
}

/* ==================================================================== *
 * Real-space correlation function for the non-linear power spectrum.	*
 * ==================================================================== */
#define EPS 1.0e-8
double xi_dm_NL(cosmo_hm *self, double a, double r, error **err)
{
  double val, dk, k;
  
  cosmo_hmANDhjmcc2 intpar; 

  // variables for integration 
  intpar.r    = r; 
  intpar.a    = a; 
  intpar.self = self;

  testErrorRetVA(r<EPS, math_infnan, "Division by zero (r=%g)", *err, __LINE__, 0.0, r);

  k   = k_min;
  val = 0;
  dk  = pi/r/2.0/50;

  while (k+dk<=k_max_HOD) {
    val += int_for_xi_dm_NL(k, (void*)&intpar, err);
    forwardError(*err, __LINE__, 0.0);
    k   += dk;
  }
  val = val*dk/(2.0*pi*pi);

  return val;
}
#undef EPS

double int_for_xi_dm_NL(double k, void *intpar, error **err)
{
  double val;
  cosmo_hmANDhjmcc2 *cANDs;
  cosmo_hm *self;
  double a, r;

  cANDs = (cosmo_hmANDhjmcc2 *)intpar;
  a     = cANDs->a;
  r     = cANDs->r;
  self  = cANDs->self;

  val   = k*k*sin(k*r)/(k*r);
  //  fprintf (stderr,"%f\n",a);fflush(stderr);
  val   *= P_NL(self->cosmo, a, k, err);

  forwardError(*err, __LINE__, 0.0);

  return val;
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

/* 1h+2h=total dm power spectrum, k in h/Mpc */
#define eps_a  1.0e-5
#define N_k_hm 20
double Pth_dm(cosmo_hm *self, double a, double k, error **err)
{
   double res;
   double dlogk, logk, da, aa, kk, logkmin, logkmax;
   int i, j;

   if (self->Pthdm==NULL) {

      logkmin = log(k_min);
      logkmax = log(k_max_HOD);
      dlogk   = (logkmax - logkmin)/(N_k_hm-1.0);
      da      = (1.0-self->cosmo->a_min)/(Na_hm-1.0);

      /* MK: TODO upper extrapolation index ??? */
      self->Pthdm = init_interTable2D(Na_hm, self->cosmo->a_min, 1.0, da, N_k_hm, logkmin, logkmax, dlogk,
				      self->cosmo->n_spec, -3.0, err);
      forwardError(*err, __LINE__, 0.0);

      for (i=0,aa=self->cosmo->a_min; i<Na_hm; i++,aa+=da) {
	 //fprintf(stderr, "%2d ", i);
	 for (j=0,logk=logkmin; j<N_k_hm; j++,logk+=dlogk) {
	    kk = exp(logk);

	    res  = P1h_dm(self, aa, kk, err);   forwardError(*err, __LINE__, 0);
	    res += P2h_dm(self, aa, kk, err);   forwardError(*err, __LINE__, 0);
	    self->Pthdm->table[i][j] = log(res);
	 }
      }
      //fprintf(stderr, "\n");
   }

   logk = log(k);
   res  = interpol2D(self->Pthdm, a, logk, err); forwardError(*err, __LINE__, 0.0);
   return exp(res);
}
#undef eps_a
#undef N_k_hm
