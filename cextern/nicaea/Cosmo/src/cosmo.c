/* ============================================================ *
 * cosmo.c							*
 * Martin Kilbinger, Karim Benabed 2006 - 2010			*
 *								*
 * Early versions of the power spectrum code (linear +          *
 * PD96-nonlinear) are based on a program by P. Schneider.	*
 * The Smith et al. (2003) non-linear power spectrum is taken   *
 * from the halofit program by R.E.Smith and J.Peacock.		*
 *								*
 * References:							*
 * - Press et al., Numerical Recipes in C (NR)			*
 * - Carroll, Press & Turner 1992 (CPT)				*
 * - Peacock & Dodds 1994, 1996 (PD2, PD)			*
 * - Bardeen, Bond, Kaiser & Szalay 1986 (BBKS)			*
 * - Bond & Efstathiou 1984					*
 * - Smith et al. 2003						*
 * - Eisenstein & Hu 1998					*
 * - Sugijyama 1995						*
 * - Heath et al. 1997						*
 * - Wetterich 2004						*
 * - Percival 2005						*
 * - Linder & Jenkins 2003					*
 * - Linder 2003                                                *
 * - Jassal, Bagla & Padmanabhan 2005				*
 * - Takahashi et al. 2012					*
 * ============================================================ */


#include "cosmo.h"


/* ============================================================ *
 * Creates a new cosmo model with all tables set to 0.		*
 * ============================================================ */

cosmo* init_parameters(double OMEGAM, double OMEGADE, double W0_DE, double W1_DE,
		       double *W_POLY_DE, int N_POLY_DE,
		       double H100, double OMEGAB, double OMEGANUMASS,
		       double NEFFNUMASS, double NORM, double NSPEC,
		       nonlinear_t NONLINEAR, transfer_t TRANSFER, growth_t GROWTH,
		       de_param_t DEPARAM, norm_t normmode, double AMIN,
		       error **err)
{
   cosmo* res;
   int N_a_min;

   res = (cosmo*) malloc_err(sizeof(cosmo), err);   forwardError(*err, __LINE__, NULL);

   res->Omega_m       = OMEGAM;
   res->Omega_de      = OMEGADE;

   if (DEPARAM != poly_DE) {
      res->w0_de         = W0_DE;
      res->w1_de         = W1_DE;
      res->w_poly_de     = NULL;
      res->N_poly_de     = 0;
   } else {
      set_w_poly_de(&res->w_poly_de, &res->N_poly_de, W_POLY_DE, N_POLY_DE, 1, err);

      /* MDKDEBUG: know what we are doing here! */
      res->w0_de = res->w_poly_de[0];
      res->w1_de = -1.0e30;
   }

   res->h_100         = H100;
   res->Omega_b       = OMEGAB;
   res->Omega_nu_mass = OMEGANUMASS;
   res->Neff_nu_mass  = NEFFNUMASS;
   res->normalization = NORM;
   res->n_spec        = NSPEC;

   res->nonlinear     = NONLINEAR;
   res->transfer      = TRANSFER;
   res->growth        = GROWTH;
   res->de_param      = DEPARAM;
   res->normmode      = normmode;
	 
   if (AMIN>0.0) res->a_min = AMIN;
   else res->a_min = a_minmin;

   /* Reset pre-computed numbers and tables */

   /* New 06/2014: Set N_a such that da >= 0.001 */
   N_a_min = (int)((1.0 - res->a_min) / 0.001) + 1;
   if ( N_a_min < _N_a) {
      res->N_a = _N_a;
   } else {
      res->N_a = N_a_min;
   }

   //printf("N_a, da = %d %g\n", res->N_a, (1.0 - res->a_min) / (res->N_a - 1)); //, N_a_min);

   res->linearGrowth         = NULL;
   res->growth_delta0        = -1;
   res->transferFct          = NULL;
   res->transfer_alpha_Gamma = -1;
   res->transfer_s           = -1;
   res->transferBE           = NULL;
   res->cmp_sigma8           = -1;
   res->P_NL                 = NULL;
   res->slope                = NULL;
   res->w                    = NULL;
   //res->k_NL                 = NULL;
   res->ystar_allz           = NULL;

   //dump_param(res, stdout);

   set_norm(res, err);                    forwardError(*err, __LINE__, NULL);

   // MKDEBUG NEW
   consistency_parameters(res, err);      forwardError(*err, __LINE__, NULL);

   return res;
}

/* ============================================================ *
 * Checks consistency of parameters. This function is called    *
 * in init_parameters and updateFrom.				*
 * ============================================================ */
void consistency_parameters(const cosmo *self, error **err)
{
   testErrorRetVA(self->h_100 <= 0.0, ce_negative,
		  "Negateive Hubble parameter h_100 = %g", *err, __LINE__,, self->h_100);
   testErrorRetVA(self->Omega_m < 0.0, ce_negative,
		  "Negative matter density parameter Omega_m = %g", *err, __LINE__,, self->Omega_m);
   testErrorRetVA(self->Omega_b < 0.0, ce_negative,
		  "Negative baryon density parameter Omega_b = %g", *err, __LINE__,, self->Omega_b);
   testErrorRetVA(self->Omega_nu_mass<0.0, ce_negative,
		  "Negative density parameter of massive neutrinos Omega_nu_mass = %g", *err,
		  __LINE__,, self->Omega_nu_mass);
   testErrorRetVA(self->normalization < 0.0, ce_negative,
		  "Negative normalization = %g", *err, __LINE__,, self->normalization);

   /* MKDEBUG: New. Without this test here, there are NaNs in Tsqr */
   testErrorRetVA(self->Omega_m < self->Omega_b, ce_omega, "Omega_m (%g) < Omega_b (%g)", *err, __LINE__,,
		  self->Omega_m, self->Omega_b);
}

cosmo* copy_parameters_only(cosmo* source, error **err)
{
   cosmo* res;

   res = init_parameters(source->Omega_m,source->Omega_de, source->w0_de, source->w1_de,
			 source->w_poly_de, source->N_poly_de,
			 source->h_100, source->Omega_b, source->Omega_nu_mass, source->Neff_nu_mass,
			 source->normalization, source->n_spec, source->nonlinear, source->transfer,
			 source->growth, source->de_param, source->normmode, source->a_min, err);
   forwardError(*err,__LINE__,NULL);

   return res;
}

cosmo* copy_parameters(cosmo* source, error **err)
{
   cosmo* res;

   res = init_parameters(source->Omega_m, source->Omega_de, source->w0_de, source->w1_de,
			 source->w_poly_de, source->N_poly_de, source->h_100,
			 source->Omega_b, source->Omega_nu_mass, source->Neff_nu_mass,
			 source->normalization, source->n_spec, source->nonlinear,source->transfer,
			 source->growth, source->de_param, source->normmode, source->a_min, err);

   forwardError(*err,__LINE__,NULL);
   res->N_a = source->N_a;
   res->linearGrowth = copy_interTable(source->linearGrowth,err);
   forwardError(*err,__LINE__,NULL);
   res->transferFct = copy_interTable(source->transferFct,err);
   forwardError(*err,__LINE__,NULL);
   res->transferBE = copy_interTable(source->transferBE,err);
   forwardError(*err,__LINE__,NULL);
   res->P_NL = copy_interTable2D(source->P_NL,err);
   forwardError(*err,__LINE__,NULL);
   res->slope = copy_interTable(source->slope,err);
   forwardError(*err,__LINE__,NULL);
   res->w = copy_interTable(source->w,err);
   forwardError(*err,__LINE__,NULL);

   //res->k_NL = copy_interTable(source->k_NL, err);
   //forwardError(*err,__LINE__,NULL);

   //  DEBUGGING THIS doesn't work with halomodel 
   //memcpy(res->ystar_allz, source->ystar_allz, fr_rs*fr_nsim);

   res->growth_delta0        = source->growth_delta0;
   res->cmp_sigma8           = source->cmp_sigma8;
   res->transfer_s           = source->transfer_s;
   res->transfer_alpha_Gamma = source->transfer_alpha_Gamma;

   return res;
}

void read_cosmological_parameters(cosmo **self, FILE *F, error **err)
{
   cosmo *tmp;
   struct {
     char snonlinear[CSLENS], stransfer[CSLENS], sgrowth[CSLENS], sde_param[CSLENS];
   } strings;
   config_element c = {0, 0.0, ""};
   int j;
   char str[128];

   tmp = set_cosmological_parameters_to_default(err);
   forwardError(*err, __LINE__,);

   CONFIG_READ(tmp, Omega_m, d, F, c, err);
   CONFIG_READ(tmp, Omega_de, d, F, c, err);
   CONFIG_READ(tmp, w0_de, d, F, c, err);
   CONFIG_READ(tmp, w1_de, d, F, c, err);
   CONFIG_READ(tmp, h_100, d, F, c, err);
   CONFIG_READ(tmp, Omega_b, d, F, c, err);
   CONFIG_READ(tmp, Omega_nu_mass, d, F, c, err);
   CONFIG_READ(tmp, Neff_nu_mass, d, F, c, err);
   CONFIG_READ(tmp, normalization, d, F, c, err);
   CONFIG_READ(tmp, n_spec, d, F, c, err);

   CONFIG_READ_S(&strings, snonlinear, s, F, c, err);
   STRING2ENUM(tmp->nonlinear, strings.snonlinear, nonlinear_t, snonlinear_t, j, Nnonlinear_t, err);

   CONFIG_READ_S(&strings, stransfer, s, F, c, err);
   STRING2ENUM(tmp->transfer, strings.stransfer, transfer_t, stransfer_t, j, Ntransfer_t, err);

   CONFIG_READ_S(&strings, sgrowth, s, F, c, err);
   STRING2ENUM(tmp->growth, strings.sgrowth, growth_t, sgrowth_t, j, Ngrowth_t, err);

   CONFIG_READ_S(&strings, sde_param, s, F, c, err);
   STRING2ENUM(tmp->de_param, strings.sde_param, de_param_t, sde_param_t, j, Nde_param_t, err);

   if (tmp->de_param == poly_DE) {
      CONFIG_READ(tmp, N_poly_de, i, F, c, err);
      tmp->w_poly_de = malloc_err(sizeof(double) * tmp->N_poly_de, err);
      forwardError(*err, __LINE__,);
      CONFIG_READ_ARR(tmp, w_poly_de, d, j, tmp->N_poly_de, str, F, c, err);
   }

   CONFIG_READ(tmp, normmode, i, F, c, err);

   CONFIG_READ(tmp, a_min, d, F, c, err);

   *self = copy_parameters_only(tmp, err);
   forwardError(*err, __LINE__,);
}

cosmo *set_cosmological_parameters_to_default(error **err)
{
   /* Parameters are:
      Om Od w0 w1 h Ob Onu Neffnu s8 ns
      nonlin transfer growth deparam norm amin
   */

   cosmo *self;

   self = init_parameters(0.25, 0.75, -1.0, 0.0, NULL, 0, 0.70, 0.044, 0.0, 0.0, 0.80, 0.96,
			  smith03, eisenhu, growth_de, linder, norm_s8, 0.0, err);
   forwardError(*err, __LINE__, NULL);

   return self;
}

cosmo *set_cosmological_parameters_to_default2(error **err)
{
   /* Parameters are:
      Om Od w0 w1 h Ob Onu Neffnu s8 ns
      nonlin transfer growth deparam norm amin
   */

   cosmo *self;

   self = init_parameters(0.25, 0.75, -1.0, 0.0, NULL, 0, 0.7, 0.04307, 0.0, 0.0, 0.80, 0.96,
			  smith03, eisenhu, growth_de, linder, norm_s8, 1.0/1211.0, err);
   forwardError(*err, __LINE__, NULL);

   return self;
}

void free_parameters(cosmo** self)
{
   cosmo *s;

   s = *self;

   if (s->de_param == poly_DE) {
      free(s->w_poly_de);
   }

   del_interTable(&s->linearGrowth);
   del_interTable(&s->transferFct);
   del_interTable(&s->transferBE);
   del_interTable2D(&s->P_NL);
   del_interTable(&s->slope);
   del_interTable(&s->w);
   //del_interTable(&s->k_NL);
   //if (&s->ystar_allz != NULL) free(&s->ystar_allz);

   free(s);
   s = NULL;
}

void updateParameters(cosmo* model, double OMEGAM, double OMEGADE, double W0_DE, double W1_DE,
		      double *W_POLY_DE, int N_POLY_DE,
		      double H100, double OMEGAB, double OMEGANUMASS,
		      double NEFFNUMASS, double NORM, double NSPEC,
		      nonlinear_t NONLINEAR, transfer_t TRANSFER, growth_t GROWTH,
		      de_param_t DEPARAM, norm_t normmode, double AMIN, error **err)
{
   cosmo* prev;

   prev = copy_parameters_only(model, err);
   forwardError(*err,__LINE__,);

   model->Omega_m   = OMEGAM;
   model->Omega_de  = OMEGADE;
   model->w0_de     = W0_DE;
   model->w1_de     = W1_DE;
   set_w_poly_de(&model->w_poly_de, &model->N_poly_de, W_POLY_DE, N_POLY_DE, 0, err);
   forwardError(*err, __LINE__,);

   model->h_100     = H100;
   model->Omega_b   = OMEGAB;
   model->Omega_nu_mass = OMEGANUMASS;
   model->Neff_nu_mass  = NEFFNUMASS;
   model->normalization = NORM;
   model->n_spec        = NSPEC;

   model->nonlinear = NONLINEAR;
   model->transfer  = TRANSFER;
   model->growth    = GROWTH;
   model->de_param  = DEPARAM;
   model->normmode  = normmode;

   if (AMIN>0.0) model->a_min     = AMIN;
   else model->a_min = a_minmin;

   updateFrom(prev, model, err); forwardError(*err,__LINE__,);
   set_norm(model, err);         forwardError(*err,__LINE__,);
   
   free_parameters(&prev);
}

/* ============================================================ *
 * Deletes pre-calculated tables in cosmo structure apres       *
 * if corresponding parameters are different compared to avant, *
 * so tables will be re-calculated when needed.                 *
 * To be used after parameters have been changed, e.g. with     *
 * updateParameters, or manually.                               *
 * ============================================================ */

void updateFrom(cosmo* avant, cosmo* apres, error **err)
{
   consistency_parameters(apres, err);    forwardError(*err, __LINE__,);

   if (change_D_plus(avant,apres)) {
      del_interTable(&(apres->linearGrowth));
      apres->growth_delta0=-1;
   }
   if (change_Delta_L_BE2(avant,apres)) {
      del_interTable(&(apres->transferBE));
   }
   if (change_Tsqr(avant,apres)) {
      del_interTable(&(apres->transferFct));
      apres->transfer_alpha_Gamma = -1;
      apres->transfer_s = -1;
   }
   if (change_sigma_8_sqr(avant,apres)) {
      apres->cmp_sigma8=-1;
   }
   if (change_P_NL(avant,apres)) {
      del_interTable2D(&(apres->P_NL));
      del_interTable(&(apres->slope));

      /* MKDEBUG: A bit of overkill; The coyote13 *
       * power spectrum depends on                *
       * fewer parameters than the general P_NL   */
      if (apres->ystar_allz != NULL) free(apres->ystar_allz);
   }
   if (change_w(avant,apres)) {
      del_interTable(&(apres->w));
   }
   if (change_norm(avant, apres)) {
      set_norm(apres, err);
      forwardError(*err,__LINE__,);
   }
}

int change_norm(cosmo* avant, cosmo* apres) {
   if (NCOEQ(avant,apres,normmode) || NCOCLOSE(avant,apres,normalization))
     return 1;
   if (change_D_plus(avant,apres) || change_sigma_8_sqr(avant,apres)) return 1;
   return 0;
}

void set_norm(cosmo* self, error **err)
{
   double d0, nrm;

   testErrorRet(self->normmode==norm_as, ce_badFormat,
		"norm_as at the moment not supported", *err, __LINE__,);
   self->sigma_8 = self->normalization;
   self->As      = 0.0;
   return;
   
   
   d0  = D_plus(self, 1.0, 0, err); forwardError(*err,__LINE__,);
   nrm = sigma_8_sqr(self, err); forwardError(*err,__LINE__,);

   if (self->normmode==norm_s8) {
     self->sigma_8 = self->normalization;
     self->As      = self->sigma_8*self->sigma_8/d0/d0/nrm;
     return;
   }
   if (self->normmode==norm_as) {
     self->As      = self->normalization;
     self->sigma_8 = sqrt(self->As * d0*d0*nrm);
     return;
   }

   *err = addError(ce_unknown, "unknown normmode", *err, __LINE__);
}

void set_w_poly_de(double **w_target, int *N_target, const double *w_source, int N_source, int check, error **err)
{
   int i;

   if (check) {
      testErrorRet(w_source==NULL, ce_de, "DE parameter array is NULL", *err, __LINE__,);
      testErrorRetVA(N_source<=0, ce_de, "Number of DE parameters (%d) not positive", *err, __LINE__,, N_source);
   }

   *w_target = malloc_err(sizeof(double) * N_source, err);
   forwardError(*err, __LINE__,);
   for (i=0; i<N_source; i++) (*w_target)[i] = w_source[i];
   *N_target = N_source;
}

void dump_param(cosmo* self, FILE *F)
{
   int i;

   if (!F) F = stderr;
   fprintf(F, "#  O_m    O_de   w0_de  w1_de  h_100  O_b    O_nu   Neffnu sig_8  n_s "
	   "nonlinear transfer growth de_param normmode N_a\n");
   fprintf(F, "# % .3f % .3f % .3f % .3f % .3f % .3f % .3f % .3f % .3f % .3f "
	   "%s(%d) %s(%d) %s(%d) %s(%d) %d %d\n",
	   self->Omega_m, self->Omega_de, self->w0_de, self->w1_de, self->h_100, self->Omega_b,
	   self->Omega_nu_mass, self->Neff_nu_mass, self->sigma_8, self->n_spec,
	   snonlinear_t(self->nonlinear), self->nonlinear, stransfer_t(self->transfer), self->transfer,
	   sgrowth_t(self->growth), self->growth, sde_param_t(self->de_param), self->de_param,
	   self->normmode, self->N_a);

   if (self->de_param == poly_DE) {
      fprintf(F, "# Npde wpde_i\n");
      fprintf(F, "    %d", self->N_poly_de);
      for (i=0; i<self->N_poly_de; i++) {
	 fprintf(F, " %.3f ", self->w_poly_de[i]);
      }
      fprintf(F, "\n");
   }
}

void dump_param2(cosmo* self, FILE *F, char* pream)
{
   if (!F) F = stderr;
   fprintf(F,"%s:: \n", pream);
   dump_param(self, F);
}


/* ============================================================ *
 * Cosmology.							*
 * ============================================================ */

/* CPT 9, no dark energy, no neutrinos */
double da_dtau(cosmo* self,double a, error **err)
{
   double res, det;

   testErrorRet(a<EPSILON, ce_underflow, "Division by 0, the scale factor a is too small",
		*err, __LINE__, 0.0);
   det = 1.0 + (self->Omega_m+self->Omega_nu_mass)*(1.0/a - 1.0);
   forwardError(*err,__LINE__,0);
   det	+= self->Omega_de*(dsqr(a) - 1.0);
   forwardError(*err,__LINE__,0);
   res = sqrt(det);
   return res;
}

double da_dtau_m3(double a, void *intpar, error **err)
{
   double res;
   if (a < EPSILON) return 0.0;
   res = da_dtau((cosmo*) intpar,a,err);
   forwardError(*err,__LINE__,0);
   testErrorRet(res<EPSILON, ce_underflow, "Division by 0", *err, __LINE__, 0.0);
   return 1.0/(res*res*res);
}

/* Wetterich 2004, Phys. Lett. B, 594, 17 */
#define Omegadeinf0 1.0e-5
double b_early(double w0, double Omegam, double Omegadeinf, error **err)
{
   double b;

   testErrorRet(Omegadeinf<=0 || Omegadeinf>=1, ce_negative, "Omegadeinf (w1_de) outside [0;1]",
		*err, __LINE__, 0.0);

   if (Omegadeinf<Omegadeinf0) {
      b = -3.0*w0/(-log(Omegadeinf) + log((1.0-Omegam)/Omegam));
   } else {
      b = -3.0*w0/(log((1.0-Omegadeinf)/Omegadeinf) + log((1.0-Omegam)/Omegam));
   }

   return b;
}

double w_de(cosmo* self, double a, error **err)
{
   double res, b;
   int i;

   res = self->w0_de;

   switch (self->de_param) {
      case jassal  :
	 if (fabs(self->w1_de)>EPSILON) res += self->w1_de*a*(1.0-a);
	 break;

      case linder  :
	 if (fabs(self->w1_de)>EPSILON) res += self->w1_de*(1.0-a);
	 break;

      case earlyDE :
	 /* TODO: neutrinos? */
	 b = b_early(self->w0_de, self->Omega_m, self->w1_de, err);
	 forwardError(*err, __LINE__, 0.0);
	 testErrorRet(fabs(1.0-b*log(a))<EPSILON, ce_singularValue, "Division by zero",
		      *err, __LINE__, 0.0);
	 res = self->w0_de/dsqr(1.0-b*log(a));
	 break;

      case poly_DE :
	 for (i=0,res=0.0,b=1.0; i<self->N_poly_de; i++) {
	    res += self->w_poly_de[i] * b;
	    b   *= a;   /* b = a^i */
	 }
	 break;

      default     : 
	 *err = addError(ce_unknown, "Unknown de_param value", *err, __LINE__);
	 return 0;
	 break;

   }

   return res;
}

/* Per05 (4) */
double f_de(cosmo* self, double a, error **err)
{
   double res, b;
   int i;

   switch (self->de_param) {

      case jassal :
	 res = -3.0*(1.0 + self->w0_de);
	 if (fabs(self->w1_de)>EPSILON)
	   res += 3.0*self->w1_de/(2.0*log(a))*(1.0-a)*(1.0-a);
	 break;

      case linder :
	 res = -3.0*(1.0 + self->w0_de);
	 if (fabs(self->w1_de)>EPSILON)
	   res += 3.0*self->w1_de*((a-1.0)/log(a) - 1.0);
	 break;

      case earlyDE :
	 /* TODO: neutrinos? */
	 b   = b_early(self->w0_de, self->Omega_m, self->w1_de, err);
	 forwardError(*err, __LINE__, 0.0);
	 testErrorRet(fabs(1.0+b*log(a))<EPSILON, ce_singularValue, "Division by zero",
		      *err, __LINE__, 0.0);
	 res = -3.0*(1.0 + self->w0_de/(1.0-b*log(a)));
	 break;

      case poly_DE :
	 res = -3.0 * (1.0 + self->w_poly_de[0]);
	 for (i=1; i<self->N_poly_de; i++) {
	    res += -3.0 / log(a) * self->w_poly_de[i] * (pow(a, i) - 1.0) / (double)i;
	 }
	 break;

      default     : 
	 *err = addErrorVA(ce_unknown, "Unknown de_param value %d", *err, __LINE__, self->de_param);
	 return 0.0;

   }

   return res;
}

int change_Esqr(cosmo *avant, cosmo *apres)
{
   if (NCOEQ(avant, apres, Omega_m) || NCOEQ(avant, apres, Omega_de) || NCOEQ(avant, apres, h_100) ||
       NCOEQ(avant, apres, Omega_nu_mass)) return 1;
   if (change_w_de(avant, apres)) return 1;
   return 0;
}

/* ============================================================ *
 * Returns E^2(a) = [H(a)/H_0]^2, see Percival (2005).		*
 * If wOmegar=1, Omega_radiation>0 is included (photons +       *
 * neutrinos), needed for high-z quantities such as the sound   *
 * horizon at the drag epoch. Note: For low redshift, wOmega=0  *
 * should be used, otherwise the nonlinear power-spectrum       *
 * fitting formulae might not work.				*
 * ============================================================ */

double Esqr(cosmo* self, double a, int wOmegar, error **err)
{
   double asqr, EE, Omega_r;

   asqr = a*a;

   if (wOmegar==1) {
      /* For early-Universe calculations like r_sound(z_*) */
      Omega_r = omega_gamma/self->h_100/self->h_100*(1.0 + 0.2271*NEFF);
   } else {
      /* Omega_r produces too much power for PD96 !!?!! */
      Omega_r = 0.0;
   }


   EE = self->Omega_m/(asqr*a)
     + (1.0 - self->Omega_m - self->Omega_de - self->Omega_nu_mass - Omega_r)/asqr
     + self->Omega_de*pow(a, f_de(self, a, err))
     + Omega_r/(asqr*asqr);
   forwardError(*err,__LINE__,0.0);

   if (self->Omega_nu_mass>EPSILON) EE += self->Omega_nu_mass/pow(a, 3.0*(1+w_nu_mass(self, a)));

   testErrorRetVA(EE<0, ce_negative, "E^2 = (H/H_0)^2 is negative for a=%g", *err, __LINE__, 0.0, a);
   testErrorRetVA(!finite(EE), ce_infnan, "E^2 not finite for a=%g", *err, __LINE__, 0.0, a);

   return EE;
}

/* Per05 (6) */
double Omega_m_a(cosmo* self, double a, double Esqrpre, error **err)
{
   double EE;
   if (Esqrpre>0) EE = Esqrpre;
   else {
      EE = Esqr(self, a, 0, err);
      forwardError(*err,__LINE__,0);
   }
   return self->Omega_m/(a*a*a)/EE;
}

/* Per05 (6) */
double Omega_de_a(cosmo* self,double a, double Esqrpre, error **err)
{
   double EE,res;

   if (Esqrpre>0) EE = Esqrpre;
   else {
      EE = Esqr(self, a, 0, err);
      forwardError(*err, __LINE__, 0);
   }
   res = self->Omega_de*pow(a, f_de(self, a, err))/EE;
   forwardError(*err,__LINE__,0);

   return res;
}

double w_nu_mass(cosmo *self, double a)
{
   double m_nu, w, m0, alf, bet;

   m_nu = 93.0*dsqr(self->h_100)*self->Omega_nu_mass/3.0;

   m0  = 0.000585;
   alf = 1.652;
   bet = 0.561;
   w = pow(1.0+pow(m_nu/m0*a, alf), -bet)/3.0;

   /* 1/3 ok for m_nu <= 0.05 eV */
   //return (1+0.618232*pow(m_nu*a,0.753005))/(1+9.58117*pow(m_nu*a,1.34954));
   return w;
}

void D_plus_derivs(double a, double *y, double *yp, void* extra, error **err)
{
   double q, r, om, ode, EE;
   cosmo* self;

   self = (cosmo*) extra;
   testErrorRet(a==0.0, ce_underflow,
		"The scale factor a has to be larger than zero. Check a_min (struct cosmo)",
		*err, __LINE__,);

   EE  = Esqr(self, a, 0, err);          forwardError(*err,__LINE__,);
   om  = Omega_m_a(self, a, EE, err);    forwardError(*err,__LINE__,);
   ode = Omega_de_a(self, a, EE, err);   forwardError(*err,__LINE__,);
   /* TODO: neutrinos */

   q = (2.0 - 0.5*(om + (1.0+3.0*w_de(self,a,err))*ode))/a;
   forwardError(*err,__LINE__,);
   r = 1.5*om/a/a;

   yp[1] = y[2];
   yp[2] = -q*y[2] + r*y[1];
}

int change_D_plus(cosmo* avant, cosmo* apres)
{
   if (NCOEQ(avant,apres,growth) || NCOEQ(avant,apres,N_a))
     return 1;
   switch (avant->growth) {
      case heath :
	 if (NCOCLOSE(avant,apres,Omega_m) || NCOCLOSE(avant,apres,Omega_de) ||
	     NCOCLOSE(avant,apres,a_min)) return 1;
      case growth_de :
         /* TODO: Neff_nu?? */
	 if (change_w(avant,apres)) return 1;
   }
   return 0;
}

#define nvar 2
#define eps  1.0e-8
double D_plus(cosmo* self, double a, int normalised, error **err)
{

   double delta, aa, res;
   int  i;
   double *table;
   double ystart[nvar+1], a1, a2, h1, hmin, a1min;
   int nok, nbad;
   double da;
   interTable *linearGrowth;

   if (self->linearGrowth == NULL) {
      da = (1.0 - self->a_min)/(self->N_a-1.0);
      linearGrowth = init_interTable(self->N_a,self->a_min,1.,da,0.0,0.0,err);
      forwardError(*err,__LINE__,0);
      table=linearGrowth->table;
      switch (self->growth) {

	 case heath :
	    self->growth_delta0 = 2.5*self->Omega_m
	      *sm2_qromberg(da_dtau_m3, (void*)self, 0.0, 1.0, 1.0e-6, err);
	    forwardError(*err,__LINE__,0);
	    for (i=0,aa=self->a_min; i<self->N_a; i++, aa+=da) {
	       delta = 2.5*self->Omega_m
		 *sm2_qromberg(da_dtau_m3, (void*)self, 0.0, aa, 1.0e-6, err);
	       forwardError(*err,__LINE__,0);
	       table[i] = da_dtau(self,aa,err)/aa*delta;
	       forwardError(*err,__LINE__,0);
	    }
	    break;
				
	 case growth_de :
	    a1min = 0.0001;
	    if (self->de_param!=earlyDE) {
	       ystart[1] = 5.0/3.0*a1min;   /* D (EdS)  */
	    } else {
	       ystart[1] = 5.0/3.0*pow(a1min, 1.0-3.0*Omega_de_a(self, a1min, -1, err)/5.0);
	       forwardError(*err, __LINE__, 0.0);
	    }
	    ystart[2] = 0.0;             /* D'(EdS)  */
	    h1   = 0.0001;
	    hmin = 0.0;

	    rkdrive_var_t bsvar;
	    init_rkdrive_var(&bsvar);
	    odeint(ystart, nvar, a1min, self->a_min, eps, h1, (void*) self, hmin, &nok, &nbad,
		    D_plus_derivs, bsstep, &bsvar, err);
	    forwardError(*err,__LINE__,0);
					
	    table[0] = ystart[1];

	    for (i=1,a1=self->a_min,a2=a1+da; i<self->N_a; i++,a1+=da,a2+=da) {
	       odeint(ystart, nvar, a1, a2, eps, h1, (void*) self, hmin, &nok, &nbad,
		       D_plus_derivs, bsstep, &bsvar, err);
	       forwardError(*err,__LINE__,0);
	       table[i] = ystart[1];
	    }
	    self->growth_delta0 = table[self->N_a-1];
	    break;

	 default : 
	    *err = addError(ce_unknown,"Unknown growth",*err,__LINE__);
	    return 0;
	    break;
      }
      self->linearGrowth = linearGrowth;
   }

   res = interpol_wr(self->linearGrowth, a, err);
   forwardError(*err,__LINE__,0);

   if (normalised==1) return res/self->growth_delta0;
   else return res;
}
#undef nvar
#undef eps

/* PD 96 - (15,16) */
void Omega_a(cosmo* self, double a, double *omega_m, double *omega_v)
{
   double f, a3;
   a3 = a*a*a;
   f = a + self->Omega_m*(1.-a) + self->Omega_de*(a3-a);
   *omega_m = self->Omega_m/f;
   *omega_v = self->Omega_de*a3/f;
}

/* CPT 29, PD 15+16, actually: Lahav 1991, Lightman&Schechter 1990 */
double g(cosmo* self, double a)
{
   double omega_m, omega_v;
   Omega_a(self,a, &omega_m, &omega_v);

   return 2.5*omega_m/(pow(omega_m,4.0/7.0) - omega_v +
		       (1.0 + 0.5*omega_m)*(1.0 + omega_v/70.0));
}


/* ============================================================ *
 * Transfer function. k in h/Mpc.				*
 * ============================================================ */
int change_Tsqr(cosmo* avant, cosmo* apres)
{
   if (NCOEQ(avant,apres,transfer))
      return 1;
   switch (apres->transfer) {
      case bbks : case eisenhu : case eisenhu_osc : case be84 :
         if (NCOCLOSE(avant,apres,Omega_m) || NCOCLOSE(avant,apres,h_100) || NCOCLOSE(avant,apres,Omega_b)
               || NCOCLOSE(avant,apres,n_spec) || NCOCLOSE(avant,apres,Omega_nu_mass) || NCOCLOSE(avant,apres,Neff_nu_mass))
            return 1;
      default :
            return 0;
         /* TODO: neutrinos, T_nu from BBKS */
	 return 0;
   }
   return 0;
}

/* EH98 (19). [k] = Mpc/h */
double T_tilde(const cosmo *self, double k, double alpha_c, double beta_c)
{
   double q, T0, L, C;

   /* EH98 (10); [q] = [k] = h/Mpc */
   q  = k * dsqr(T_CMB/2.7) / (self->Omega_m * self->h_100);
   L  = log(M_E + 1.8 * beta_c * q);
   C  = 14.2 / alpha_c + 386.0 / (1.0 + 69.9 * pow(q, 1.08));
   T0 = L/(L + C*q*q);

   return T0;
}

/* ============================================================ *
 * Returns sqaure of transfer function times k^n, [k] = h/Mpc.	*
 * ============================================================ */
double Tsqr_one(cosmo *self, double k, double Gamma_eff, error **err)
{
   double f1, f2, q, res, L, C;
   double Tc, Tb, f, a1, a2, om, fb, fc,b1, b2, alpha_c, beta_c, T_2_7_sqr, z_eq, z_d,
     beta_node, beta_b, alpha_b, tilde_s, R_d, k_eq, s;

   switch (self->transfer) {

      case bbks :
 
	 q = k/Gamma_eff;
         f1 = log(1 + 2.34*q)/(2.34*q);
         f2 = 1 + q*(3.89 + q*(259.21 + q*(162.771336 + q*2027.16958081)));
         res = dsqr(f1)/sqrt(f2)*pow(k, self->n_spec);
         break;

      case eisenhu : 

	 /* [q] = [kk] = h/Mpc. EH98 (28) */
	 q = k * dsqr(T_CMB/2.7)/Gamma_eff;

	 /* W.Hu's problem set 4 */
	 /* L = log(2.71828 + 1.84*q*alpha_Gamma); */
	 /* C = 14.4 + 325./(1. + 60.5*pow(q, 1.11)); */

	 /* EH98 (29) */
	 L = log(2.*2.71828 + 1.8*q);
	 C = 14.2 + 731.0/(1.0 + 62.5*q);

	 res = dsqr(L/(L + C*q*q))*pow(k, self->n_spec);
	 break;

      case eisenhu_osc :

	 om      = self->Omega_m * self->h_100 * self->h_100;
	 fb      = self->Omega_b / self->Omega_m;
	 fc      = (self->Omega_m - self->Omega_b) / self->Omega_m;

	 if (self->transfer_s<0) {
	    self->transfer_s = r_sound_drag_analytical(self, err);
	    forwardError(*err, __LINE__, 0.0);
	 }
	 s = self->transfer_s;

	 /* Cold dark matter transfer function */

	 /* EH98 (11, 12) */
	 a1      = pow(46.9*om, 0.670) * (1.0 + pow(32.1*om, -0.532));
	 a2      = pow(12.0*om, 0.424) * (1.0 + pow(45.0*om, -0.582));
	 alpha_c = pow(a1, -fb) * pow(a2, -DCUB(fb));
	 b1      = 0.944 / (1.0 + pow(458.0*om, -0.708));
	 b2      = pow(0.395*om, -0.0266);
	 beta_c  = 1.0 + b1*(pow(fc, b2) - 1.0);
	 beta_c  = 1.0 / beta_c;

	 /* EH98 (17, 18) */
	 f  = 1.0 / (1.0 + dsqr(dsqr(k * s / 5.4)));
	 Tc = f * T_tilde(self, k, 1, beta_c) + (1.0 - f) * T_tilde(self, k, alpha_c, beta_c);

	 /* Baryon transfer function */

	 /* EH98 (14, 21) */
	 z_d       = z_drag(self);
	 R_d       = ratio_b_gamma(self, 1.0 / (z_d + 1.0));
	 T_2_7_sqr = dsqr(T_CMB/2.7);
	 k_eq      = 7.46e-2 * om / T_2_7_sqr / self->h_100;   /* [h/Mpc] */
	 z_eq      = 2.50e4 * om / T_2_7_sqr / T_2_7_sqr;
	 alpha_b   = 2.07 * k_eq * s * pow(1.0 + R_d, -0.75)
	              * G_EH98((1.0 + z_eq) / (1.0 + z_d));

	 beta_node = 8.41 * pow(om, 0.435);
	 tilde_s   = s / pow(1.0 + DCUB(beta_node / (k * s)), 1.0/3.0);
	 beta_b    = 0.5 + fb + (3.0 - 2.0 * fb) * sqrt(dsqr(17.2 * om) + 1.0);
	 /* [tilde_s] = Mpc/h */

	 Tb = (T_tilde(self, k, 1.0, 1.0) / (1.0 + dsqr(k * s / 5.2))
	       + alpha_b / (1.0 + DCUB(beta_b/(k * s))) * exp(-pow(k / k_silk(self), 1.4)) )
	   * sinc(k*tilde_s);

	 /* Total transfer function */
	 res = dsqr(fb * Tb + fc * Tc) * pow(k, self->n_spec);
	 break;

	 //case be84 :
	 //Delta_L = Delta_L_BE2(self, k, err);   forwardError(*err, __LINE__, 0.0);
	 //break;

      default :
	 *err = addError(ce_unknown, "Unknown transfer function", *err, __LINE__);
	 return 0; 
   }

   return res;
}

double Gamma_Sugiyama(cosmo *self)
{
   double Gamma;
   Gamma = self->h_100*self->Omega_m*exp(-self->Omega_b - sqrt(self->h_100/0.5)*self->Omega_b/self->Omega_m);
   return Gamma;
}

double int_for_r_sound(double a, void *intpar, error **err)
{
   cosmo *self;
   double E, res, h2, R;

   if (a==0) return 0.0;

   self = (cosmo*)intpar;
   E    = sqrt(Esqr(self, a, 1, err));
   forwardError(*err, __LINE__, -1.0);

   h2  = self->h_100*self->h_100;
   R   = a*0.75*self->Omega_b/(omega_gamma/h2);
   res = a*a*E*sqrt(1.0 + R);

   if (res<=0) return 0.0;

   return 1.0/res;

}

/* Comoving sound horizon [Mpc/h], K08 */
double r_sound_integral(cosmo *self, double a, error **err)
{
   double rs;

   rs = sm2_qromberg(int_for_r_sound, (void*)self, 0.0, a, 1.0e-7, err);
   forwardError(*err, __LINE__, -1.0);

   return rs*R_HUBBLE/sqrt(3.0);
}

/* Eisenstein & Hu (1998), eq. 26. Sound horizon at drag epoch, [Mpc/h] */
double r_sound_drag_fit(cosmo *model, error **err)
{
   double om, ob, h2, r_s;

   h2   = model->h_100 * model->h_100;
   om   = model->Omega_m*h2;
   ob   = model->Omega_b*h2;

   testErrorRetVA(ob < 0.0125, ce_range, "Omega_b h^2 = %g has to be larger than %g (see Eisenstein & Hu 1998, eq. (26)",
		  *err, __LINE__, 0.0, ob, 0.0125);
   testErrorRetVA(om < 0.025, ce_range, "Omega_m h^2 = %g has to be larger than %g (see Eisenstein & Hu 1998, eq. (26)",
		  *err, __LINE__, 0.0, om, 0.025);
   testErrorRetVA(om > 0.5, ce_range, "Omega_m h^2 = %g has to be larger than %g (see Eisenstein & Hu 1998, eq. (26)",
		  *err, __LINE__, 0.0, om, 0.5);

   r_s  = 44.5*log(9.83/om) / sqrt(1.0 + 10.0*pow(ob, 0.75));

   return r_s * model->h_100;
}

/* Eisenstein & Hu (1998) eq. (6). Sound horizon at drag epoch, [Mpc/h] */
double r_sound_drag_analytical(cosmo *self, error **err)
{
   double R_d, R_eq, k_eq, z_eq, z_d, omega_m, T_2_7_sqr, s;

   T_2_7_sqr = dsqr(T_CMB/2.7);
   omega_m   = self->Omega_m*self->h_100*self->h_100;
   k_eq      = 7.46e-2*omega_m/T_2_7_sqr;              /* [1/Mpc] */
   z_eq      = 2.50e4*omega_m/T_2_7_sqr/T_2_7_sqr;

   z_d       = z_drag(self);
   R_d       = ratio_b_gamma(self, 1.0/(z_d+1.0));
   R_eq      = ratio_b_gamma(self, 1.0/(z_eq+1.0));

   testErrorRetVA(R_eq<EPSILON1, ce_infnan,
		  "Division by zero: Probably Omega_b=%g is too small to determine the sound horizon at drag redshift",
		  *err, __LINE__, 0.0, self->Omega_b);

   s = 2.0/(3.0*k_eq) * sqrt(6.0/R_eq) * log( (sqrt(1.0+R_d) + sqrt(R_eq+R_d))/(1.0 + sqrt(R_eq)) );

   return s * self->h_100;
}

/* Silk damping scale, in h/Mpc, fit from EH98 (7) */
double k_silk(const cosmo *model)
{
   double om, ob, k_s;

   om  = model->Omega_m * model->h_100 * model->h_100;
   ob  = model->Omega_b * model->h_100 * model->h_100;
   k_s = 1.6 * pow(ob, 0.52) * pow(om, 0.73) * (1.0 + pow(10.4*om, -0.95));

   return k_s / model->h_100;
}

/* Eisenstein & Hu (1998). Baryon-to-photon ratio */
double ratio_b_gamma(cosmo *self, double a)
{
   double R, z, T_2_7_sqr, h2;

   h2  = self->h_100*self->h_100;
   R   = a*0.75*self->Omega_b/(omega_gamma/h2);
   return R;


   T_2_7_sqr = dsqr(T_CMB/2.7);
   R  = 31.5*self->Omega_b*h2;
   R /= T_2_7_sqr/T_2_7_sqr;

   z  = 1.0/a - 1.0;
   R /= z/1000.0;

   return R;
}

/* Drag epoch, Eisenstein & Hu (1998) */
double z_drag(cosmo *self)
{
   double zd, om, ob, h2, b1, b2, fnu;

   h2 = self->h_100*self->h_100;
   ob = self->Omega_b*h2;

   /* Neutrino fraction, e.g. cosmomc:bao.f90. Neutrinos are relativistic at the drag epoch. *
    * In Percival (2007,2009), fnu=0 is used. */
   //fnu = (21.0/8.0)*pow(4.0/11.0, 4.0/3.0);
   fnu = 0.0;

   om = self->Omega_m*h2  - fnu * (self->Omega_m - self->Omega_b)*h2;

   b1 = 0.313*pow(om, -0.419)*(1.0 + 0.607*pow(om, 0.674));
   b2 = 0.238*pow(om, 0.223);
   zd  = 1291.0*pow(om, 0.251)/(1.0 + 0.659*pow(om, 0.828))
     *(1.0 + b1*pow(ob, b2));

   return zd;
}

/* Eisenstein & Hu (1998) eq. (15) */
double G_EH98(double y)
{
   double G, x;

   x = sqrt(1.0 + y);
   G = y * (-6.0 * x + (2.0 + 3.0*y) * log((x + 1.0) / (x - 1.0)));

   return G;
}

/* ============================================================ *
 * Returns T^2(k) * k^n, [k] = h/Mpc.				*
 * ============================================================ */
#define Nk 500
double Tsqr(cosmo* self, double k, error **err)
{

   double dlogk, logkmin, logkmax;
   double *table;
   interTable *transferFct;
   double kk, f1, Gamma_eff, omhh, f_b;
   int i;

   if (self->transferFct==NULL) {

      logkmin = log(k_min);
      logkmax = log(k_max);
      dlogk   = (logkmax-logkmin)/(Nk-1.0);
      transferFct = init_interTable(Nk,logkmin,logkmax,dlogk,self->n_spec,self->n_spec-4.,err);
      forwardError(*err,__LINE__,0);
      table = transferFct->table;

      testErrorRet(self->h_100<0.0 || self->Omega_m<0.0 || self->Omega_nu_mass<0.0 || self->Omega_b<0.0,
		   ce_negative,
		   "h_100 or density parameter negative", *err, __LINE__, 0.0);

      if (self->transfer==eisenhu || self->transfer==eisenhu_osc) {
	 /* Sound horizon in Mpc/h, EH98 (26) */
	 if (self->transfer_s<0) {
	    self->transfer_s = r_sound_drag_analytical(self, err);
	    forwardError(*err, __LINE__, 0.0);
	 }
      }
      
      switch (self->transfer) {
	 case bbks:
	    /* TODO: neutrinos */
	    Gamma_eff = Gamma_Sugiyama(self);
	    for (i=0,kk=k_min; i<Nk; i++,kk*=exp(dlogk)) {
	       table[i] = log(Tsqr_one(self, kk, Gamma_eff, err));
	       forwardError(*err,__LINE__,0);
	    }
	    break;

	 case eisenhu:
	    /* TODO: neutrinos */
	    omhh = self->Omega_m*self->h_100*self->h_100;
	    f_b  = self->Omega_b/self->Omega_m;

  	    /* EH98 (31) */
	    if (self->transfer_alpha_Gamma<0)
	      self->transfer_alpha_Gamma = 1 - 0.328*log(431*omhh)*f_b 
		+ 0.38*log(22.3*omhh)*f_b*f_b;

	    for (i=0,kk=k_min; i<Nk; i++,kk*=exp(dlogk)) {

	       /* EH98 (30) */
	       Gamma_eff = self->Omega_m*self->h_100*(self->transfer_alpha_Gamma + 
		      (1.0-self->transfer_alpha_Gamma)/(1.0+dsqr(dsqr(0.43*kk*self->transfer_s))));
	       table[i] = log(Tsqr_one(self, kk, Gamma_eff, err));
	       forwardError(*err, __LINE__, 0);

	    }
	    break;

	 case eisenhu_osc :

	    for (i=0,kk=k_min; i<Nk; i++,kk*=exp(dlogk)) {
	       table[i] = log(Tsqr_one(self, kk, -1.0, err));
	       forwardError(*err, __LINE__, 0.0);
	    }
	    break;

	 case be84 :
	    break;

	 default :
	    *err = addErrorVA(ce_unknown, "Unknown transfer type %d", *err, __LINE__, self->transfer);
	    return 0;
      }

      self->transferFct = transferFct;

   }

   if (k>k_max) {

      /* k out of range: calculate transfer function 'by hand' */
      switch (self->transfer) {
	 case bbks :
	    Gamma_eff = Gamma_Sugiyama(self);
	    f1 = Tsqr_one(self, k, Gamma_eff, err);
	    forwardError(*err,__LINE__,0);
	    break;
	 case eisenhu :
	    testErrorRet(self->transfer_alpha_Gamma<0 || self->transfer_s<0,
			 ce_negative, "Transfer function variables not initialised", *err, __LINE__, 0.0);
	    Gamma_eff = self->Omega_m*self->h_100*(self->transfer_alpha_Gamma + 
						   (1.0-self->transfer_alpha_Gamma)/(1.0+dsqr(dsqr(0.43*k*self->transfer_s))));
	    f1 = Tsqr_one(self, k, Gamma_eff, err);
	    forwardError(*err,__LINE__,0);
	    break;
	 case eisenhu_osc :
	    testErrorRet(self->transfer_s<0, ce_negative, "Transfer function variable not initialised", *err, __LINE__, 0.0);
	    f1 = Tsqr_one(self, k, -1.0, err);
	    break;
	 default :
	    *err = addErrorVA(ce_unknown, "Unknown transfer type %d", *err, __LINE__, self->transfer);
	    return 0;
      }

   } else {

      /* Interpolate transfer function */
      f1 = interpol_wr(self->transferFct, log(k), err);
      forwardError(*err,__LINE__,0);
      f1 = exp(f1);

   }

   return f1;
}
#undef Nk

/* ============================================================ *
 * dfridr.c                                                     *
 * NR page 188. Returns derivate of func at x, initial step is  *
 * h. Error estimate in err.                                    *
 * Modified! func depends on two double! (like P_L)             *
 * ============================================================ */

#define CON 1.4
#define CON2 (CON*CON)
#define BIG 1.0e30
#define NTAB 10
#define SAFE 2.0
double sm2_dfridr(double (*func)(cosmo*,double,double,error **), double x, double h,
                  double *errn, double aa, cosmo* self, error **err)
{
   int i,j;
   double errt,fac,hh,**a,ans;

   ans = 1e30; /* dummy initialization */
   testErrorRet(h==0.0, ce_wrongValue, "h has to be larger than zero", *err, __LINE__, 0.0);
   a=sm2_matrix(1,NTAB,1,NTAB,err);
   forwardError(*err,__LINE__,0);

   hh=h;
   a[1][1]=(*func)(self,aa,x+hh,err);
   forwardError(*err,__LINE__,0);
   a[1][1]=(a[1][1]-(*func)(self,aa,x-hh,err))/(2.0*hh);
   forwardError(*err,__LINE__,0);

   *errn=BIG;
   for (i=2;i<=NTAB;i++) {
      hh /= CON;
      a[1][i]=(*func)(self,aa,x+hh,err);
      forwardError(*err,__LINE__,0);
      a[1][i]=(a[1][i]-(*func)(self,aa,x-hh,err))/(2.0*hh);
      forwardError(*err,__LINE__,0);
      fac=CON2;
      for (j=2;j<=i;j++) {
         a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
         fac=CON2*fac;
         errt=fmax(fabs(a[j][i]-a[j-1][i]),fabs(a[j][i]-a[j-1][i-1]));
         if (errt <= *errn) {
            *errn=errt;
            ans=a[j][i];
         }
      }
      if (fabs(a[i][i]-a[i-1][i-1]) >= SAFE*(*errn)) break;
   }
   sm2_free_matrix(a,1,NTAB,1,NTAB);
   return ans;
}
#undef CON
#undef CON2
#undef BIG
#undef NTAB
#undef SAFE

double W_tophat(double x)
{
   if (x<EPSILON) {
      double xsqr;
      xsqr = x*x;
      return 1.0 - xsqr/10.0 + xsqr*xsqr/280.0;
   }
   return 3./(x*x*x)*(sin(x) - x*cos(x));
}

double int_for_sigma_R(double logk, void *intpar, error **err)
{
   double res, k, x, R, tt;
   cosmoANDdouble *extra;
   cosmo *self;
   
   k     = exp(logk);
   extra = (cosmoANDdouble*)intpar;
   self  = extra->self;
   R     = extra->r;
   x     = W_tophat(k*R)/3.0;

   tt  = P_L_nonorm(self, 1.0, k, err);
   forwardError(*err,__LINE__,0);
   tt /= dsqr(self->sigma_8);   /* Undo sigma_8 from P_L_nonorm */

   res  = k*dsqr(k)*tt*x*x;
   return res;
}

int change_sigma_8_sqr(cosmo* avant, cosmo* apres)
{
   return change_Tsqr(avant,apres);
}

/* PD2 42, for a=1, so D_+=1. Normalization of the power spectrum */
double sigma_8_sqr(cosmo* self, error **err)
{
   double integral;
   cosmoANDdouble cANDd;

   if (self->cmp_sigma8<0) {
      cANDd.self = self;
      cANDd.r    = 8.0;    /* [Mpc/h] */
      integral = sm2_qromberg(int_for_sigma_R, (void*)&cANDd, log(k_min), log(k_max), 1.0e-6, err);
      forwardError(*err,__LINE__,0);
      /* The prefactor is 9/(2*pi*pi). The '9' is from the top-hat window (the '3' from W_tophat is *
       * divided out again in int_for_sigma_8. The '2pi^2' is 4pi/(2pi)^3 */
      self->cmp_sigma8 = 4.5/pi_sqr*integral;
   }
   if (self->cmp_sigma8<0) {
      *err = addError(ce_negative,"sigma8 can't be negative",*err,__LINE__);
      return 0;
   }

   return self->cmp_sigma8;
}

/* ============================================================ *
 * Returns the linear power spectrum.				*
 * ============================================================ */
double P_L(cosmo* self, double a, double k, error **err)
{
   double pl, s;

   pl = P_L_nonorm(self, a, k, err);
   forwardError(*err,__LINE__,0.0);

   s = sigma_8_sqr(self, err);
   forwardError(*err,__LINE__,0.0);

   return pl / s;
}

/* ============================================================ *
 * Returns sigma_8^2 * D+^2(a) * k^n * T^2(k).			*
 * The normalised linear power spectrum is obtained by dividing *
 * by int dk k P_L_nonorm(k) FT[W_tophat](k 8 Mpc/h).           *
 * ============================================================ */
double P_L_nonorm(cosmo* self, double a, double k, error **err)
{
   double d, tt;

   if (k<0) {
      *err = addError(ce_negative,"k Negative !",*err,__LINE__);
      return 0;
   }

   switch (self->transfer) {
      case bbks : case eisenhu : case eisenhu_osc :
	 d    = D_plus(self, a, 1, err);          forwardError(*err,__LINE__,0);
	 tt   = Tsqr(self, k, err);               forwardError(*err,__LINE__,0);
	 break;
      default :
	 *err = addError(ce_transfer, "Wrong transfer type", *err, __LINE__);
	 return 0;
   }

   /* MKDEBUG New: Moved to P_L (this function got renamed to P_L_nonorm */
   //s = sigma_8_sqr(self, err);           forwardError(*err,__LINE__,0);
   //return dsqr(self->sigma_8)*ddtt/s;

   return dsqr(self->sigma_8 * d) * tt;
}

/* PD 22 */
double n_L(cosmo* self, double a, double k, error **err)
{
   double diff, hh, errn, n;

   hh   = k/20.0;
   diff = sm2_dfridr(P_L, 0.5*k, hh, &errn, a, self, err);
   forwardError(*err,__LINE__,0);
   n    = 0.5*k/P_L(self,a, 0.5*k, err)*diff;
   forwardError(*err,__LINE__,0);
   return n;
}

/* PD 21, 23-27 */
double f_NL(cosmo* self, double x, double a, double k, error **err)
{
   double A, B, alpha, beta, V;
   double c, gg, f0, f1, f2, f3, f4;

   c = 1.0 + n_L(self,a,k,err)/3.0;
   forwardError(*err,__LINE__,0);
   testErrorRet(c<0, ce_negative, "Spectral index n_spec too small for PD fitting formula",
		*err, __LINE__, 0);

   A = 0.482/pow(c, 0.947);
   B = 0.226/pow(c, 1.778);
   alpha = 3.31/pow(c, 0.244);
   beta = 0.862/pow(c, 0.287);
   V = 11.55/pow(c, 0.423);

   switch (self->transfer) {
      case bbks : case eisenhu : case eisenhu_osc :
	 gg = D_plus(self, a, 0, err)/a;
	 forwardError(*err,__LINE__,0);
	 break;
      default :
	 *err = addError(ce_transfer, "Wrong transfer type", *err, __LINE__);
	 return 0;
   }
   
   f0 = pow(A*x, alpha);
   f1 = 1.0 + B*beta*x + pow(f0, beta);
   f2 = f0*gg*gg*gg/(V*sqrt(x));
   f3 = 1.0 + pow(f2, beta);
   f4 = x*pow(f1/f3, 1.0/beta);
	 
   if (!finite(f4)) {
      *err = addError(ce_infnan, "inf or nan encountered", *err, __LINE__);
      return 0;
   }
   return f4;
}

double sm2_transfer(cosmo* self, double k, error **err)
{
   double tt;

   /* New: Don't use k_eff with Peacock Gamma any more (caused 2% bias on linear scales wrt CAMB */
   tt  = P_L(self, 1.0, k, err);            forwardError(*err, __LINE__, 0.0);
   tt *= k * k * k / (2 * pi_sqr);

   return tt;
}

/* ============================================================ *
 * Bond&Efstathiou 1984 approximation to Delta_L		*
 * ============================================================ */

#define N_kk 500
/* in units of H_0/c */
#define kmin 1.0e-5
#define kmax 6000.0
int change_Delta_L_BE2(cosmo* avant, cosmo* apres)
{
   if (NCOCLOSE(avant,apres,Omega_m) || NCOCLOSE(avant,apres,h_100) || NCOCLOSE(avant,apres,Omega_b)
       || NCOCLOSE(avant,apres,n_spec) || NCOCLOSE(avant,apres,normalization)
       || NCOCLOSE(avant,apres,normalization))
     return 1;
   return 0;
}

double Delta_L_BE2(cosmo* self, double k, error **err)
{
   double keff, q8=0.0, logkmin, logkmax, dlogk, tk8=0.0, Gamma=0.0;
   double q, tk;
   double klog, kk, res;
   int j;
   double* table; interTable *transferBE;

   if (self->transferBE==NULL) {
      Gamma      = Gamma_Sugiyama(self);
      keff       = 0.172 + 0.011*dsqr(log(Gamma/0.36));
      q8         = 1.0e-20 + keff/Gamma;
      tk8        = 1.0/pow(1+pow(6.4*q8+pow(3.0*q8,1.5)+dsqr(1.7*q8),1.13),1.0/1.13);
      logkmin    = log(kmin);
      logkmax    = log(kmax);
      dlogk      = (logkmax - logkmin)/(N_kk-1.);
      transferBE = init_interTable(N_kk, logkmin, logkmax, dlogk, 0.0, 0.0,err);
      forwardError(*err,__LINE__,0);
		
      table = transferBE->table;
      for (j=0,klog=logkmin; j<N_kk; j++,klog+=dlogk) {
	 kk = exp(klog);
	 q  = 1e-20 + kk/Gamma;
	 tk = 1/pow(1+pow(6.4*q+pow(3.0*q,1.5)+dsqr(1.7*q),1.13),1/1.13);
	 table[j] = (3+self->n_spec)*log(q/q8) + 2*log(self->sigma_8*tk/tk8);
      }
      self->transferBE = transferBE;
   }

   kk = log(k);
   if (kk<self->transferBE->a || kk>self->transferBE->b) {
      q  = 1e-20 + k/Gamma;
      tk = 1/pow(1+pow(6.4*q+pow(3.0*q,1.5)+dsqr(1.7*q),1.13),1/1.13);
      res = self->sigma_8*self->sigma_8*pow(q/q8,3+self->n_spec)*tk*tk/tk8/tk8;
   } else {
      res = exp(interpol_wr(self->transferBE, kk,err));
      forwardError(*err,__LINE__,0);
   }
   return res;
}
#undef N_kk
#undef kmin
#undef kmax

/* ============================================================ *
 * Calculates k_NL, n_eff, n_cur.				*
 * ============================================================ */ 
double int_for_wint2_knl(double logk, void *intpar, error **err)
{
   double krsqr, k, r;
   cosmoAND2double* extra;

   extra = (cosmoAND2double*)intpar;
   r     = extra->r;
   k     = exp(logk);
   krsqr = dsqr(k*r);
   r     = sm2_transfer(extra->self, k, err) * exp(-krsqr);
   forwardError(*err,__LINE__,0);

   return r;
}

double int_for_wint2_neff(double logk, void *intpar, error **err)
{
   double krsqr, k, r;
   cosmoAND2double* extra;

   extra = (cosmoAND2double*) intpar;
   r     = extra->r; 
   k     = exp(logk);
   krsqr = dsqr(k*r);
   r     = sm2_transfer(extra->self, k, err) * 2.0 * krsqr * exp(-krsqr);
   forwardError(*err,__LINE__,0);
   return r;
}

double int_for_wint2_ncur(double logk, void *intpar,error **err)
{
   double krsqr, k, r;
   cosmoAND2double* extra;
	
   extra = (cosmoAND2double*) intpar;
   r     = extra->r;
   k     = exp(logk);
   krsqr = dsqr(k*r);
   r     = sm2_transfer(extra->self, k, err) * 4.0 * krsqr * (1.0-krsqr) * exp(-krsqr);
   forwardError(*err,__LINE__,0);
   return r;
}

#define kmin 1.e-2
#define kmaxdefault 1.e-8
void wint2(cosmo* self, double r, double *sig, double *d1, double *d2, double a, int onlysig,
	   error **err, double precision)
{
   double kmax, logkmin, logkmax, s1, s2, s3, amp;
   cosmoAND2double intpar;

   /* Choose upper integration limit to where filter function drops
    * substantially */
   kmax  = sqrt(10.0*log(10.0))/r;
   if (kmax<kmaxdefault) kmax = kmaxdefault;

   logkmin = log(kmin);
   logkmax = log(kmax);
   intpar.r = r;
   intpar.a = a;                 /* parameter is not needed */
   intpar.self = self;

   switch (self->transfer) {
      case bbks : case eisenhu : case eisenhu_osc :
         amp = D_plus(self, a, 1, err);   forwardError(*err,__LINE__,);
         break;
      default :
         *err = addError(ce_transfer, "Wrong transfer type", *err, __LINE__);
         return;
   }

   if (onlysig==1) {
      s1   = sm2_qromberg(&int_for_wint2_knl, (void*)&intpar, logkmin, logkmax, precision, err);
      forwardError(*err, __LINE__,);
      *sig = amp*sqrt(s1);
   } else s1 = dsqr(1.0/amp);   /* sigma = 1 */

   if (onlysig==0) {
      s2  = sm2_qromberg(int_for_wint2_neff, (void*)&intpar, logkmin, logkmax, 1.0e-6, err);
      forwardError(*err, __LINE__,);
      s3  = sm2_qromberg(int_for_wint2_ncur, (void*)&intpar, logkmin, logkmax, 1.0e-6, err);
      forwardError(*err, __LINE__,);
      *d1 = -s2/s1;
      *d2 = -dsqr(*d1) - s3/s1;
   }
}
#undef kmin
#undef kmaxdefault

/* Slope in the highly nonlinear regime, c.f. Smith et al (2003) eq. (61) */
double slope_NL(double n, double ncur, double om_m, double om_v)
{
   double gam, f1a, f1b, frac, f1;

   gam = 0.86485 + 0.2989*n + 0.1631*ncur;
   if(fabs(1-om_m)>0.01) {
      f1a = pow(om_m,(-0.0732));
      f1b = pow(om_m,(-0.0307));
      frac = om_v/(1.0-om_m);  
      f1 = frac*f1b + (1-frac)*f1a;
   } else {
      f1 = 1.0;
   }

   return 3.0*(f1-1.0) + gam - 3.0;
}

void halofit(double k, double n, double ncur, double knl, double plin, 
	     double om_m, double om_v, double *pnl, nonlinear_t nonlinear, double aa, cosmo *self,
	     error **err)
{
   double gam,a,b,c,xmu,xnu,alpha,beta,f1,f2,f3;
   double y, ysqr;
   double f1a,f2a,f3a,f1b,f2b,f3b,frac,pq,ph;
   double nsqr, w;

   nsqr  = n*n;
   f1b   = pow(om_m, -0.0307);
   f2b   = pow(om_m, -0.0585);
   f3b   = pow(om_m, 0.0743);

   if (nonlinear == smith03_revised) {
      /* Takahashi et al. 2012 */
      
      w = w_de(self, aa, err);
      forwardError(*err, __LINE__,);

      a = 1.5222 + 2.8553*n + 2.3706*nsqr + 0.9903*n*nsqr
         + 0.2250*nsqr*nsqr - 0.6038*ncur + 0.1749*om_v*(1.0 + w);
      a = pow(10.0, a);
      b = pow(10.0, -0.5642 + 0.5864*n + 0.5716*nsqr - 1.5474*ncur + 0.2279*om_v*(1.0 + w));
      c = pow(10.0, 0.3698 + 2.0404*n + 0.8161*nsqr + 0.5869*ncur);
      gam = 0.1971 - 0.0843*n + 0.8460*ncur;
      alpha = fabs(6.0835 + 1.3373*n - 0.1959*nsqr - 5.5274*ncur);
      beta  = 2.0379 - 0.7354*n + 0.3157*nsqr + 1.2490*n*nsqr + 0.3980*nsqr*nsqr - 0.1682*ncur;
      xmu   = 0.0;
      xnu   = pow(10.0, 5.2105 + 3.6902*n);

      f1    = f1b;
      f2    = f2b;
      f3    = f3b;

   } else {

      gam   = 0.86485 + 0.2989*n + 0.1631*ncur;
      a     = 1.4861 + 1.83693*n + 1.67618*nsqr + 0.7940*n*nsqr
         + 0.1670756*nsqr*nsqr - 0.620695*ncur;
      a     = pow(10,a);
      b     = pow(10,(0.9463+0.9466*n+0.3084*nsqr-0.940*ncur));
      c     = pow(10,(-0.2807+0.6669*n+0.3214*nsqr-0.0793*ncur));
      xmu   = pow(10,(-3.54419+0.19086*n));
      xnu   = pow(10,(0.95897+1.2857*n));
      alpha = 1.38848+0.3701*n-0.1452*nsqr;
      beta  = 0.8291+0.9854*n+0.3400*nsqr;

      if(fabs(1-om_m)>0.01) {
         f1a = pow(om_m,(-0.0732));
         f2a = pow(om_m,(-0.1423));
         f3a = pow(om_m,(0.0725));
  
         /* Original halofit */
         frac = om_v/(1.0-om_m);

         if (nonlinear==smith03_de) {

            /* icosmo.org 'Pfusch' */
            double we, wde;
            wde = w_de(self, aa, err);
            forwardError(*err, __LINE__,);
            /* Interpolate between this w and w=-1/3 (behaves like OCDM) */
            we   = frac*wde + (1.0-frac)*(-1.0/3.0);
            frac = -(3.0*we+1.0)/2.0;

         }

         f1 = frac*f1b + (1-frac)*f1a;
         f2 = frac*f2b + (1-frac)*f2a;
         f3 = frac*f3b + (1-frac)*f3a;
      } else {      /* EdS Universe */
         f1 = f2 = f3 = 1.0;
      }

   }

   y = k/knl;
   ysqr = y*y;
   ph = a*pow(y,f1*3)/(1+b*pow(y,f2)+pow(f3*c*y,3-gam));
   ph = ph/(1+xmu/y+xnu/ysqr);
   pq = plin*pow(1+plin,beta)/(1+plin*alpha)*exp(-y/4.0-ysqr/8.0);
   *pnl = pq + ph;

   testErrorRet(!finite(*pnl), ce_overflow, "Power spectrum pnl not finite", *err, __LINE__,);
   testErrorRet(*pnl<0, ce_negative,
         "Negative power spectrum pnl. This can happen if n_eff is very negative and plin large.",
         *err, __LINE__,);
}

double dlog(double x)
{
   return log(x)/log(10.0);
}

int change_w_de(cosmo *avant, cosmo * apres)
{
   int i;

   if (NCOEQ(avant, apres, w0_de) || NCOEQ(avant, apres, w1_de)) return 1;

   if (NCOEQ(avant, apres, N_poly_de)) return 1;
   for (i=0; i<avant->N_poly_de; i++) {
      if (NCOCLOSE(avant, apres, w_poly_de[i])) return 1;
   }
   return 0;
}

int change_P_NL(cosmo* avant, cosmo* apres) {

   if (NCOEQ(avant,apres,growth) || NCOEQ(avant,apres,transfer) || NCOEQ(avant,apres,de_param) || 
       NCOEQ(avant,apres,nonlinear) || NCOEQ(avant,apres,N_a))
     return 1;
   if (NCOCLOSE(avant,apres,Omega_m) || NCOCLOSE(avant,apres,Omega_de) || NCOCLOSE(avant,apres,w0_de) ||
       NCOCLOSE(avant,apres,w1_de) || NCOCLOSE(avant,apres,a_min) || NCOEQ(avant,apres,de_param) || 
       NCOCLOSE(avant,apres,h_100) || NCOCLOSE(avant,apres,Omega_b) || NCOCLOSE(avant,apres,n_spec) ||
       NCOCLOSE(avant,apres,normalization) || 
       NCOCLOSE(avant,apres,Omega_nu_mass) || NCOCLOSE(avant,apres,Neff_nu_mass))
     return 1;

   if (change_w_de(avant, apres)) return 1;

   return 0;
}

double P_NL(cosmo *self, double a, double k, error **err)
{
   double p_cb;

   testErrorRet(self->nonlinear == pd96 && self->transfer == eisenhu_osc, ce_transfer,
         "Nonlinear mode 'pd96' not possible with transfer function 'eisenhu_osc'\n"
         "(Choose e.g. 'eisenhu' instead)",
         *err, __LINE__, -1.0);

   /* Matter power spectrum */
   switch (self->nonlinear) {

      case linear : case pd96 : case smith03 : case smith03_de : case smith03_revised :
         p_cb = P_NL_fitting(self, a, k, err);
         forwardError(*err, __LINE__, -1.0);
         break;

      case coyote10 : case coyote13 :
         p_cb = P_NL_coyote(self, a, k, err);
         forwardError(*err, __LINE__, -1.0);
         break;

      default :
         *err = addErrorVA(ce_unknown, "Unknown nonlinear flag %d", *err, __LINE__, self->nonlinear);
         return -1.0;

   }

   return p_cb;
}

/* ============================================================ *
 * Dark (+ baryonic) matter power spectrum, for			*
 * snonlinear = linear, pd96, smith03,smith03_de, and           *
 * smith03_revised.						*
 * For a general power spectrum, use cosmo.c:P_NL.       	*
 * ============================================================ */
#define itermax       25
#define logrmin       -6.0
#define logrmax       5.0
#define eps_a         1.0e-5
#define epsi          1.0e-3 // 1e-5
#define precision_min 1.0e-6
#define precision_max 1.01e-11
double P_NL_fitting(cosmo* self, double a, double k_NL, error **err)
{
   double table_k[N_k], table_P[N_k], y2[N_k];
   double *table_slope = 0;
   double logkmin = 0.0, logkmax = 0.0, dk = 0.0, da = 0.0;
   double Delta_NL, Delta_L, k_L, lnk_NL=0.0;
   double precision = precision_min;

   double omm=0.0, omv=0.0, amp=0.0;
   double logr1, logr2 ,diff, rmid, sig, d1, d2;
   double rknl=0.0, rneff=0.0, rncur=0.0;

   double aa, klog, val, logrmid, EE, r, upper=0.0;
   int i,j, iter, golinear=-1;
   interTable2D *lP_NL;

   testErrorRetVA(k_NL<=0, ce_negative, "Negative k=%g", *err, __LINE__, -1.0, k_NL);
   testErrorRetVA(a<self->a_min || a>1.0,
         ce_wrongValue, "a=%g out of range", *err, __LINE__, -1.0, a);


   if (self->P_NL == NULL) {
      testErrorRetVA(!(self->nonlinear==linear || self->nonlinear==pd96 || self->nonlinear==smith03
               || self->nonlinear==smith03_de || self->nonlinear==smith03_revised),
            ce_unknown, "Unknown nonlinear mode %d",
            *err, __LINE__, 0, self->nonlinear);

      testErrorRet(self->nonlinear == coyote10, ce_wrongValue, "nonlinear mode 'coyote10' not valid in P_NL",
            *err, __LINE__, 0.0);
      testErrorRet(self->nonlinear == coyote13, ce_wrongValue, "nonlinear mode 'coyote13' not valid in P_NL",
            *err, __LINE__, 0.0);

      /* upper = (dlnP/dlnk)_{k=kmax}, for splines & extrapolation */
      /* Note that the in the range considered here the linear power
       * spectrum is still a bit shallower than k^(n-4) since T(k) has
       * not yet reached its asymptotic limit. 
       */

      da = (1.0 - self->a_min)/(self->N_a-1.0);
      aa = self->a_min;
      logkmin = log(k_min);
      logkmax = log(k_max);
      dk = (logkmax - logkmin)/(N_k-1.0);

      if (self->nonlinear==linear) {
         upper = self->n_spec-4.0;
      } else if (self->nonlinear==pd96) {
         upper = -2.5;
      } else if (self->nonlinear==smith03 || self->nonlinear==smith03_de || self->nonlinear==smith03_revised) {
         self->slope = init_interTable(self->N_a, self->a_min, 1.0, da, 1.0e31, 1.0e31, err);
         forwardError(*err,__LINE__,0);
         table_slope = self->slope->table;
      }
      lP_NL = init_interTable2D(self->N_a, self->a_min, 1.0, da, N_k, logkmin, logkmax, 
            dk, self->n_spec, upper, err);
      forwardError(*err,__LINE__,0);

      for (i=0; i<self->N_a; i++,aa+=da) {

         if (self->nonlinear==smith03 || self->nonlinear==smith03_de || self->nonlinear==smith03_revised) {
            /* Smith et al. (2003) */
            EE = Esqr(self, aa, 0, err);           forwardError(*err,__LINE__,0);
            omm = Omega_m_a(self, aa, EE, err);    forwardError(*err,__LINE__,0);
            omv = Omega_de_a(self, aa, EE, err);   forwardError(*err,__LINE__,0);
            golinear = 0;

            /* Find non-linear scale with iterative bisection */
iter_start:
            logr1 = logrmin;
            logr2 = logrmax;
            iter = 0;
            do {

               logrmid = (logr2+logr1)/2.0;
               rmid    = pow(10,logrmid);
               wint2(self,rmid, &sig, NULL, NULL, aa, 1, err, precision);
               if (isError(*err))
                  printError(stderr, *err);
               forwardError(*err, __LINE__, 0);
               diff = sig - 1.0;
               if(diff>epsi)
                  logr1 = dlog(rmid);
               if(diff<-epsi)
                  logr2 = dlog(rmid);
            } while (fabs(diff)>=epsi && ++iter<itermax);

            /* Non-linear scale not found */
            if (iter>=itermax) {
               /* Non-linear scale very very small, all scales basically
                * linear: set flag golinear and continue
                */
               if (logrmid-logrmin<EPSILON) {
                  /* Scale at low end */
                  golinear = 1;
                  upper    = table_slope[i] = self->n_spec-4.0;
               } else if (logrmax-logrmid<EPSILON) {
                  /* scale at high end */
                  *err = addError(ce_nonlin, "All scales are nonlinear. Something went wrong in P_NL",
                        *err,__LINE__);
                  return 0.0;
               } else {
                  /* Somewhere in between */
                  if (precision>precision_max) {
                     precision /= 10.0;
                     goto iter_start;
                  } else {
                     *err = addErrorVA(ce_noknl, "Nonlinear scale not found for a=%g, precision=%g", *err, __LINE__, aa, precision);
                     return 0.0;
                  }
               }
            } else {
               /* Spectral index & curvature at non-linear scale */
               wint2(self,rmid, &sig, &d1, &d2, aa, 0, err, precision_min);
               forwardError(*err,__LINE__,0);
               rknl  = 1.0/rmid;
               rneff = -3.0-d1;
               rncur = -d2;
               upper = table_slope[i] = slope_NL(rneff, rncur, omm, omv);
               forwardError(*err,__LINE__,0);
            }
         }
         klog = logkmin;
         for (j=0; j<N_k; j++, klog+=dk) {
            k_L = exp(klog);
            if (self->nonlinear==linear) {

               Delta_NL = P_L(self, aa, k_L, err)*k_L*k_L*k_L/(2.0*pi_sqr);
               forwardError(*err, __LINE__, 0);
               lnk_NL   = klog;

            } else if (self->nonlinear==pd96) {

               Delta_L  = P_L(self,aa,k_L,err)*k_L*k_L*k_L/(2.0*pi_sqr);
               forwardError(*err,__LINE__,0);
               Delta_NL = f_NL(self,Delta_L, aa, k_L, err);
               forwardError(*err,__LINE__,0);
               lnk_NL   = klog + 1.0/3.0*log(1.0 + Delta_NL);  /* PD (5) */

            } else if (self->nonlinear==smith03 || self->nonlinear==smith03_de ||
                  self->nonlinear==smith03_revised) {

               switch (self->transfer) {
                  case bbks : case eisenhu : case eisenhu_osc :
                     amp = D_plus(self, aa, 1, err);        forwardError(*err,__LINE__,0);
                     break;
                  default :
                     *err = addError(ce_transfer, "wrong transfer type", *err, __LINE__);
                     return 0;
               }

               /* MKDEBUG NEW: Argument aa in sm2_transfer removed. This was ignored anyway. *
                * D+(a) = amp, so probably ok */
               Delta_L = amp*amp*sm2_transfer(self, k_L, err);
               forwardError(*err, __LINE__, 0);
               if (golinear==0) {
                  halofit(k_L, rneff, rncur, rknl, Delta_L, omm, omv, &Delta_NL, self->nonlinear, aa, self, err);
                  forwardError(*err, __LINE__, 0);
               } else if (golinear==1) {
                  Delta_NL = Delta_L;
               }
               lnk_NL = klog;

            }

            table_k[j] = lnk_NL;
            table_P[j] = log(2*pi_sqr*Delta_NL) - 3.0*lnk_NL; /* PD (3) */

         }
         sm2_spline(table_k-1, table_P-1, N_k, self->n_spec, upper, y2-1, err);
         forwardError(*err,__LINE__,0);
         klog = logkmin;
         for (j=0; j<N_k; j++, klog += dk) {
            sm2_splint(table_k-1, table_P-1, y2-1, N_k, klog, &val,err);
            forwardError(*err,__LINE__,0);
            lP_NL->table[i][j] = val;
         }
      }
      self->P_NL = lP_NL;
   }

   klog = log(k_NL);

#ifdef fastxi
   if (self->nonlinear==smith03 || self->nonlinear==smith03_de || self->nonlinear==smith03_revised) {
      upper = interpol_wr(self->slope,a,err);
      forwardError(*err, __LINE__, 0);
      self->P_NL->upper = upper;
   }

   r = (a-self->a_min)/self->P_NL->dx1;   /* da */
   i = (int)(r+0.5);                      /* Round to next integer */
   if (fabs(r-i)<eps_a) {
      val = sm2_interpol(self->P_NL->table[i], self->P_NL->n2, self->P_NL->a2, self->P_NL->b2,
            self->P_NL->dx2, klog, self->n_spec, self->P_NL->upper, err);
      forwardError(*err, __LINE__, 0);
      testErrorRet(!finite(val), ce_infnan, "Power spectrum P_NL not finite", *err, __LINE__, -1);
      val = exp(val);
      return val;
   }
   else
#endif
      val = interpol2D(self->P_NL, a, klog, err);
   forwardError(*err, __LINE__, 0);
   val = exp(val);
   return val;
}
#undef itermax 
#undef logrmin 
#undef logrmax 
#undef eps_a  
#undef epsi
#undef precision_min
#undef precision_max

/* ============================================================ *
 * Calls the Coyote emulator and returns the non-linear power   *
 * spectrum.							*
 * ============================================================ */

void set_H0_Coyote(cosmo *self, error **err)
{
   cosmo *tmp;
   double h_100_WMAP5;

   tmp         = copy_parameters_only(self, err);   forwardError(*err, __LINE__,);
   h_100_WMAP5 = getH0fromCMB(self->Omega_m, self->Omega_b, self->w0_de, 0);
   self->h_100 = h_100_WMAP5;
   updateFrom(tmp, self, err);                      forwardError(*err, __LINE__,);
   free_parameters(&tmp);
}

double P_NL_coyote(cosmo *self, double a, double k, error **err)
{
   double K, h2, h_100_WMAP7, val;
   int i_max;

   /* Check flatness and w(z) = const */
   K = self->Omega_m + self->Omega_nu_mass + self->Omega_de - 1.0;
   testErrorRetVA(fabs(K)>EPSILON, coyote_flat, "Universe is not flat (K=%g)", *err, __LINE__, -1.0, K);
   testErrorRetVA(fabs(self->w1_de)>EPSILON, ce_de, "Dark-energy equation of state not constant (w_1=%g)",
         *err, __LINE__, -1.0, self->w1_de);
   testErrorRet(self->de_param == poly_DE, ce_de, "Dark-energy mode not consistent (no constant w)",
         *err, __LINE__, -1.0);

   h2  = self->h_100 * self->h_100;

   if (self->nonlinear == coyote10) {

      /* Check whether Hubble constant is consistent with CMB constraints */
      h_100_WMAP7 = getH0fromCMB(self->Omega_m * h2, self->Omega_b * h2, self->w0_de, 1);
      //h_100_WMAP7 = getH0fromCMB(self->Omega_m, self->Omega_b, self->w0_de, 0);

      /* Note: for some (numerical?) reasons, precision should be set lower than 1.0e-3 (in getH0fromCMB) */
      testErrorRetVA(fabs(h_100_WMAP7 - self->h_100) > 0.005, coyote_h,
            "Input Hubble parameter %g not consistent with CMB constraints (h=%g)",
            *err, __LINE__, -1.0, self->h_100, h_100_WMAP7);

      i_max = nsim - 1;

   } else {

      i_max = fr_nsim - 1;

   }

   /* Tabulated Coyote k's are [1/Mpc] */
   if (self->nonlinear == coyote10) {

      if (k < ksim[0] / self->h_100) {

         /* Large scales -> linear power spectrum */
         val = P_L(self, a, k, err);   forwardError(*err, __LINE__, -1.0);

      } else if (k > ksim[i_max] / self->h_100) {

         /* Small scales, option 1 -> Set power to zero */
         val = 0;

         /* Small scales, option 2-> Glue Smith03 P(k) a la Eifler 2011 */
         /*
            static cosmo *tmp = NULL;
            static double ratio;
            if (tmp == NULL) {
            tmp = copy_parameters_only(self, err);   forwardError(*err, __LINE__, -1.0);
            tmp->nonlinear = smith03_de;
            updateFrom(self, tmp, err);              forwardError(*err, __LINE__, -1.0);

            ratio  = P_NL_coyote5(self->Omega_m * h2, self->Omega_b * h2, self->n_spec, self->sigma_8,
            self->w0_de, self->h_100, a, ksim[i_max] / self->h_100, err);
            forwardError(*err, __LINE__, -1.0);
            ratio /= P_NL_fitting(tmp, a, ksim[i_max] / self->h_100, err);
            forwardError(*err, __LINE__, -1.0);
            printf("ratio = %g\n", ratio);
            }
            val  = P_NL_fitting(tmp, a, k, err);                 forwardError(*err, __LINE__, -1.0);
            val *= ratio;
          */

      } else {

         /* Coyote v1 */
         val = P_NL_coyote5(self->Omega_m * h2, self->Omega_b * h2, self->n_spec, self->sigma_8,
               self->w0_de, self->h_100, a, k, err);
         forwardError(*err, __LINE__, -1.0);

      }

   } else if (self->nonlinear == coyote13) {

      if (k < fr_ksim[0] / self->h_100) {

         /* Large scales -> linear power spectrum */
         val = P_L(self, a, k, err);   forwardError(*err, __LINE__, -1.0);

      } else if (k > fr_ksim[i_max] / self->h_100) {

         /* Small scales, option 1 -> Set power to zero */
         val = 0;

      } else {

         /* Coyote v2 */
         val = P_NL_coyote6(self->Omega_m * h2, self->Omega_b * h2, self->n_spec, self->sigma_8,
               self->w0_de, self->h_100, a, k, &(self->ystar_allz), err);
         forwardError(*err, __LINE__, -1.0);

      }

   } else {

      val = -1;
      *err = addError(ce_unknown, "Invalid nonlinear mode (should be coyote10 or coyote13 here)", *err, __LINE__);

   }

   return val;
}

double int_for_w(double a, void *intpar, error **err)
{
   double asqr, d;
   cosmoANDint *ci;

   ci = (cosmoANDint*)intpar;
   asqr = dsqr(a);
   d    = asqr*asqr*Esqr(ci->self, a, ci->i, err);
   forwardError(*err,__LINE__,-1);

   return 1.0/sqrt(d);
}

int change_w(cosmo *avant, cosmo *apres) {
   if (NCOEQ(avant,apres,de_param) || NCOEQ(avant,apres,N_a)) 
      return 1;
   if (NCOCLOSE(avant,apres,Omega_m) || NCOCLOSE(avant,apres,Omega_de) || NCOCLOSE(avant,apres,w0_de) ||
         NCOCLOSE(avant,apres,w1_de) || NCOCLOSE(avant,apres,a_min) || NCOCLOSE(avant,apres,Omega_nu_mass))
      return 1;
   /* MKDEBUG TODO: check N_poly_de */

   if (change_w_de(avant, apres)) return 1;

   return 0;

   /* To be strict: With wOmega=1 in Esqr, there is a dependence on h_100 */
}

/* Comoving distance to an object at scale factor a */
/* BS01 2.41, with a(z_2) = 0, a(z_1) = z, [w] = Mpc/h */
double w(cosmo* self, double a, int wOmegar, error **err)
{
   double *table;
   double da;
   double aa;
   int    i;
   double res;
   interTable *lw;
   cosmoANDint intpar;

   if (self->w==NULL) {
      da = (1.0-self->a_min)/(self->N_a-1.0);
      aa = self->a_min;
      lw = init_interTable(self->N_a, self->a_min, 1.0, da, 0.0, 0.0, err);
      forwardError(*err,__LINE__,0);
      table = lw->table;
      intpar.self = self;
      intpar.i    = wOmegar;

      for (i=0; i<self->N_a-1; i++, aa+=da) {
         table[i] = R_HUBBLE*sm2_qromberg(int_for_w, (void*)(&intpar), aa, 1.0, 1.0e-6, err);
         forwardError(*err,__LINE__,0);
      }

      table[self->N_a-1] = 0.0;
      self->w = lw;
   }

   res = interpol_wr(self->w, a, err);
   forwardError(*err, __LINE__, 0.0);

   return res;
}

/* dw/da in Mpc/h */
double dwoverda(cosmo *self, double a, error **err)
{
   double res;
   cosmoANDint intpar;

   intpar.self = self;
   intpar.i    = 0;    /* wOmegar */
   res = R_HUBBLE*int_for_w(a, (void*)&intpar, err);
   forwardError(*err,__LINE__,0);
   return res;
}

/* dr/dz, r=comoving distance, see Hamana et al. (2004), eq. (10) */
double drdz(cosmo *self, double a, error **err)
{
   double res;
   int wOmegar = 0;   

   res  = 1.0/sqrt(Esqr(self, a, wOmegar, err));
   forwardError(*err,__LINE__,-1);
   res  *= R_HUBBLE;

   return res;
}

/* dVol/dz */
double dvdz(cosmo *self, double a, error **err)
{
   double res;
   int wOmegar;

   wOmegar = 0;
   res   = dsqr(w(self, a, wOmegar, err));     forwardError(*err, __LINE__, -1.0);
   res  *= drdz(self, a, err);                 forwardError(*err, __LINE__, -1.0);

   return res;
}

/* BS01 2.4, 2.30 */
double f_K(cosmo* self, double w, error **err)
{
   double K, K_h, f;

   testErrorRet(w<0, ce_negative, "Comoving distance is negative", *err, __LINE__, -1.0);

   /* Curvature without radiation (photons + massless neutrinos) density */
   K = self->Omega_m + self->Omega_nu_mass + self->Omega_de - 1.0;

   if (K>EPSILON) {             /* closed */
      K_h = sqrt(K)/R_HUBBLE;
      testErrorRet(K_h*w>pi, ce_zmax, "z larger than zmax in closed Universe",
            *err, __LINE__, -1.0);
      f = 1.0/K_h*sin(K_h*w);
   } else if (K<-EPSILON) {     /* open */
      K_h = sqrt(-K)/R_HUBBLE;
      f = 1.0/K_h*sinh(K_h*w);
   } else {                     /* flat */
      f = w;
   }

   return f;
}

/* Luminosity distance in Mpc/h */
double D_lum(cosmo *self, double a, error **err)
{
   double dlum, ww;
   ww   = w(self, a, 0, err);          forwardError(*err, __LINE__, 0);
   dlum = f_K(self, ww, err)/a;        forwardError(*err, __LINE__, 0);
   return dlum;
}

/* Returns 1 if w(a) is outside of the conservative de prior range *
 * and 0 if it is inside.					   */

int test_range_de_conservative(cosmo *model, error **err)
{
   double w;

   /* Prior:
      1. w(a) >= -1
      2. w(a) <= -1/3 for a>=a_acc=2/3 (z=1/5)
    */

   testErrorRet(model->de_param!=linder, ce_unknown,
         "Conservative de prior at the moment only defined for de_param=linder",
         *err, __LINE__, 1);

   w = w_de(model, 1.0, err);   forwardError(*err, __LINE__, 0.0);
   if (w<-1 || w>-1.0/3.0) return 1;

   w = w_de(model, model->a_min, err);   forwardError(*err, __LINE__, 0.0);
   if (w<-1) return 1;

   w = w_de(model, a_acc, err);   forwardError(*err, __LINE__, 0.0);
   if (w>-1.0/3.0) return 1;

   return 0;
}

/* ============================================================ *
 * Returns the index of the (i_bin j_bin)-correlation if these  *
 * are ordered lexically as					*
 * 0 0								*
 * 0 1								*
 * ...								*
 * 0 n-1							*
 * 1 1								*
 * ...								*
 * n-1 n-1							*
 * ============================================================ */
int idx_zz(int i_bin, int j_bin, int Nzbin, error **err)
{
   testErrorRetVA(i_bin >= Nzbin || j_bin >= Nzbin || i_bin < 0 || j_bin < 0, ce_range,
		  "Redshift bin(s) (%d,%d) out of range [0; %d]",
		  *err, __LINE__, -1, i_bin, j_bin, Nzbin);

   return i_bin*Nzbin - i_bin*(i_bin-1)/2 + j_bin - i_bin;
}

/* ============================================================ *
 * Returns the index of the (i_bin j_bin k_bin)-correlation     *
 * if these are ordered lexically.				*
 * ============================================================ */
int idx_zzz(int i_bin, int j_bin, int k_bin, int Nzbin)
{
   int i, j, k, index;

   for (i=0,index=0; i<Nzbin; i++) {
      for (j=i; j<Nzbin; j++) {
	 for (k=j; k<Nzbin; k++,index++) {
	    if (i==i_bin && j==j_bin && k==k_bin) return index;
	 }
      }
   }

   return -1;
}
