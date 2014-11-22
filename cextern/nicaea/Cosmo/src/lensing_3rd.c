 /* =========================================================== *
 * lensing_3rd.c						*
 * Bispectrum from Hyper-Extended Perturbation Theory (HEPT)    *
 * See Scoccimarro&Couchman 2001				*
 * Old code from 3.cf/num/ps.c: Martin Kilbinger 2006 (later	*
 * hept.c).							*
 * Integrated into cosmo_pmc: Martin Kilbinger, Liping Fu 2010	*
 * ============================================================ */


#include "lensing_3rd.h"

#define fastxi



cosmo_3rd *init_parameters_3rd(double OMEGAM, double OMEGADE, double W0_DE, double W1_DE, 
			       double *W_POLY_DE, int N_POLY_DE,
			       double H100, double OMEGAB, double OMEGANUMASS,
			       double NEFFNUMASS, double NORM, double NSPEC,
			       int Nzbin, const int *Nnz, const nofz_t *nofz, double *par_nz,
			       nonlinear_t NONLINEAR, transfer_t TRANSFER,
			       growth_t GROWTH, de_param_t DEPARAM,
			       norm_t NORMMODE , bispmode_t BISPMODE, error **err)
{
   cosmo_3rd *res;
   int m;

   testErrorRet(BISPMODE==PT && NONLINEAR!=linear, ce_nonlin, "bispmode=PT and nonlinear!=linear not consistent",
		*err, __LINE__, NULL);

   res = malloc_err(sizeof(cosmo_3rd), err); forwardError(*err, __LINE__, NULL);

   res->lens = init_parameters_lens(OMEGAM, OMEGADE, W0_DE, W1_DE, W_POLY_DE, N_POLY_DE,
				    H100, OMEGAB, OMEGANUMASS,
				    NEFFNUMASS, NORM, NSPEC, Nzbin, Nnz, nofz, par_nz,
				    NONLINEAR, TRANSFER, GROWTH, DEPARAM,
				    NORMMODE, tomo_all, reduced_none, 0.0, err);
   forwardError(*err, __LINE__, NULL);

   res->bispmode = BISPMODE;

   /* Reset pre-computed tables */
   res->k_NL          = NULL;
   res->scale_NL_amin = -1.0;
   for (m=0; m<=2; m++)
     res->B_kappa[m] = NULL;

   return res;
}



cosmo_3rd *copy_parameters_3rd_only(cosmo_3rd *source, error **err)
{
   cosmo_3rd *res;

   res = init_parameters_3rd(source->lens->cosmo->Omega_m,source->lens->cosmo->Omega_de,
			     source->lens->cosmo->w0_de, source->lens->cosmo->w1_de,
			     source->lens->cosmo->w_poly_de, source->lens->cosmo->N_poly_de,
			     source->lens->cosmo->h_100, source->lens->cosmo->Omega_b,
			     source->lens->cosmo->Omega_nu_mass, source->lens->cosmo->Neff_nu_mass,
			     source->lens->cosmo->normalization, source->lens->cosmo->n_spec,
			     source->lens->redshift->Nzbin, source->lens->redshift->Nnz,
			     source->lens->redshift->nofz, source->lens->redshift->par_nz,
			     source->lens->cosmo->nonlinear, source->lens->cosmo->transfer,
			     source->lens->cosmo->growth, source->lens->cosmo->de_param, 
			     source->lens->cosmo->normmode, source->bispmode, err);
   forwardError(*err, __LINE__, NULL);

   return res;
}



void updateFrom_3rd(cosmo_3rd *avant, cosmo_3rd *apres, error **err)
{
   int m, Nzbin, Nzcorr;

   updateFrom_lens(avant->lens, apres->lens, err);
   forwardError(*err, __LINE__,);

   /* TODO: change_.. functions */
   del_interTable(&(apres->k_NL));
   apres->scale_NL_amin = -1.0;

   Nzbin  = apres->lens->redshift->Nzbin;
   Nzcorr = Nzbin*(Nzbin+1)*(Nzbin+2)/6;
   for (m=0; m<=2; m++) del_interTable2D_arr(&(apres->B_kappa[m]), Nzcorr);
}


#define NZBIN 1
#define NNZ 5
cosmo_3rd *set_cosmological_parameters_to_default_lens_3rd(error **err)
{
   cosmo_3rd *res;

   int    Nnz[NZBIN]        = {NNZ};
   double par_nz[NZBIN*NNZ] = {0.0, 6.0, 0.612, 8.125, 0.62};
   nofz_t nofz[NZBIN]       = {ymmk};

   res = init_parameters_3rd(0.25, 0.75, -1.0, 0.0, NULL, 0, 0.70, 0.044, 0.0, 0.0, 0.80, 1.0,
			     NZBIN, Nnz, nofz, par_nz, smith03, eisenhu, growth_de, linder,
			     norm_s8, SCOCOU, err);
   forwardError(*err, __LINE__, NULL);

   return res;
}
#undef NZBIN
#undef NNZ

void read_cosmological_parameters_lens_3rd(cosmo_3rd **self, FILE *F, error **err)
{
   cosmo_3rd *tmp;
   struct { char cosmo_lens_file[128]; char sbispmode[128]; } tmp2;
   config_element c = {0, 0.0, ""};
   int j;
   FILE *FD;

   tmp = set_cosmological_parameters_to_default_lens_3rd(err);
   forwardError(*err, __LINE__,);

   CONFIG_READ_S(&tmp2, cosmo_lens_file, s, F, c, err);
   if (strcmp(tmp2.cosmo_lens_file, "-")!=0) {
      FD = fopen_err(tmp2.cosmo_lens_file, "r", err);
      forwardError(*err, __LINE__,);
   } else {
      FD = F;
   }
   read_cosmological_parameters_lens(&tmp->lens, FD, err);
   forwardError(*err, __LINE__,);
   if (strcmp(tmp2.cosmo_lens_file, "-")!=0) fclose(FD);

   /* Now read 3rd-order parameters */
   CONFIG_READ_S(&tmp2, sbispmode, s, F, c, err);
   STRING2ENUM(tmp->bispmode, tmp2.sbispmode, bispmode_t, sbispmode_t, j, Nbispmode_t, err);

   *self = copy_parameters_3rd_only(tmp, err);
   forwardError(*err, __LINE__,);

   free_parameters_3rd(&tmp);
}

void free_parameters_3rd(cosmo_3rd **self)
{
   int m, Nzbin, Nzcorr;
   cosmo_3rd *s;

   s = *self;
   del_interTable(&s->k_NL);
   Nzbin  = s->lens->redshift->Nzbin;
   Nzcorr = Nzbin*(Nzbin+1)*(Nzbin+2)/6;
   for (m=0; m<=2; m++) del_interTable2D_arr(&(s->B_kappa[m]), Nzcorr);

   free_parameters_lens(&s->lens);
   free(s);
   s = NULL;
}


void dump_param_3rd(cosmo_3rd* self, FILE *F, error **err)
{
   dump_param_lens(self->lens, F, 1, err);
   forwardError(*err, __LINE__,);
}


double dcub(double a)
{
   return a*a*a;
}

/* See also cosmo:n_L */

double n_eff(cosmo *self, double k, int check, error **err)
{
   double hh, error, n;

   hh = k/20.0;

   /* the slope of P_L only depends on k, not on a */
   n = sm2_dfridr(P_L, k, hh, &error, 0.9, self, err); forwardError(*err,__LINE__,0);
   n *= k/P_L(self, 0.9, k, err);                      forwardError(*err,__LINE__,0);

   if (check==1 && (n<-N_EFF_MIN || n>0)) {
      fprintf(stderr, "# n=%f at k=%f, ", n, k);
      fprintf(stderr, "n out of range in n_eff\n");
   }
   return n;
}


/* cf. Sco+Cou (9) */

double Q3(double n)
{
   return (4.0-pow(2.0,n))/(1.0+pow(2.0,n+1.0));
}



/* used for root finding in scale_NL */

double temp_NL(double a, double x, cosmo *self, error **err)
{
   double f;
   /* see DHS 03 */
   f = x*x*x/(2.0*pi_sqr)*P_L(self,a,x,err) - 1.0;
   forwardError(*err,__LINE__,0);
   return f;
}


/* ============================================================ *
 * Returns k_NL where 4 pi k_NL^3 P_L(k_NL)/(2pi)^3 = 1         *
 * ============================================================ */

double scale_NL(cosmo_3rd *self, double a, error **err)
{
   double k_NL, xacc, kmin, kmax, da, aa;
   int i;
   /* FILE *F; */

   if (self->k_NL==NULL) {

      /* F = fopen("scale_NL", "w"); */
      da = (1.0 - self->lens->cosmo->a_min)/(self->lens->cosmo->N_a-1.0);
      self->k_NL = init_interTable(self->lens->cosmo->N_a, self->lens->cosmo->a_min, 1.0, da, 0.0, 0.0, err);
      forwardError(*err,__LINE__,0);

      self->scale_NL_amin = 1.0;
      for (i=0,aa=self->lens->cosmo->a_min; i<self->lens->cosmo->N_a; i++,aa+=da) {

	 kmin = k_min;
	 kmax = k_max;

	 /* look for interval where sign changes */
	 while (temp_NL(aa, kmin, self->lens->cosmo, err)*temp_NL(aa, kmax, self->lens->cosmo, err)>=0) {
	    forwardError(*err,__LINE__,0);
	    kmin *= 10.0;
	    kmax *= 10.0;
	    if (kmax>1.e5) {
	       self->k_NL->table[i] = -1.0;
	       /* goto next_i; */
	       continue;
	    }
	 }
	 xacc = kmin/100.0;
	 self->k_NL->table[i] = hept_rtbis(temp_NL, kmin, kmax, xacc, aa, self->lens->cosmo, err);
	 if (aa<self->scale_NL_amin) self->scale_NL_amin = aa;
	 forwardError(*err,__LINE__,0);

	 /* next_i: */
	 /* fprintf(F, "%d %f %e\n", i, aa, self->k_NL->table[i]); */

      }
      /* fclose(F); */

   }

   if (a<self->scale_NL_amin) return -1;

   k_NL = interpol_wr(self->k_NL, a, err);
   forwardError(*err,__LINE__,0);

   return k_NL;
}


/* ============================================================ *
 * The following three functions are the fitting formulae for   *
 * the non-constant coefficients for the quasilinear bispectrum *
 * from Scoccimarro&Couchman 2001, equations (6)-(8).		*
 * k must be given in h times inverse Hubble radius		*
 * ============================================================ */


double ascocou(cosmo_3rd *self, double a, double k, double *n, error **err)
{
   double  q3, k_NL, q4, s8, erg;


   if (a<A_NL_MIN) return 1.0;     /* PT */

   if (*n<N_EFF_MIN) {
      *n = n_eff(self->lens->cosmo, k, 0, err);
      forwardError(*err,__LINE__,0);
   }

   q3   = sqrt(0.7*Q3(*n));
   k_NL = scale_NL(self, a, err); forwardError(*err,__LINE__,0);
   if (k_NL<0) q4 = 0.0;
   else q4   = pow(k/(4.0*k_NL), *n + 3.5);

   /* sigma_8 = sqrt(sigma_8_sqr(8./3000.)) */
   /* BUG in earlier version??! */
   //s8   = dsqr(D_plus(self->lens->cosmo, a, 1.0, err))*self->lens->cosmo->sigma_8;
   s8   = D_plus(self->lens->cosmo, a, 1.0, err)*self->lens->cosmo->sigma_8;
   forwardError(*err,__LINE__,0);

   erg  = (1. + pow(s8, -0.2)*q3*q4)/(1.0 + q4);
   testErrorRet(!finite(erg), ce_infnan, "inf or nan encountered", *err, __LINE__, 0.0);

   return erg;
}



double bscocou(cosmo_3rd *self, double a, double k, double *n, error **err)
{
   double k_NL, q, erg;


   if (a<A_NL_MIN) return 1.0;     /* PT */

   k_NL = scale_NL(self, a, err); forwardError(*err,__LINE__,0);

   if (k_NL<0) q = 0.0;
   else q    = k/k_NL;

   if (*n<N_EFF_MIN) {
      *n = n_eff(self->lens->cosmo, k, 0, err);
      forwardError(*err,__LINE__,0);
   }

   erg  = (1.0 + 0.4*(*n+3.0)*pow(q, *n+3.0))/(1.0 + pow(q, *n+3.5));
   testErrorRet(!finite(erg), ce_infnan, "inf or nan encountered", *err, __LINE__, 0.0);

   return erg;
}



double cscocou(cosmo_3rd *self, double a, double k, double *n, error **err)
{
   double k_NL, twoq, erg;


   if (a<A_NL_MIN) return 1.0;     /* PT */

   if (*n<N_EFF_MIN) {
      *n = n_eff(self->lens->cosmo, k, 0, err);
      forwardError(*err,__LINE__,0);
   }

   k_NL = scale_NL(self, a, err); forwardError(*err,__LINE__,0);
   if (k_NL<0) twoq = 0;
   else twoq = 2.*k/k_NL;

   erg  = (1.0 + 4.5/(1.5 + dsqr((*n+3.0)*(*n+3.0)))
	   *pow(twoq, *n+3.0))/(1. + pow(twoq, *n+3.5));
   testErrorRet(!finite(erg), ce_infnan, "inf or nan encountered", *err, __LINE__, 0.0);

   return erg;
}

double int_for_B_kappa_bar0(double a, void *intpar, error **err)
{
   int n_bin[3], i;
   double res, ggg, fK, ell1, ell2, dwda, ww, z, p1, p2, lens_eff;
   cosmo3ANDtomo *cANDr12;
   cosmoANDint ci;
   cosmo_3rd *self;

   cANDr12 = (cosmo3ANDtomo*)intpar;
   self = cANDr12->self;
   ell1 = cANDr12->r1;
   ell2 = cANDr12->r2;
   for (i=0; i<3; i++) n_bin[i] = cANDr12->n_bin[i];

   ci.self = self->lens->cosmo;
   ci.i    = 0;
   dwda    = R_HUBBLE*int_for_w(a, (void*)(&ci), err);        forwardError(*err, __LINE__, 0.0);
   for (i=0,ggg=1.0; i<3; i++) {
      ggg *= G(self->lens, a, n_bin[i], err);
      forwardError(*err,__LINE__,0.0);
      if (ggg==0.0) return 0.0;
   }

   ww = w(self->lens->cosmo, a, 0, err);       forwardError(*err, __LINE__, 0.0);
   fK = f_K(self->lens->cosmo, ww, err);       forwardError(*err, __LINE__, 0.0);
   z  = 1.0/a - 1.0;

   lens_eff = ggg;

   p1 = P_NL(self->lens->cosmo, a, ell1/fK, err);   forwardError(*err, __LINE__, 0.0);
   p2 = P_NL(self->lens->cosmo, a, ell2/fK, err);   forwardError(*err, __LINE__, 0.0);

   res = lens_eff * dwda/fK * p1 * p2;

   return res;
}



double int_for_B_kappa_bar1(double a, void *intpar, error **err)
{
   double erg, fK, ell1, ell2, ww, n;
   cosmo3ANDtomo *cANDr12;
   cosmo_3rd *self;

   cANDr12 = (cosmo3ANDtomo*)intpar;
   self = cANDr12->self;
   ell1 = cANDr12->r1;
   ell2 = cANDr12->r2;

   ww = w(self->lens->cosmo, a, 0, err);        forwardError(*err, __LINE__, 0.0);
   fK = f_K(self->lens->cosmo, ww, err);        forwardError(*err,__LINE__,0);

   n = N_EFF_MIN-1.0;
   erg  = ascocou(self, a, ell1/fK, &n, err);   forwardError(*err,__LINE__,0);
   n = N_EFF_MIN-1.0; /* Different k -> recalculate n */
   erg *= ascocou(self, a, ell2/fK, &n, err);   forwardError(*err,__LINE__,0); 
   testErrorRet(!finite(erg), ce_infnan, "inf or nan encountered", *err, __LINE__, 0.0);

   erg *= int_for_B_kappa_bar0(a, intpar, err); forwardError(*err,__LINE__,0);
   testErrorRet(!finite(erg), ce_infnan, "inf or nan encountered", *err, __LINE__, 0.0);

   return erg;
}



double int_for_B_kappa_bar2(double a, void *intpar, error **err)
{
   double erg, fK, ell1, ell2, ww, n;
   cosmo3ANDtomo *cANDr12;
   cosmo_3rd *self;


   cANDr12 = (cosmo3ANDtomo*)intpar;
   self = cANDr12->self;
   ell1 = cANDr12->r1;
   ell2 = cANDr12->r2;

   ww = w(self->lens->cosmo, a, 0, err);      forwardError(*err, __LINE__, 0.0);
   fK = f_K(self->lens->cosmo, ww, err);      forwardError(*err,__LINE__,0);

   n = N_EFF_MIN-1.0;
   erg  = bscocou(self, a, ell1/fK, &n, err); forwardError(*err,__LINE__,0);
   n = N_EFF_MIN-1.0;
   erg *= bscocou(self, a, ell2/fK, &n, err); forwardError(*err,__LINE__,0);
   testErrorRet(!finite(erg), ce_infnan, "inf or nan encountered", *err, __LINE__, 0.0);

   erg *= int_for_B_kappa_bar0(a, intpar, err); forwardError(*err,__LINE__,0);
   testErrorRet(!finite(erg), ce_infnan, "inf or nan encountered", *err, __LINE__, 0.0);

   return erg;
}



double int_for_B_kappa_bar3(double a, void *intpar, error **err)
{
   double erg, fK, ell1, ell2, ww, n;
   cosmo3ANDtomo *cANDr12;
   cosmo_3rd *self;


   cANDr12 = (cosmo3ANDtomo*)intpar;
   self = cANDr12->self;
   ell1 = cANDr12->r1;
   ell2 = cANDr12->r2;

   ww = w(self->lens->cosmo, a, 0, err);     forwardError(*err, __LINE__, 0.0);
   fK = f_K(self->lens->cosmo, ww, err);  forwardError(*err,__LINE__,0.0);

   n = N_EFF_MIN-1.0;
   erg  = cscocou(self, a, ell1/fK, &n, err); forwardError(*err,__LINE__,0.0);
   n = N_EFF_MIN-1.0;
   erg *= cscocou(self, a, ell2/fK, &n, err); forwardError(*err,__LINE__,0.0);
   testErrorRet(!finite(erg), ce_infnan, "inf or nan encountered", *err, __LINE__, 0.0);

   erg *= int_for_B_kappa_bar0(a, intpar, err); forwardError(*err,__LINE__,0);
   testErrorRet(!finite(erg), ce_infnan, "inf or nan encountered", *err, __LINE__, 0.0);

   return erg;
}


/* ============================================================ *
 * cf. Sco&Cou (5)						*
 * ============================================================ */


#define eps 1.0e-10
double F2eff(cosmo_3rd *self, double a, double k1, double k2, double cosphi, error **err)
{
   double f, x, n1, n2;

   switch (self->bispmode) {

      case SCOCOU :
	 n1 = n2 = N_EFF_MIN-1.0;

	 x  = ascocou(self,a,k1,&n1,err);       forwardError(*err, __LINE__, 0.0);
	 if (fabs(k1-k2)>eps) {
	    x *= ascocou(self,a,k2,&n2,err);    forwardError(*err, __LINE__, 0.0);
	 } else {
	    x *= x;
	 }
	 f  = 5.0/7.0*x;

	 x  = bscocou(self,a,k1,&n1,err);       forwardError(*err, __LINE__, 0.0);
	 if (fabs(k1-k2)>eps) {
	    x *= bscocou(self,a,k2,&n2,err);    forwardError(*err, __LINE__, 0.0);
	 } else {
	    x *= x;
	 }
	 f += 1.0/2.0*(k1/k2+k2/k1)*cosphi*x;

	 x  = cscocou(self,a,k1,&n1,err);       forwardError(*err, __LINE__, 0.0);
	 if (fabs(k1-k2)>eps) { 
	    x *= cscocou(self,a,k2,&n2,err);    forwardError(*err, __LINE__, 0.0);
	 } else {
	    x *= x;
	 }
	 f += 2.0/7.0*dsqr(cosphi)*x;
	 break;

      case PT :
	 f  = 5.0/7.0 + 1.0/2.0*(k1/k2+k2/k1)*cosphi + 2.0/7.0*dsqr(cosphi);
	 break;

      default :
	 *err = addErrorVA(ce_wrongmode, "Wrong mode %d", *err, __LINE__, self->bispmode);
	 return 0.0;

   }

   return f;
}

#undef eps


double F2(int q_kappa_term, double l1, double l2, error **err)
{
   if (q_kappa_term == 0) return 10.0/7.0;
   if (q_kappa_term == 1) return l1/l2 + l2/l1;
   if (q_kappa_term == 2) return 4.0/7.0;

   *err = addError(ce_wrongmode, "Wrong term", *err, __LINE__);
   return 0.0;

}



double F2cos(int q_kappa_term, double cosphi, error **err)
{
   if (q_kappa_term == 0) return 1.0;
   if (q_kappa_term == 1) return cosphi;
   if (q_kappa_term == 2) return dsqr(cosphi);

   *err = addErrorVA(ce_wrongmode, "Wrong term %d", *err, __LINE__, q_kappa_term);
   return 0.0;
}



/* ============================================================ *
 * cf. Sco&Cou (4)						*
 * ============================================================ */


double Q_123(cosmo_3rd *self, double a, double K1, double K2, double K3, error **err)
{
   double f2, Q, cosphi, P1, P2, P3;
   typedef double frr(cosmo *, double, double, error **);
   frr *PS[] = {&P_NL, &P_L};

   if (fabs(K1)<EPSILON2 || fabs(K2)<EPSILON2) {
      *err = addError(ce_negative,"K1 or K2 Negative !", *err, __LINE__);
      return 0;
   }

   P1 = PS[self->bispmode](self->lens->cosmo, a, K1, err); forwardError(*err,__LINE__,0);
   P2 = PS[self->bispmode](self->lens->cosmo, a, K2, err); forwardError(*err,__LINE__,0);
   P3 = PS[self->bispmode](self->lens->cosmo, a, K3, err); forwardError(*err,__LINE__,0);

   cosphi = (K1*K1 + K2*K2 - K3*K3)/(-2.0*K1*K2);
   testErrorRetVA(cosphi<-1-EPSILON || cosphi>1+EPSILON, ce_overflow, "|cosphi=%.10f|>1",
		  *err, __LINE__, 0.0, cosphi);

   f2  = F2eff(self, a, K1, K2, cosphi, err)*P1*P2;
   forwardError(*err,__LINE__,0);

   cosphi = (K2*K2 + K3*K3 - K1*K1)/(-2.0*K2*K3);
   testErrorRetVA(cosphi<-1-EPSILON || cosphi>1+EPSILON, ce_overflow, "|cosphi=%.10f|>1",
		  *err, __LINE__, 0.0, cosphi);

   f2  += F2eff(self, a, K2, K3, cosphi, err)*P2*P3;
   forwardError(*err,__LINE__,0);

   cosphi = (K3*K3 + K1*K1 - K2*K2)/(-2.*K3*K1);
   testErrorRetVA(cosphi<-1-EPSILON || cosphi>1+EPSILON, ce_overflow, "|cosphi=%.10f|>1",
		  *err, __LINE__, 0.0, cosphi);

   f2 += F2eff(self, a, K3, K1, cosphi, err)*P3*P1;
   forwardError(*err,__LINE__,0);

   Q = 2.0*f2/(P1*P2 + P2*P3 + P3*P1);

   return Q;
}



/* ============================================================ *
 * Returns b(k1, k2, cosphi) from eq. (21/PT) or (27/HEPT)	*
 * ============================================================ */

double bb(cosmo_3rd *self, double k1, double k2, double cosphi, int i_bin, int j_bin, int k_bin, error **err)
{
   int m;
   double b, btmp;

   if (self->bispmode==PT) {

      for (m=0,b=0.0; m<=2; m++) {
	 btmp  = F2(m, k1, k2, err);     forwardError(*err, __LINE__, 0.0);
	 btmp *= F2cos(m, cosphi, err);  forwardError(*err,__LINE__,0);
	 b += btmp;
      }
      b = b*B_kappa_bar(self, k1, k2, 0, i_bin, j_bin, k_bin, err);
      forwardError(*err,__LINE__,0);

   } else {

      for (m=0,b=0; m<=2; m++) {
	 btmp  = F2(m, k1, k2, err);                             forwardError(*err, __LINE__, 0.0);
	 btmp *= F2cos(m, cosphi, err);	                         forwardError(*err, __LINE__, 0.0);
	 btmp *= B_kappa_bar(self, k1, k2, m+1, i_bin, j_bin, k_bin, err);
	 forwardError(*err, __LINE__, 0.0);
	 b += btmp;
      }

   }

   return b;
}

/* ============================================================ *
 * Calculates the bispectrum of the convergence.                *
 * scocou=0 gives the bispectrum				*
 * with a=b=c=1,scocou=1,2,3 returns the integrated (bispectrum *
 * times a*a, b*b, c*c). bispmode=0, 1,   determines whether    *
 * SCO&COU or PT method is used.				*
 * The full convergence bispectrum is:				*
 * for HEPT:							*
 *  b(s1, s2, phi) = sum_{m=0}^2 F2(m, s1, s2)*F2cos(m, cosphi) *
 *                   *B_kappa_bar(s1, s2, m+1, SCOCOU)          *
 *								*
 * for perturbation theory:					*
 *  b(s1, s2, phi) = sum_{m=0}^2 F2(m, s1, s2)*F2cos(m, cosphi) *
 *                   *B_kappa_bar(s1, s2, 0, PT) + cycl		*
 *								*
 * <kappa(s1) kappa(s2) kappa(s3)> = (2pi)^2 [ b(s1, s2, phi12) *
 *                   + b(s2, s3, phi23) + b(s3, s1, phi31)]	*
 *								*
 * There occur small rounding errors when data is written       *
 * to a file. The returned value may thus depend on whether the *
 * bispectrum is calculated or read from file.			*
 * ============================================================ */


double B_kappa_bar(cosmo_3rd *self, double s1, double s2, int abc, int i_bin, int j_bin, int k_bin, error **err)
{
   double ds, logsmin, logsmax;
   double ss1, ss2, s1log, s2log, f1, f2, a, da;
   int    i, j, m, Nzcorr, Nzbin, ii, jj, kk, index;

   funcwithpars int_for_B[] = {&int_for_B_kappa_bar0, &int_for_B_kappa_bar1,
			       &int_for_B_kappa_bar2, &int_for_B_kappa_bar3};
   cosmo3ANDtomo intpar;

   testErrorRetVA(abc<0 || abc>3, ce_wrongmode, "Wrong term abc=%d", *err, __LINE__, 0.0, abc);
   testErrorRetVA(self->bispmode!=PT && self->bispmode!=SCOCOU, ce_wrongmode, "Wrong bispectrum mode %d",
		  *err, __LINE__, 0.0, self->bispmode);
   testErrorRetVA((abc==0 && self->bispmode==SCOCOU) || (abc>0 && self->bispmode==PT),
		  ce_wrongmode, "wrong combination of bispmode (%d) and term (abc=%d)", *err, __LINE__, 0.0,
		  self->bispmode, abc);

   Nzbin  = self->lens->redshift->Nzbin;
   Nzcorr = Nzbin*(Nzbin+1)*(Nzbin+2)/6;
   if (self->B_kappa[0]==NULL) {

      logsmin = log(s2_min);
      logsmax = log(s2_max);
      ds = (logsmax - logsmin)/((double)N_s2 - 1.0);

      /* if bispmode=PT, B_kappa[0] is used. Else, B_kappa[0,1,2] are used. */
      for (m=0; m<=2; m++) {
	 self->B_kappa[m] = init_interTable2D_arr(Nzcorr, N_s2, logsmin, logsmax, ds, N_s2, logsmin, logsmax, ds,
						  0.0, 0.0, err);
	 forwardError(*err,__LINE__,0.0);
      }

      da = (1.0 - self->lens->cosmo->a_min)/(self->lens->cosmo->N_a-1.0);

      intpar.self = self;

      for (ii=0; ii<Nzbin; ii++) {
	 intpar.n_bin[0] = ii;
	 for (jj=ii; jj<Nzbin; jj++) {
	    intpar.n_bin[1] = jj;
	    for (kk=jj; kk<Nzbin; kk++) {
	       intpar.n_bin[2] = kk;

	       index = idx_zzz(ii, jj, kk, Nzbin);

	       for (i=0,s1log=logsmin; i<N_s2; i++,s1log+=ds) {
		  ss1 = exp(s1log);
		  intpar.r1 = ss1;
		  for (j=0,s2log=logsmin; j<N_s2; j++,s2log+=ds) {
		     ss2 = exp(s2log);
		     intpar.r2 = ss2;

		     if (self->bispmode==PT) {

#ifndef fastxi
			if (self->lens->cosmo->a_min<0.5) {
			   f1 = sm2_qromberg(int_for_B[0], (void*)&intpar, self->lens->cosmo->a_min, 0.5, 1.0e-6, err);
			   forwardError(*err,__LINE__,0);
			   f2 = sm2_qrombergo(int_for_B[0], (void*)&intpar, 0.5, 1.0, sm2_midpntberg, 1.0e-7, err);
			   forwardError(*err,__LINE__,0);
			} else {
			   f1 = 0.0;
			   f2 = sm2_qrombergo(int_for_B[0], (void*)&intpar, self->lens->cosmo->a_min, 1.0,
					      sm2_midpntberg, 1.0e-7, err);
			   forwardError(*err,__LINE__,0);
			}
#else
			for (a=self->lens->cosmo->a_min,f1=0.0; a<1.0; a+=da) {
			   f1 += int_for_B[0](a, (void*)&intpar, err);
			   forwardError(*err, __LINE__, 0.0);
			}
			forwardError(*err,__LINE__,0);
			f1 *= da;
			f2  = 0.0;
#endif
			self->B_kappa[0][index]->table[i][j] = log(f1+f2);

		     } else { /* HEPT */

			for (m=0; m<=2; m++) {

#ifndef fastxi
			   if (self->lens->cosmo->a_min<0.5) {
			      f1 = sm2_qromberg(int_for_B[m+1], (void*)&intpar, self->lens->cosmo->a_min, 0.5, 1.0e-6, err);
			      forwardError(*err,__LINE__,0);
			      f2 = sm2_qrombergo(int_for_B[m+1], (void*)&intpar, 0.5, 1.0, sm2_midpntberg, 1.0e-7, err);
			      forwardError(*err,__LINE__,0);
			   } else {
			      f1 = 0.0;
			      f2 = sm2_qrombergo(int_for_B[m+1], (void*)&intpar, self->lens->cosmo->a_min, 1.0,
						 sm2_midpntberg, 1.0e-7, err);
			      forwardError(*err,__LINE__,0);
			   }
#else
			   for (a=self->lens->cosmo->a_min,f1=0.0; a<1.0; a+=da) {
			      f1 += int_for_B[m+1](a, (void*)&intpar, err);
			      forwardError(*err,__LINE__,0);
			   }
			   f1 *= da;
			   f2 = 0.0;
#endif
			   self->B_kappa[m][index]->table[i][j] = log(f1+f2);

			}
		     }
		  }
	       }
	    }
	 }
      }
   }

   s1log = log(s1);
   s2log = log(s2);
   index = idx_zzz(i_bin, j_bin, k_bin, Nzbin);
   f1 = interpol2D(self->B_kappa[self->bispmode==PT? 0: abc-1][index], s1log, s2log, err);
   forwardError(*err,__LINE__,0);

   //fprintf(stderr, "B_kappa: %g %g  %g\n", s1, s2, exp(f1));
   return exp(f1);
}


/* ============================================================ *
 * mathmatical functions (mainly from Num Rec)			*
 * ============================================================ */

#define JMAX 40

double hept_rtbis(double (*func)(double,double,cosmo*,error**),
		  double x1, double x2, double xacc, double a,
		  cosmo *self, error **err)
{
   int j;
   double dx,f,fmid,xmid,rtb;

   f=(*func)(a, x1, self, err);    forwardError(*err,__LINE__,0);
   fmid=(*func)(a, x2, self, err); forwardError(*err,__LINE__,0);
   if (f*fmid >= 0.0) {
      *err = addError(ce_rootbracket, "root must be bracketed for bisection", *err, __LINE);
      return 0.0;
   }
   rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
   for (j=1;j<=JMAX;j++) {
      fmid = (*func)(a, xmid=rtb+(dx *= 0.5), self, err); forwardError(*err,__LINE__,0);
      if (fmid <= 0.0) rtb=xmid;
      if (fabs(dx) < xacc || fmid == 0.0) return rtb;
   }
   *err = addError(ce_tooManySteps,"Too many bisections", *err, __LINE__);
   return 0.0;
}

#undef JMAX


/* ============================================================ *
 * 3d matter bispectrum. *Not* used for parallel sampling,      *
 * contains static variable!					*
 * ============================================================ */
#define eps 1.0e-10
#define N_K2 N_s2
#define K_MIN (s_min)
#define K_MAX (s_max)
double B_delta(cosmo_3rd *self, double k1, double k2, double cosbeta, double a,
	       error **err)
{
   double res, dk, k1log, k2log, kk1, kk2, n1, n2, x, ww, fK;
   static interTable2D *Bdelta[3] = {NULL, NULL, NULL};
   static double this_a=-1.0;
   int i, j, m;

   testErrorRetVA(fabs(cosbeta)>1.0+eps, math_overflow, "cosbeta=%.25f out of range",
		  *err, __LINE__, 0.0, cosbeta);
   assert(k1>0); assert(k2>0);

   if (cosbeta>1.0)  cosbeta = 1.0;
   if (cosbeta<-1.0) cosbeta = -1.0;


   if (fabs(this_a-a)>eps) {

      //fprintf(stderr, "Tabulating B_delta for a=%g\n", a);

      for (m=0; m<=2; m++) {
	 if (Bdelta[m]!=NULL) del_interTable2D(&(Bdelta[m]));
      }
      dk     = (log(K_MAX)-log(K_MIN))/((double)N_K2-1.0);
      for (m=0; m<=2; m++) {
	 Bdelta[m] = init_interTable2D(N_K2, log(K_MIN), log(K_MAX), dk, N_K2, log(K_MIN),
				       log(K_MAX), dk, 0.0, 0.0, err);
	 forwardError(*err, __LINE__, 0.0);
      }

      ww = w(self->lens->cosmo, a, 0, err);                      forwardError(*err, __LINE__, 0.0);
      fK = f_K(self->lens->cosmo, ww, err);                      forwardError(*err, __LINE__, 0.0);

      for (i=0,k1log=log(K_MIN); i<N_K2; i++,k1log+=dk) {
	 kk1 = exp(k1log);
	 for (j=0,k2log=log(K_MIN); j<N_K2; j++,k2log+=dk) {
	    kk2 = exp(k2log);


	    x = P_NL(self->lens->cosmo, a, kk1/fK, err);      forwardError(*err, __LINE__, 0.0);
	    if (i!=j) {
	       x *= P_NL(self->lens->cosmo, a, kk2/fK, err);  forwardError(*err, __LINE__, 0.0);
	    } else {
	       x *= x;
	    }

	    for (m=0; m<=2; m++) {
	       Bdelta[m]->table[i][j]  = x;
	    }

	    if (self->bispmode==SCOCOU) {
	       n1 = n2 = N_EFF_MIN-1.0;

	       x = ascocou(self, a, kk1/fK, &n1, err);     forwardError(*err, __LINE__, 0.0);
	       if (i!=j) {
		  x *= ascocou(self, a, kk2/fK, &n2, err); forwardError(*err, __LINE__, 0.0);
	       } else {
		  x *= x;
	       }

	       Bdelta[0]->table[i][j] *= x;

	       x = bscocou(self, a, kk1/fK, &n1, err);     forwardError(*err, __LINE__, 0.0);
	       if (i!=j) {
		  x *= bscocou(self, a, kk2/fK, &n2, err); forwardError(*err, __LINE__, 0.0);
	       } else {
		  x *= x;
	       }
	       Bdelta[1]->table[i][j] *= x;

	       x = cscocou(self, a, kk1/fK, &n1, err);     forwardError(*err, __LINE__, 0.0);
	       if (i!=j) {
		  x *= cscocou(self, a, kk2/fK, &n2, err); forwardError(*err, __LINE__, 0.0);
	       } else {
		  x *= x;
	       }
	       Bdelta[2]->table[i][j] *= x;
	    } else {
	       /* PT-term = 1 */
	    }

	 }
      }
      for (i=0; i<N_K2; i++) {
	 for (j=0; j<N_K2; j++) {
	    for (m=0; m<3; m++) {
	       Bdelta[m]->table[i][j] = log(Bdelta[m]->table[i][j]);
	    }
	 }
      }

      this_a = a;	 
   }

   if (k1>K_MAX || k2>K_MAX) {
      return 0.0;
   }
   if (k1<K_MIN || k2<K_MIN) {
      return 0.0;

      res = 2.0*F2eff(self, a, k1, k2, cosbeta, err);
      forwardError(*err, __LINE__, 0);

      ww = w(self->lens->cosmo, a, 0, err);                      forwardError(*err, __LINE__, 0.0);
      fK = f_K(self->lens->cosmo, ww, err);                      forwardError(*err, __LINE__, 0.0);

      res  = P_NL(self->lens->cosmo, a, k1/fK, err);  forwardError(*err, __LINE__, 0.0);
      if (fabs(k1-k2)>eps) {
	 res *= P_NL(self->lens->cosmo, a, k2/fK, err);  forwardError(*err, __LINE__, 0.0);
      } else {
	 res *= res;
      }
      return res;
   }

   for (m=0,res=0.0; m<=2; m++) {
      x    = interpol2D(Bdelta[m], log(k1), log(k2), err); forwardError(*err, __LINE__, 0.0);
      x    = exp(x);
      x   *= F2(m, k1, k2, err);        forwardError(*err, __LINE__, 0.0);
      x   *= F2cos(m, cosbeta, err);    forwardError(*err, __LINE__, 0.0);
      res += x;
   }


   return res;
}

#undef eps
#undef N_K2
#undef K_MIN
#undef K_MAX

/* ============================================================ *
 * Aperture-mass statistics					*
 * ============================================================ */

/* ============================================================ *
 * FT of the Map-filter.					*
 * ============================================================ */
double Uhat_one(double x, filter_t wfilter)
{
   const double etamin = 0.1;
   double uhat;

   switch (wfilter) {
      case fgauss :
	 uhat = x*x/2.0*exp(-x*x/2.0);
	 return uhat;

      case fpoly :
	 if (x<etamin) uhat = x*x/16.0*(1.0 - x*x/20.0);
	 else uhat = 24.0*gsl_sf_bessel_Jn(4, x)/(x*x);
	 return uhat;

      case ftophat :
	 if (x<etamin) uhat = 1.0-x*x/8.0;
	 else uhat = 2.0*gsl_sf_bessel_J1(x)/x;
	 return uhat;

      case fdelta :
	 return 1.0;
	 
      default :
	 assert(0);
   }
}

/* For Gaussian filter: Interpolation was not precise enough for chi^2 ! */
#define Neta 10000
#define etamax 4.0
double Uhat(double eta, filter_t wfilter)
{
   double x;
   x = Uhat_one(eta, wfilter);
   return x;
}
#undef Neta
#undef etamax

/* cyclically permutes the vector x[i] i=offset..offset+3 */
void permute3(double *x, int offset)
{
   int i;
   double tmp;

   for (tmp=x[2+offset],i=2+offset; i>offset; i--) {
      x[i] = x[i-1];
   }
   x[offset] = tmp;
}


double int_for_map3_3d(double x[], size_t dim, void *intpar)
{
   const double eps = 1.0e-15;

   double R[3], res, l1, l2, phi, cp, det;
   filter_t wfilter;
   int i, m, n_bin[3];
   cosmo_3rd *self;
   error **err;
   cosmo3ANDmorestuff *extra;

   //assert(intbeta==1);

   l1 = exp(x[0]); l2 = exp(x[1]); phi = x[2];

   extra = (cosmo3ANDmorestuff*)intpar;
   self  = extra->self;
   for (i=0; i<3; i++) R[i] = extra->R[i];
   for (i=0; i<3; i++) n_bin[i] = extra->n_bin[i];
   wfilter = extra->wfilter;
   m       = extra->m;
   err     = extra->err;

   res  = l1*l1*Uhat(l1*R[0], wfilter);

   /* MKDEBUG TODO: Get rid of inf tests */
   testErrorRet(!finite(res), ce_infnan, "inf or nan", *err, __LINE__, 0.0);

   res *= l2*l2*Uhat(l2*R[1], wfilter);

   testErrorRet(!finite(res), ce_infnan, "inf or nan", *err, __LINE__, 0.0);

   cp = cos(phi);
   det = dsqr(l1*R[2]) + dsqr(l2*R[2]) + 2.0*l1*l2*dsqr(R[2])*cp;
   if (det<eps) return 0.0;

   res *= Uhat(sqrt(det), wfilter);
   testErrorRet(!finite(res), ce_infnan, "inf or nan", *err, __LINE__, 0.0);

   for (i=0; i<m; i++) res *= cp;

   res *= F2(m, l1, l2, err);
   forwardError(*err, __LINE__, 0.0);
   testErrorRet(!finite(res), ce_infnan, "inf or nan", *err, __LINE__, 0.0);

   res *= B_kappa_bar(self, l1, l2, self->bispmode==PT?0:m+1, n_bin[0], n_bin[1], n_bin[2], err);
   forwardError(*err, __LINE__, 0.0);

   testErrorRet(!finite(res), ce_infnan, "inf or nan", *err, __LINE__, 0.0);

   return res;
}



double map3_perm(cosmo_3rd *self, double R[3], int i_bin, int j_bin, int k_bin, filter_t wfilter, error **err)
{
   const int calls = 5000;

   double xl[3], xu[3], m3, res, abserr, chisqr, dim;
   int i, m;
   cosmo3ANDmorestuff intpar;

   gsl_monte_vegas_state *mc_state;
   gsl_monte_function F;
   const gsl_rng_type *T;
   gsl_rng *r;

   /* initialize workspace */
   //if (intbeta==0) dim = 2;
   //else
   dim = 3;
   mc_state = gsl_monte_vegas_alloc(dim);
   gsl_monte_vegas_init(mc_state);

   /* integrand function and parameters */
   //if (intbeta==0) F.f = int_for_map3;
   //else
   F.f      = int_for_map3_3d;
   F.dim    = dim;
   F.params = (void*)&intpar;

   /* random number generator */
   gsl_rng_env_setup();
   T = gsl_rng_default;
   r = gsl_rng_alloc(T);

   /* additional parameters */
   intpar.self = self;
   for (i=0; i<3; i++) intpar.R[i] = R[i];
   intpar.n_bin[0] = i_bin;
   intpar.n_bin[1] = j_bin;
   intpar.n_bin[2] = k_bin;
   intpar.wfilter  = wfilter;
   intpar.err  = err;

   /* l-integration boundaries */
   xl[0] = xl[1] = log(s2_min);
   xu[0] = xu[1] = log(s2_max);

   /* phi-integration boundaries (unused for wfilter=gauss) */
   xl[2] = 0.0; xu[2] = twopi;

   /* Vegas MC integration */

   for (m=0,res=0.0; m<=2; m++) {
      intpar.m = m;
      gsl_monte_vegas_integrate(&F, xl, xu, dim, calls, r, mc_state, &m3, &abserr);
      forwardError(*err,__LINE__,0);
      chisqr = mc_state->chisq;

      if (0) {
	 if (abserr>fabs(m3)) {
	    fprintf(stderr, "Warning: |m3|<abserr (%e<%e)\n", fabs(m3), abserr);
	 } else if (chisqr>5.0) {
	    fprintf(stderr, "Warning: chi^2 (%f)>2, (m3,abserr)=(%e,%e)\n",
		    chisqr, m3, abserr);
	 }
      }

      res += m3;
   }

   testErrorRet(!finite(res), ce_infnan, "inf or nan", *err, __LINE__, -1.0);

   gsl_monte_vegas_free(mc_state);
   gsl_rng_free(r);

   res = res/(twopi*twopi*twopi);
   return res;
}



double map3(cosmo_3rd *self, double R[3], int i_bin, int j_bin, int k_bin, filter_t wfilter, error **err)
{
   const double eps = 1.0e-3*arcmin;
   double myR[3], m3;
   int i;

   for (i=0; i<3; i++) myR[i] = R[i];

   /* Diagonal, all angles equal */
   if (fabs(R[0]-R[1])<eps && fabs(R[1]-R[2])<eps) {
      m3 = 3.0*map3_perm(self, R, i_bin, j_bin, k_bin, wfilter, err);
      forwardError(*err,__LINE__,0);
      return m3;
   }

   for (i=0,m3=0.0; i<3; i++) {
      m3 += map3_perm(self, myR, i_bin, j_bin, k_bin, wfilter, err);
      forwardError(*err,__LINE__,0);
      permute3(myR, 0);
   }
   return m3;
}


/* ============================================================ *
 * Third-order lensing functions. So far no tomography.         *
 * ============================================================ */
double lensing_signal_3rd(cosmo_3rd *self, double theta[3], int i_bin, int j_bin, int k_bin, error **err)
{
   double res;

   res = map3(self, theta, i_bin, j_bin, k_bin, fgauss, err);
   forwardError(*err, __LINE__, 0.0);

   testErrorRet(!finite(res), ce_infnan, "inf or nan", *err, __LINE__, -1.0);

   return res;
}

double chi2_lensing_3rd(cosmo_3rd *self, datcov *dc, error **err)
{
   double *data_minus_model, theta[3], model, res;
   int j, k, Nzbin, in, i, i_bin, j_bin, k_bin;
   gsl_matrix_view A;
   gsl_vector_view x;


   Nzbin = self->lens->redshift->Nzbin;
   testErrorRetVA(Nzbin!=dc->Nzbin, redshift_Nzbin,
		  "Number of redshift bins for model (%d) inconsistent with data (%d)",
		  *err, __LINE__, 0, Nzbin, dc->Nzbin);

   data_minus_model = malloc_err(sizeof(double)*dc->n, err);
   forwardError(*err, __LINE__, 0);

   switch (dc->type) {

      case map3gauss_diag:
	 for (i_bin=0,in=0; i_bin<Nzbin; i_bin++) {
	    for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
	       for (k_bin=j_bin; k_bin<Nzbin; k_bin++) {
		  for (j=0; j<dc->Ntheta; j++,in++) {

		     theta[0] = theta[1] = theta[2] = dc->theta[j];
		     model = lensing_signal_3rd(self, theta, i_bin, j_bin, k_bin, err);
		     forwardError(*err, __LINE__, -1.0);
		     data_minus_model[in] = dc->data[in] - model;

		     //printf("%d %d %d  %d   %g %g\n", i_bin, j_bin, k_bin, j, dc->data[in], model);

		  }
	       }
	    }
	 }
	 break;

      case map2gauss_map3gauss_diag:
	 /* <Map^2> */
	 for (i_bin=0,in=0; i_bin<Nzbin; i_bin++) {
	    for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
	       for (j=0; j<dc->Ntheta; j++,in++) {
		  model = lensing_signal(self->lens, dc->theta[j], i_bin, j_bin, map2gauss, decomp_eb_none, err);
		  forwardError(*err, __LINE__, -1.0);
		  data_minus_model[in] = dc->data[in] - model;
	       }
	    }
	 }
	 /* <Map^3_diag> */
	 for (i_bin=0,in=0; i_bin<Nzbin; i_bin++) {
	    for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
	       for (k_bin=j_bin; k_bin<Nzbin; k_bin++) {
		  for (j=0; j<dc->Ntheta2; j++,in++) {
		     theta[0] = theta[1] = theta[2] = dc->theta2[j];
		     model = lensing_signal_3rd(self, theta, i_bin, j_bin, k_bin, err);
		     forwardError(*err, __LINE__, -1.0);
		     data_minus_model[in] = dc->data[in] - model;
		  }
	       }
	    }
	 }
	 break;

      case map2gauss_map3gauss:
	 testErrorRetVA(Nzbin!=1, redshift_Nzbin,
			"Number of redshift bins (%d) cannot be larger than 1 for general third-order aperture mass",
			*err, __LINE__, 0, Nzbin);

	 /* <Map^2> */
	 for (j=0,in=0; j<dc->Ntheta; j++,in++) {
	    model = lensing_signal(self->lens, dc->theta[j], 0, 0, map2gauss, decomp_eb_none, err);
	    forwardError(*err, __LINE__, -1.0);
	    data_minus_model[in] = dc->data[in] - model;
	 }
	 /* <Map^3> */
	 for (i=0; i<dc->Ntheta2; i++) {
	    theta[0] = dc->theta2[i];
	    for (j=i; j<dc->Ntheta2; j++) {
	       theta[1] = dc->theta2[j];
	       for (k=j; k<dc->Ntheta2; k++,in++) {
		  theta[2] = dc->theta2[k];
		  model = lensing_signal_3rd(self, theta, 0, 0, 0, err);
		  forwardError(*err, __LINE__, -1.0);
		  data_minus_model[in] = dc->data[in] - model;
	       }
	    }
	 }
	 break;

      case map3gauss:

	 testErrorRetVA(Nzbin!=1, redshift_Nzbin,
			"Number of redshift bins (%d) cannot be larger than 1 for general third-order aperture mass",
			*err, __LINE__, 0, Nzbin);

	 for (i=0,in=0; i<dc->Ntheta; i++) {
	    theta[0] = dc->theta[i];
	    for (j=i; j<dc->Ntheta; j++) {
	       theta[1] = dc->theta[j];
	       for (k=j; k<dc->Ntheta; k++,in++) {
		  theta[2] = dc->theta[k];
		  model = lensing_signal_3rd(self, theta, 0, 0, 0, err);
		  forwardError(*err, __LINE__, -1.0);
		  data_minus_model[in] = dc->data[in] - model;
	       }
	    }
	 }
	 break;

      default:
	 *err = addError(ce_unknown, "lensdata type 'map3gauss' not yet supported", *err, __LINE__);
	 return 0.0;
   }

   x = gsl_vector_view_array(data_minus_model, dc->n);
   A = gsl_matrix_view_array(dc->cov, dc->n, dc->n);

   /* Calculate L^{-1} * (data-model) */
   gsl_blas_dtrsv(CblasLower, CblasNoTrans, CblasNonUnit, &A.matrix, &x.vector);
   /* Square of the above */
   for (i=0,res=0.0; i<dc->n; i++) {
      res += dsqr(gsl_vector_get(&x.vector, i));
   }

   testErrorRetVA(res<0.0, math_negative, "Negative chi^2 %g. Maybe the covariance matrix is not positive",
		  *err, __LINE__, -1.0, res);

   testErrorRet(!finite(res), ce_infnan, "inf or nan", *err, __LINE__, -1.0);
   free(data_minus_model);

   /* det C ... */
   res = -0.5*res;

   return res;
}

