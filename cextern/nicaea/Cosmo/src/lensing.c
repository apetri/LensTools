/* ============================================================ *
 * lensing.c							*
 *								*
 * Martin Kilbinger, Karim Benabed 2006 - 2012			*
 * With many thanks to P. Schneider, J. Hartlap and P. Simon.	*
 * ============================================================ */

#include "lensing.h"


/* TODO: Not hard-code these halomodel parameters */
#define C0        9.0
#define ALPHA_NFW 1.0
#define BETA_NFW  0.13
#define MASSFCT   st2
#define HALO_BIAS halo_bias_sc


/* ============================================================ *
 * Creates and returns a new cosmo_lens structure with          *
 * parameters given in function call.                           *
 * ============================================================ */
cosmo_lens *init_parameters_lens(double OMEGAM, double OMEGADE, double W0_DE, double W1_DE,
				 double *W_POLY_DE, int N_POLY_DE,
				 double H100, double OMEGAB, double OMEGANUMASS,
				 double NEFFNUMASS, double NORM, double NSPEC,
				 int Nzbin, const int *Nnz, const nofz_t *nofz, double *par_nz,
				 nonlinear_t NONLINEAR, transfer_t TRANSFER,
				 growth_t GROWTH, de_param_t DEPARAM,
				 norm_t NORMMODE, tomo_t TOMO, reduced_t REDUCED, double Q_MAG_SIZE,
				 ia_t IA, ia_terms_t IA_TERMS, double A_IA, error **err)
{
   cosmo_lens *res;
   double amin;

   res = malloc_err(sizeof(cosmo_lens), err);  forwardError(*err, __LINE__, NULL);

   res->redshift = init_redshift(Nzbin, Nnz, nofz, par_nz, NULL, err); forwardError(*err, __LINE__, NULL);

   amin = get_amin(res->redshift, err);
   forwardError(*err, __LINE__, NULL);

   res->cosmo = init_parameters(OMEGAM, OMEGADE, W0_DE, W1_DE, W_POLY_DE, N_POLY_DE,
				H100, OMEGAB, OMEGANUMASS,
				NEFFNUMASS, NORM, NSPEC, NONLINEAR, TRANSFER, GROWTH, DEPARAM,
				NORMMODE, amin, err);
   forwardError(*err, __LINE__, NULL);

   /* Lensing parameters */
   res->tomo       = TOMO;
   res->reduced    = REDUCED;
   res->q_mag_size = Q_MAG_SIZE;

   res->ia         = IA;
   res->ia_terms   = IA_TERMS;
   res->A_ia       = A_IA;

   if (NONLINEAR==halodm) {
      res->hm = init_parameters_hm(res->cosmo->Omega_m, res->cosmo->Omega_de, res->cosmo->w0_de, res->cosmo->w1_de,
				   res->cosmo->w_poly_de, res->cosmo->N_poly_de,
				   res->cosmo->h_100, res->cosmo->Omega_b, res->cosmo->Omega_nu_mass,
				   res->cosmo->Neff_nu_mass, res->cosmo->normalization, res->cosmo->n_spec,
				   res->redshift->Nzbin, res->redshift->Nnz, res->redshift->nofz,
				   res->redshift->par_nz,
				   //[jean]---------------------------------------------------------------
				   //This is for the comoving volume in hod_FFTLog.c
				   -1, -1,
				   //---------------------------------------------------------------------
				   halodm, res->cosmo->transfer, res->cosmo->growth, res->cosmo->de_param,
				   res->cosmo->normmode,
				   C0, ALPHA_NFW, BETA_NFW, MASSFCT, HALO_BIAS, 0.0, 0.0, 0.0, 0.0, 0.0,
				   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,-1.0, 0.0, 0.0, 0.0,
				   hod_none, 60.0, err);
      
      forwardError(*err, __LINE__, NULL);
   } else {
      res->hm = NULL;
   }

   /* Reset pre-computed tables */
   res->g_source  = NULL;
   res->Pshear    = NULL;
   res->Pg1       = NULL;
   res->xiP       = NULL;
   res->xiM       = NULL;
   res->gamma     = NULL;
   res->map_gauss = NULL;
   res->map_poly  = NULL;
   res->c_cosebi  = NULL;
   res->psimin_cosebi = res->psimax_cosebi = 0.0; 
   res->N_cosebi  = 0;
 
   consistency_parameters_lens(res, err);    forwardError(*err, __LINE__, NULL);

   return res;
}

/* ============================================================ *
 * Checks for consistent parameters in structure cosmo_lens.    *
 * ============================================================ */
void consistency_parameters_lens(const cosmo_lens *self, error **err)
{
   testErrorRet(self->ia == ia_none && self->ia_terms != ia_undef, lensing_ia,
		"IA terms should be 'ia_undef' for no intrinsic alignment",
		*err, __LINE__,);
   testErrorRet(self->ia != ia_none && self->ia_terms == ia_undef, lensing_ia,
		"IA terms cannot be 'ia_undef' for intrinsic alignment",
		*err, __LINE__,);
}

/* ============================================================ *
 * Creates and returns a new cosmo_lens structure with          *
 * parameters copied from source.                               *
 * ============================================================ */
cosmo_lens* copy_parameters_lens_only(cosmo_lens* source, sm2_error **err)
{
   cosmo_lens *res;

   res = init_parameters_lens(source->cosmo->Omega_m,source->cosmo->Omega_de, source->cosmo->w0_de,
			      source->cosmo->w1_de, source->cosmo->w_poly_de, source->cosmo->N_poly_de,
			      source->cosmo->h_100, source->cosmo->Omega_b,
			      source->cosmo->Omega_nu_mass, source->cosmo->Neff_nu_mass,
			      source->cosmo->normalization, source->cosmo->n_spec,
			      source->redshift->Nzbin, source->redshift->Nnz, source->redshift->nofz,
			      source->redshift->par_nz,
			      source->cosmo->nonlinear, source->cosmo->transfer,
			      source->cosmo->growth, source->cosmo->de_param, 
			      source->cosmo->normmode, source->tomo, source->reduced, source->q_mag_size,
			      source->ia, source->ia_terms, source->A_ia, err);
   forwardError(*err, __LINE__, NULL);

   return res;
}

cosmo_lens* copy_parameters_lens(cosmo_lens* source, sm2_error **err)
{
   cosmo_lens *res;
   int Nzbin, Ncoeff, Nzcorr;

   res = copy_parameters_lens_only(source, err);             forwardError(*err, __LINE__, NULL);

   /* Reset cosmo and redshift */
   free_parameters(&res->cosmo);
   res->cosmo = copy_parameters(source->cosmo, err);         forwardError(*err, __LINE__, NULL);

   free_redshift(&res->redshift);
   res->redshift  = copy_redshift(source->redshift, err);    forwardError(*err, __LINE__, NULL);

   Nzbin = res->redshift->Nzbin;
   Nzcorr = Nzbin * (Nzbin+1) / 2;

   res->tomo       = source->tomo;
   res->reduced    = source->reduced;
   res->q_mag_size = source->q_mag_size;

   res->ia         = source->ia;
   res->ia_terms   = res->ia_terms;
   res->A_ia       = source->A_ia;

   if (source->cosmo->nonlinear==halodm) {
      res->hm = copy_parameters_hm(source->hm, err);          forwardError(*err, __LINE__, NULL);
   } else {
      res->hm = NULL;
   }

   /* TODO: Nzbin^2 -> Nzcorr should also work */
   res->Pshear = copy_interTable_arr(source->Pshear, Nzcorr, err);
   forwardError(*err, __LINE__, NULL);
   res->Pg1    = copy_interTable_arr(source->Pg1, Nzcorr, err);
   forwardError(*err, __LINE__, NULL);
   res->xiP    = copy_interTable_arr(source->xiP, Nzcorr, err);
   forwardError(*err, __LINE__, NULL);
   res->xiM    = copy_interTable_arr(source->xiM, Nzcorr, err);
   forwardError(*err, __LINE__, NULL);
   res->gamma  = copy_interTable_arr(source->gamma, Nzcorr, err);
   forwardError(*err, __LINE__, NULL);
   res->map_poly = copy_interTable_arr(source->map_poly, Nzcorr, err);
   forwardError(*err, __LINE__, NULL);
   res->map_gauss = copy_interTable_arr(source->map_gauss, Nzcorr, err);
   forwardError(*err, __LINE__, NULL);

   Ncoeff = NMAX_COSEBI * (NMAX_COSEBI + 5) / 2;
   memcpy(res->c_cosebi, source->c_cosebi, sizeof(double) * Ncoeff);
   res->psimin_cosebi = source->psimin_cosebi;
   res->psimax_cosebi = source->psimax_cosebi;
   res->N_cosebi      = source->N_cosebi;

   return res;
}

void read_cosmological_parameters_lens(cosmo_lens **self, FILE *F, error **err)
{
   cosmo_lens *tmp;
   struct { char cosmo_file[128], nofz_file[128], stomo[128], sreduced[128],
	sia[128], sia_terms[128]; } tmp2;
   config_element c = {0, 0.0, ""};
   int j;
   FILE *FD;

   tmp = set_cosmological_parameters_to_default_lens(err);
   forwardError(*err, __LINE__,);


   /* Cosmological parameters */
   CONFIG_READ_S(&tmp2, cosmo_file, s, F, c, err);
   if (strcmp(tmp2.cosmo_file, "-")!=0) {
      FD = fopen_err(tmp2.cosmo_file, "r", err);
      forwardError(*err, __LINE__,);
   } else {
      FD = F;
   }
   read_cosmological_parameters(&tmp->cosmo, FD, err);
   forwardError(*err, __LINE__,);
   if (strcmp(tmp2.cosmo_file, "-")!=0) fclose(FD);


   /* Redshift parameters */
   CONFIG_READ_S(&tmp2, nofz_file, s, F, c, err);
   if (strcmp(tmp2.nofz_file, "-")!=0) {
      FD = fopen_err(tmp2.nofz_file, "r", err);
      forwardError(*err, __LINE__,);
   } else {
      FD = F;
   }
   read_redshift_info(&(tmp->redshift), FD, err);
   forwardError(*err, __LINE__,);
   if (strcmp(tmp2.nofz_file, "-")!=0) fclose(FD);


   /* Lensing parameters */
   CONFIG_READ_S(&tmp2, stomo, s, F, c, err);
   STRING2ENUM(tmp->tomo, tmp2.stomo, tomo_t, stomo_t, j, Ntomo_t, err);

   CONFIG_READ_S(&tmp2, sreduced, s, F, c, err);
   STRING2ENUM(tmp->reduced, tmp2.sreduced, reduced_t, sreduced_t, j, Nreduced_t, err);
   if (tmp->reduced==reduced_K10) {
      CONFIG_READ(tmp, q_mag_size, d, F, c, err);
   }

   CONFIG_READ_S(&tmp2, sia, s, F, c, err);
   STRING2ENUM(tmp->ia, tmp2.sia, ia_t, sia_t, j, Nia_t, err);
   switch (tmp->ia) {

      case ia_HS04 :
       CONFIG_READ_S(&tmp2, sia_terms, s, F, c, err);
       STRING2ENUM(tmp->ia_terms, tmp2.sia_terms, ia_terms_t, sia_terms_t, j, Nia_terms_t, err);
       CONFIG_READ(tmp, A_ia, d, F, c, err);
       break;

      default :
       tmp->ia_terms = ia_undef;
       tmp->A_ia     = 0.0;
       break;

   }


   *self = copy_parameters_lens_only(tmp, err);
   forwardError(*err, __LINE__,);

   free_parameters_lens(&tmp);
}

/* ============================================================ *
 * Updates cosmo_lens structure apres from avant: Deletes pre   *
 * calculated tables if corresponding parameters have changed,  *
 * so tables will be re-calculated when needed.                 *
 * ============================================================ */
void updateFrom_lens(cosmo_lens *avant, cosmo_lens *apres, error **err)
{
   int Nzbin, Nzcorr;

   Nzbin = apres->redshift->Nzbin;
   Nzcorr = Nzbin*(Nzbin+1)/2;

   if (change_g_source(avant,apres)) {
      del_interTable_arr(&(apres->g_source), Nzbin);
   }

   if (change_Pshear(avant,apres)) {
      del_interTable_arr(&(apres->Pshear), Nzcorr);
      del_interTable_arr(&(apres->Pg1), Nzcorr);
      del_interTable_arr(&(apres->xiP), Nzcorr);
      del_interTable_arr(&(apres->xiM), Nzcorr);
      del_interTable_arr(&(apres->gamma), Nzcorr);
      del_interTable_arr(&(apres->map_poly), Nzcorr);
      del_interTable_arr(&(apres->map_gauss), Nzcorr);
      free(apres->c_cosebi); apres->c_cosebi = NULL;
   }

   updateFrom(avant->cosmo, apres->cosmo, err);
   forwardError(*err, __LINE__,);

   updateFrom_redshift(avant->redshift, apres->redshift);

   if (avant->cosmo->nonlinear==halodm) {
      updateFrom_hm(avant->hm, apres->hm, err);
      forwardError(*err, __LINE__,);
   }
}


/* Copies parameters from model->cosmo to model->hm->cosmo if   *
 * cosmo->nonlinear==halomodel					*/
void copy_parameters_lenshm_cosmo(cosmo_lens *model, error **err)
{
   cosmo *new;

   if (model->cosmo->nonlinear!=halodm) return;

   new = copy_parameters(model->cosmo, err);       forwardError(*err, __LINE__,);
   free_parameters(&model->hm->cosmo);
   model->hm->cosmo = new;
}

#define NZBIN 1
#define NNZ 5
cosmo_lens *set_cosmological_parameters_to_default_lens(error **err)
{
   int    Nnz[NZBIN]        = {NNZ};
   double par_nz[NZBIN*NNZ] = {0.0, 6.0, 0.612, 8.125, 0.62};
   nofz_t nofz[NZBIN]       = {ymmk};
   cosmo_lens *res;

   res = init_parameters_lens(0.25, 0.75, -1.0, 0.0, NULL, 0, 0.70, 0.044, 0.0, 0.0, 0.80, 1.0,
			      NZBIN, Nnz, nofz, par_nz, smith03, eisenhu, growth_de, linder,
			      norm_s8, tomo_all, reduced_none, 0.0, ia_none, ia_undef, 0.0, err);
   forwardError(*err, __LINE__, NULL);

   return res;
}
#undef NNZ
#undef NZBIN

void free_parameters_lens(cosmo_lens** self)
{
   cosmo_lens *s;
   int Nzbin, Nzcorr;

   s = *self;

   Nzbin = s->redshift->Nzbin;
   Nzcorr = Nzbin*(Nzbin+1)/2;
   del_interTable_arr(&s->g_source, Nzbin);

   del_interTable_arr(&s->Pshear, Nzcorr);
   del_interTable_arr(&s->Pg1, Nzcorr);
   del_interTable_arr(&s->xiP, Nzcorr);
   del_interTable_arr(&s->xiM, Nzcorr);
   del_interTable_arr(&s->gamma, Nzcorr);
   del_interTable_arr(&s->map_poly, Nzcorr);
   del_interTable_arr(&s->map_gauss, Nzcorr);
   free(s->c_cosebi); s->c_cosebi = NULL;

   if (s->hm!=NULL) {
      free_parameters_hm(&s->hm);
   }

   free_redshift(&s->redshift);
   free_parameters(&s->cosmo);

   free(s);
   s = NULL;
}

void dump_param_lens(cosmo_lens* self, FILE *F, int wnofz, error **err)
{
   dump_param(self->cosmo, F);
   if (wnofz) dump_redshift(self->redshift, F, err);
   forwardError(*err, __LINE__,);
   fprintf(F, "# (s)tomo = (%s)%d (s)reduced=(%s)%d q_mag_size=%g (s)ia=(%s)%d (s)ia_terms=(%s)%d A_ia=%g\n",
	   stomo_t(self->tomo), self->tomo, sreduced_t(self->reduced), self->reduced, self->q_mag_size,
	   sia_t(self->ia), self->ia, sia_terms_t(self->ia_terms), self->ia_terms,
	   self->A_ia);
   if (self->cosmo->nonlinear==halodm) dump_param_only_hm(self->hm, F);
}

double int_for_g(double aprime, void *intpar, error **err)
{
   cosmo_lensANDintANDdouble* cANDdANDe;
   cosmo_lens* self;
   double ww, wprime, res, fKwp, a;
   int n_bin;
   
   cANDdANDe = (cosmo_lensANDintANDdouble*)intpar;
   self      = cANDdANDe->self;
   a         = cANDdANDe->r;
   n_bin     = cANDdANDe->i;

   if (aprime>=a) return 0.0;

   ww     = w(self->cosmo, a, 0, err);                        forwardError(*err, __LINE__, 0);
   wprime = w(self->cosmo, aprime, 0, err);                   forwardError(*err, __LINE__, 0);
   res    = prob(self->redshift, 1.0/aprime-1.0, n_bin, err); forwardError(*err, __LINE__, 0);
   res   *= f_K(self->cosmo, wprime-ww, err)/dsqr(aprime);    forwardError(*err, __LINE__, 0);
   fKwp   = f_K(self->cosmo, wprime, err);                    forwardError(*err, __LINE__, 0);

   return res/fKwp;
}

int change_g_source(cosmo_lens* avant, cosmo_lens* apres)
{
   if (change_w(avant->cosmo, apres->cosmo) || NCOEQ(avant->cosmo, apres->cosmo, N_a))
     return 1;
   if (change_prob(avant->redshift, apres->redshift))
     return 1;
   return 0;
}

/* ============================================================ *
 * See S98 2.9.							*
 * Returns integral (dimensionless)				*
 * int dw' n(w') f(w'-w)/f(w') = int da' n(a') f(w'-w)/f(w')    *
 *  = int da' n(z') / a'^2 f(w'-w)/f(w')			*
 * ============================================================ */
double g_source(cosmo_lens* self, double a, int n_bin, error **err)
{
   double *table, res;
   double da = 0.0;
   double aa, wdelta, fKz0, wa;
   int    i, nn;
   cosmo_lensANDintANDdouble intpar;
	 
   if (self->redshift->nofz[n_bin]==single) {
      aa = 1.0/(1.0+get_zmin(self->redshift, n_bin));       /* zmin = zmax = z0 */
      if (aa>=a) return 0.0;
      wdelta = w(self->cosmo, aa, 0, err);           forwardError(*err,__LINE__,0);
      wa     = w(self->cosmo, a, 0, err);            forwardError(*err,__LINE__,0);
      res    = f_K(self->cosmo, wdelta-wa, err);     forwardError(*err,__LINE__,0);
      fKz0   = f_K(self->cosmo, wdelta, err);        forwardError(*err,__LINE__,0);

      return res/fKz0;
   }

   if (self->g_source==NULL) {
      da = (1.0-self->cosmo->a_min)/(self->cosmo->N_a-1.0);
      self->g_source = init_interTable_arr(self->redshift->Nzbin, self->cosmo->N_a,
      				  self->cosmo->a_min, 1.0, da, 0.0, 0.0, err);
      forwardError(*err,__LINE__,0);
      for (nn=0; nn<self->redshift->Nzbin; nn++) {
          table       = self->g_source[nn]->table;
          table[0]    = 0.0;
          aa          = self->cosmo->a_min+da;
          intpar.self = self;
          intpar.i    = nn;

          for (i=1;i<self->cosmo->N_a-1;i++,aa+=da) {
             intpar.r  = aa;
             /* Precision decreased from 1e-6 to 1e-5 */
             table[i]  = sm2_qromberg(int_for_g, (void*)&intpar, self->cosmo->a_min, aa, 1.0e-5, err);
             forwardError(*err, __LINE__, 0.0);
          }

          table[self->cosmo->N_a-1] = 1.0;
      }
   }

   res = interpol_wr(self->g_source[n_bin], a, err);
   forwardError(*err,__LINE__,0);

   return res;
}

/* lens efficiency */
double G(cosmo_lens* self, double a, int n_bin, error **err)
{
   double res;
   res  = 1.5/dsqr(R_HUBBLE)*(self->cosmo->Omega_m+self->cosmo->Omega_nu_mass)/a;
   res *= g_source(self, a, n_bin, err);
   forwardError(*err, __LINE__, 0);
   return res;	
}

/* ============================================================ *
 * dP_kappa/da. Integrand for P_kappa.				*
 * ============================================================ */
double int_for_p_2(double a, void *intpar, error **err)
{
   double hoverh0, asqr, s, fKw, f, res, wa, gg;
   int i_bin, j_bin;
   cosmo_lensANDiid* cANDiid;
   cosmo_lens* self;


   cANDiid = (cosmo_lensANDiid*)intpar;
   self    = cANDiid->self;
   s       = cANDiid->r;
   i_bin   = cANDiid->i;
   j_bin   = cANDiid->j;
  
   testErrorRet(a>=1.0, ce_overflow, "Scale factor a>=1", *err, __LINE__, -1);
 
   asqr    = dsqr(a);
   wa      = w(self->cosmo, a, 0, err);          forwardError(*err, __LINE__, -1);
   fKw     = f_K(self->cosmo, wa, err);          forwardError(*err, __LINE__, -1);
   f       = s/fKw;
  
   hoverh0 = Esqr(self->cosmo, a, 0, err);       forwardError(*err, __LINE__,-1);
   hoverh0 = sqrt(hoverh0);

   gg   = g_source(self, a, i_bin, err);         forwardError(*err, __LINE__, -1);
   gg  *= g_source(self, a, j_bin, err);         forwardError(*err, __LINE__, -1);
   if (fabs(gg) < EPSILON1) return 0.0;

   res  = gg/(asqr*asqr)/hoverh0*R_HUBBLE;
   res *= P_NL_tot(self, a, f, err);             forwardError(*err,__LINE__, -1.0);

   testErrorRetVA(!finite(res), ce_overflow, "Integrand not finite at a=%g", *err, __LINE__, -1.0, a);

   return res;
}

/* ============================================================ *
 * Non-linear intrinsic alignment model                         *
 *                                                              *
 * C1 factor included in Pkappa with IGfact                     *
 * Includes extra powers of a outlined in Joachimi et al 2010   *
 * and Hirata and Seljak erratum 2010                           *
 * ============================================================ */


/* ============================================================ *
 * Integrand for GI power spectrum, HS04 model. BK08 eqs 4, 12	*
 * ============================================================ */
double int_for_p_GI(double a, void *intpar, error **err)
{
   double hoverh0, asqr, s, fKw, f, res, wa, n, g, gn, Om_o_D;
   int i_bin, j_bin;
   cosmo_lensANDiid* cANDiid;
   cosmo_lens* self;

   cANDiid = (cosmo_lensANDiid*)intpar;
   self    = cANDiid->self;
   s       = cANDiid->r;
   i_bin   = cANDiid->i;
   j_bin   = cANDiid->j;
  
   testErrorRet(a>=1.0, ce_overflow, "Scale factor a>=1", *err, __LINE__, -1);
 
   asqr    = dsqr(a); /* a*a */
   wa      = w(self->cosmo, a, 0, err);          forwardError(*err, __LINE__, -1);
   fKw     = f_K(self->cosmo, wa, err);          forwardError(*err, __LINE__, -1);
   f       = s/fKw;
  
   hoverh0 = Esqr(self->cosmo, a, 0, err);       forwardError(*err, __LINE__,-1);
   hoverh0 = sqrt(hoverh0);

   Om_o_D = self->cosmo->Omega_m/D_plus(self->cosmo, a, 1, err);     forwardError(*err,__LINE__,0)

   g   = g_source(self, a, i_bin, err);         forwardError(*err, __LINE__, -1);
   n   = prob(self->redshift, 1.0/a-1.0, j_bin, err);   forwardError(*err, __LINE__, -1);

   gn = g*n;

   g   = g_source(self, a, j_bin, err);         forwardError(*err, __LINE__, -1);
   n   = prob(self->redshift, 1.0/a-1.0, i_bin, err);   forwardError(*err, __LINE__, -1);

   res = (gn + g*n)/fKw;

   res *= Om_o_D;

   //HS04
   //res *= 1.0/a/asqr/hoverh0*R_HUBBLE; 
   //HS10 ->  * a^2
   //res *= 1.0/a/hoverh0*R_HUBBLE; 
   //BJ10
   res *= 1.0/a/asqr;

   res *= P_NL_tot(self, a, f, err);             forwardError(*err,__LINE__, -1);

   testErrorRetVA(!finite(res), ce_overflow, "Integrand not finite at a=%g", *err, __LINE__, -1.0, a);
	
   return res;
}

/* ============================================================ *
 * Integrand for II power spectrum, HS04 model. BK08 eqs. 5, 7, *
 * BJ10 eq B7.							*
 * ============================================================ */
double int_for_p_II(double a, void *intpar, error **err)
{
   double hoverh0, asqr, s, fKw, f, res, wa, nn, Om_o_D;
   int i_bin, j_bin;
   cosmo_lensANDiid* cANDiid;
   cosmo_lens* self;

   cANDiid = (cosmo_lensANDiid*)intpar;
   self    = cANDiid->self;
   s       = cANDiid->r;
   i_bin   = cANDiid->i;
   j_bin   = cANDiid->j;
  
   testErrorRet(a>=1.0, ce_overflow, "Scale factor a>=1", *err, __LINE__, -1);
 
   asqr    = dsqr(a); /* a*a */
   wa      = w(self->cosmo, a, 0, err);          forwardError(*err, __LINE__, -1);
   fKw     = f_K(self->cosmo, wa, err);          forwardError(*err, __LINE__, -1);
   f       = s/fKw;
  
   hoverh0 = Esqr(self->cosmo, a, 0, err);       forwardError(*err, __LINE__,-1);
   Om_o_D  = self->cosmo->Omega_m/D_plus(self->cosmo, a, 1, err);     forwardError(*err,__LINE__,0);

   hoverh0 = sqrt(hoverh0);

   nn   = prob(self->redshift, 1.0/a-1.0, i_bin, err);   forwardError(*err, __LINE__, -1);
   nn  *= prob(self->redshift, 1.0/a-1.0, j_bin, err);   forwardError(*err, __LINE__, -1);

   res  = Om_o_D/fKw ;
   res  = dsqr(res);
     
   // HS04:
   //   res *= nn/asqr/hoverh0*R_HUBBLE;
   // HS10 ->  * a^4
   //res *= nn*asqr/hoverh0*R_HUBBLE;
   // BJ10
   
   res *= nn/asqr*hoverh0/R_HUBBLE;
   res *= P_NL_tot(self, a, f, err);             forwardError(*err,__LINE__, -1);

   testErrorRet(!finite(res), ce_overflow, "Value not finite", *err, __LINE__, -1);

   return res;
}


/* ============================================================ *
 * Intrinsic alignment model for <Map^3>, Semboloni et al.      *
 * 2008, 2010.                                                  *
 * ============================================================ */


/* ============================================================ *
 * Total (dark+baryon+neutrino) power spectrum.			*
 * ============================================================ */
double P_NL_tot(cosmo_lens *self, double a, double k, error **err)
{
   double p_cb, p_tot;

   /*
   // **** MKDEBUG ****
   double K_MAX = (2.0*pi/147.0);
   fprintf(stderr, "K_MAX = %g h/Mpc\n", K_MAX);
   if (k>K_MAX) return 0; // MK >
   */


   /* Matter power spectrum */
   switch (self->cosmo->nonlinear) {

      case linear : case pd96 : case smith03 : case smith03_de : case coyote10 :
      case coyote13 : case smith03_revised :
	      p_cb = P_NL(self->cosmo, a, k, err);
         forwardError(*err, __LINE__, -1.0);
         break;

      case halodm :
         p_cb = P1h_dm(self->hm, a, k, err);
         forwardError(*err, __LINE__, -1.0);
         p_cb += P2h_dm(self->hm, a, k, err);
         forwardError(*err, __LINE__, -1.0);
         break;

      default :
         *err = addErrorVA(ce_unknown, "Unknown nonlinear flag %d", *err, __LINE__, self->cosmo->nonlinear);
	      return -1.0;

   }

   switch (self->cosmo->transfer) {
      case bbks: case eisenhu: case eisenhu_osc:
	 p_tot = p_cb;
	 break;

      default :
	 *err = addError(ce_transfer, "wrong transfer type", *err, __LINE__);
	 return 0;
   }

   return p_tot;
}

int change_Pshear(cosmo_lens* avant, cosmo_lens* apres)
{
   if (change_w(avant->cosmo, apres->cosmo) || NCOEQ(avant->cosmo, apres->cosmo, N_a)) return 1;
   if (change_Tsqr(avant->cosmo, apres->cosmo)) return 1;
   if (change_prob(avant->redshift, apres->redshift)) return 1;
   if (change_P_NL(avant->cosmo, apres->cosmo)) return 1;
   if (NCOEQ(avant, apres, tomo)) return 1;
   if (NCOEQ(avant, apres, reduced)) return 1;
   if (NCOEQ(avant, apres, ia)) return 1;
   if (NCOEQ(avant, apres, ia_terms)) return 1;
   if (NCOCLOSE(avant, apres, A_ia)) return 1;

   return 0;
}

/* ============================================================ *
 * Returns the shear (cross-)power spectrum for Fourier         *
 * scale s redshift bins (i_bin,j_bin). Note that i_bin<=j_bin. *
 * See S98, eq. (3.4). Extrapolates outside of [s_min; s_max].	*
 * If self->reduced==K10, adds the first-order reduced-shear    *
 * correction according to K10. This correction is zero outside *
 * of [ELL_MIN_REDUCED=0.1; ELL_MAX_REDUCED=2e5].		*
 * ============================================================ */
double Pshear(cosmo_lens* self, double s, int i_bin, int j_bin, error **err)
{
   double *table;
   double ds, logsmin, logsmax, prefactor;
   double ss, slog, f1, f2;
   int    i, Nzbin, Nzcorr, ii, jj;
   cosmo_lensANDiid intpar;

   if (s<s_min || s>s_max) return 0.0;


   testErrorRetVA(s<=0, ce_negative, "Negative or zero 2d Fourier vector l=%g",* err, __LINE__, -1, s);
   testErrorRet(i_bin>j_bin, lensing_tomoij, "Pshear_ij defined for i<=j", *err, __LINE__, -1);

   Nzbin  = self->redshift->Nzbin;
   Nzcorr = Nzbin*(Nzbin+1)/2;
   if (self->Pshear==NULL) {

      logsmin   = log(s_min);
      logsmax   = log(s_max);
      ds        = (logsmax - logsmin)/(N_s - 1.0);
      prefactor = 9.0/4.0*dsqr((self->cosmo->Omega_m+self->cosmo->Omega_nu_mass)/R_HUBBLE/R_HUBBLE);

      self->Pshear = init_interTable_arr(Nzcorr, N_s, logsmin, logsmax, ds, 1.0, -3.0, err);
      forwardError(*err, __LINE__, -1.0);

      for (ii=0; ii<Nzbin; ii++) {
       for (jj=ii; jj<Nzbin; jj++) {

          if (self->tomo==tomo_auto_only && ii!=jj) continue;
          if (self->tomo==tomo_cross_only && ii==jj) continue;

          table = self->Pshear[idx_zz(ii,jj,Nzbin,err)]->table;
          intpar.self = self;
          intpar.i    = ii;
          intpar.j    = jj;
          for (i=0,slog=logsmin; i<N_s; i++,slog+=ds) {
             ss = exp(slog);
             intpar.r = ss;

             table[i] = 0.0;

             table[i] += prefactor * int_over_P_kappa(self, int_for_p_2, (void*)&intpar, err);
             forwardError(*err, __LINE__, -1.0);

             if (self->ia != ia_none && self->ia_terms != ia_only_II) {
                table[i] += -1.0 * self->A_ia * sqrt(prefactor) * ia_c1_rho_crit
                            * int_over_P_kappa(self, int_for_p_GI, (void*)&intpar, err);
                forwardError(*err, __LINE__, -1.0);
             }
             if (self->ia != ia_none && self->ia_terms != ia_only_GI) {
                table[i] += dsqr(self->A_ia * ia_c1_rho_crit)
                            * int_over_P_kappa(self, int_for_p_II, (void*)&intpar, err);
                forwardError(*err, __LINE__, -1.0);
             }

             testErrorRetVA(!finite(table[i]), ce_overflow, "Power spectrum P(l=%g)^{%d%d} not finite",
                            *err, __LINE__, -1.0, ss, ii, jj);

          }

       }
      }

   }

   testErrorRetVA(self->tomo==tomo_auto_only && i_bin!=j_bin, lensing_tomoij,
		  "Cross-correlation (bins # %d,%d) not valid if tomo=tomo_auto_only",
		  *err, __LINE__, -1.0, i_bin, j_bin);
   testErrorRetVA(self->tomo==tomo_cross_only && i_bin==j_bin, lensing_tomoij,
		  "Cross-correlation (bin #%d) not valid if tomo=tomo_auto_only",
		  *err, __LINE__, -1.0, i_bin);

   slog = log(s);
   i    = idx_zz(i_bin,j_bin,Nzbin, err);
   forwardError(*err, __LINE__, -1.0);
   f1 = interpol_wr(self->Pshear[i], slog, err);
   forwardError(*err, __LINE__, -1.0);


   /* Reduced-shear correction */
   if (self->reduced==reduced_K10) {
      f2 = Pg1(self, s, i_bin, j_bin, err);

      /* If error: return only P_kappa */
      forwardError(*err, __LINE__, f1);

      /* Add reduced-shear correction to P_kappa */
      f1 += f2;
   }

   return f1;
}

/* Used by Hankel transform functions, here and in halomodel. Returns Pshear^{ij}(l). *
 * If l outside of reduced-shear range, no error is forwarded.			      */
double P_projected_kappa(void *self, double l, int i_bin, int j_bin, error **err)
{
   double res;
   cosmo_lens *model;

   model = (cosmo_lens*)self;

   res = Pshear(model, l, i_bin, j_bin, err);

   /* If reduced-shear spectrum out of range: purge error. Return value is Pkappa. */
   if (model->reduced==reduced_K10 && getErrorValue(*err)==reduced_fourier_limit) {
      purgeError(err);
   }

   forwardError(*err, __LINE__, 0);

   return res;
}

double int_over_P_kappa(cosmo_lens *self, funcwithpars int_for_p, void *intpar, error **err)
{
   double f1, f2, a, da;


#ifndef fastxi

   *err = addError(lensing_fastxi, "The macro cosmo.h:fastxi is not defined", *err, __LINE__);
   return -1.0;

   /* Romberg-integration (slow) */
   if (self->cosmo->a_min<0.7) {
      f1 = sm2_qromberg(int_for_p, intpar, self->cosmo->a_min, 0.7, 1.0e-6, err);
      forwardError(*err, __LINE__, -1.0);
      f2 = sm2_qrombergo(int_for_p, intpar, 0.7, 1.0, sm2_midpntberg, 1.0e-7, err);
      forwardError(*err, __LINE__, -1.0);
   } else {
      f1 = 0.0;
      f2 = sm2_qrombergo(int_for_p, intpar, self->cosmo->a_min, 1.0, sm2_midpntberg, 1.0e-7, err);
      forwardError(*err, __LINE__, -1.0);
   }

#else

   /* Riemann summation (fast) */
   da = (1.0 - self->cosmo->a_min)/(self->cosmo->N_a-1.0);
   for (a=self->cosmo->a_min,f1=0.0; a<1.0; a+=da) {
     f1 += int_for_p(a, intpar, err);
     forwardError(*err, __LINE__, -1.0);
   }
   f1 = f1*da;
   f2 = 0.0;

#endif

   return f1+f2;
}

/* ============================================================ *
 * Functions for reduced-shear correction (K10).		*
 * ============================================================ */
const int parameter[M_PAR] = {p_dummy, p_Omegam, p_Omegade, p_w0de, p_Omegab, p_h100, p_sigma8, p_ns};


/* Fiducial cosmology for reduced-shear */
cosmo_lens *set_cosmological_parameters_lens_to_WMAP7(const redshift_t *nofz, tomo_t tomo, error **err)
{
   /* Parameters are:
      Om Od w0 w1 h Ob Onu Neffnu s8 ns
      nonlin transfer growth deparam norm amin
   */

   cosmo_lens *self;

   self = init_parameters_lens(0.27, 0.73, -1.0, 0.0, NULL, 0, 0.71, 0.045, 0.0, 0.0, 0.8, 0.96,
			       nofz->Nzbin, nofz->Nnz, nofz->nofz, nofz->par_nz,
			       smith03, eisenhu, growth_de, linder, norm_s8, tomo, reduced_none, 0.0,
			       ia_none, ia_undef, 0.0, err);
   forwardError(*err, __LINE__, NULL);

   return self;
}

/* Returns a pointed to the parameter corresponding to type par. *
 * Note: don't dereference the result before forwardError, in    *
 * case of error a NULL pointer is returned.			 */
double *par_to_pointer(cosmo *self, par_t par, error **err)
{
   double *p;

   switch (par) {
      case p_Omegam :
	 p = &(self->Omega_m);
	 break;
      case p_sigma8 :
	 p = &(self->normalization);
	 break;
      case p_Omegade :
	 p = &(self->Omega_de);
	 break;
      case p_w0de :
	 p = &(self->w0_de);
	 break;
      case p_h100 :
	 p = &(self->h_100);
	 break;
      case p_ns :
	 p = &(self->n_spec);
	 break;
      case p_Omegab :
	 p = &(self->Omega_b);
	 break;
      default :
	 *err = addErrorVA(math_unknown, "Unknown par_t %d", *err, __LINE__, par);
	 return NULL;
   }

   return p;
}

#define EPS 1.0e-6
void fill_dpar(cosmo *model, cosmo *wmap7, double *dpar, error **err)
{
   int alpha;
   double *tmp;

   testErrorRetVA(fabs(wmap7->w1_de-model->w1_de)>EPS, reduced_par,
		  "Parameter w1de (= %g) different from fiducial WMAP7 value (= %g)",
		  *err, __LINE__,, model->w1_de, wmap7->w1_de);

   testErrorRetVA(fabs(wmap7->Omega_nu_mass-model->Omega_nu_mass)>EPS, reduced_par,
		  "Parameter w1de (= %g) different from fiducial WMAP7 value (= %g)",
		  *err, __LINE__,, model->Omega_nu_mass, wmap7->Omega_nu_mass);

   testErrorRetVA(fabs(wmap7->Neff_nu_mass-model->Neff_nu_mass)>EPS, reduced_par,
		  "Parameter w1de (= %g) different from fiducial WMAP7 value (= %g)",
		  *err, __LINE__,, model->Neff_nu_mass, wmap7->Neff_nu_mass);

   dpar[0] = 0.0; /* Unused */
   for (alpha=1; alpha<M_PAR; alpha++) {

      tmp         = par_to_pointer(model, parameter[alpha], err);
      forwardError(*err, __LINE__,);
      dpar[alpha] = *tmp;

      tmp         = par_to_pointer(wmap7, parameter[alpha], err);
      forwardError(*err, __LINE__,);
      dpar[alpha] = dpar[alpha] - *tmp;

   }
}
#undef EPS

double Fbar(cosmo_lens *self, double a, int m_bin, int n_bin, error **err)
{
   double dwda, gg_m, gg_n, ww, fK, res;
   cosmoANDint ci;

   ci.self = self->cosmo;
   ci.i    = 0;

   dwda = R_HUBBLE*int_for_w(a, (void*)(&ci), err); forwardError(*err, __LINE__, 0.0);
   gg_m = G(self, a, m_bin, err);                   forwardError(*err, __LINE__, 0.0);
   gg_n = G(self, a, n_bin, err);                   forwardError(*err, __LINE__, 0.0);
   ww   = w(self->cosmo, a, 0, err);                forwardError(*err, __LINE__, 0.0);
   fK   = f_K(self->cosmo, ww, err);                forwardError(*err, __LINE__, 0.0);

   res = dwda/fK*a*a*gg_m*gg_n*(gg_n + gg_m)/2.0;

   return res;

}

void fill_Fbar_array(cosmo_lens *self, double *fbar, int m_bin, int n_bin, double amin, int N_a,
		     double da, error **err)
{
   int k;
   double a;

   for (k=0,a=amin; k<N_a; k++,a+=da) {
      fbar[k]  = Fbar(self, a, m_bin, n_bin, err);
      forwardError(*err, __LINE__,);
   }
}

void fill_dFbar_dp_array(cosmo_lens *self, par_t par, double *dfbar_dp, int m_bin, int n_bin, double amin,
			 int N_a, double da, error **err)
{
   const int pm[2] = {+1, -1};
   double *param, orig, *fmn[2], h;
   cosmo_lens *self_mod;
   int j, k;

   self_mod = copy_parameters_lens_only(self, err);           forwardError(*err, __LINE__,);

   param  = par_to_pointer(self_mod->cosmo, par, err);        forwardError(*err, __LINE__,);
   orig   = *param;
   h      = FH*orig;

   for (j=0; j<2; j++) {
      *param = orig + pm[j]*h;
      updateFrom_lens(self, self_mod, err);                   forwardError(*err, __LINE__,);
      fmn[j] = malloc_err(N_a*sizeof(double), err);           forwardError(*err, __LINE__,);
      fill_Fbar_array(self_mod, fmn[j], m_bin, n_bin, amin, N_a, da, err);
      forwardError(*err, __LINE__,);
   }

   for (k=0; k<N_a; k++) {
      dfbar_dp[k] = (fmn[0][k] - fmn[1][k])/(2.0*h);
   }

   free_parameters_lens(&self_mod);
   free(fmn[0]);
   free(fmn[1]);
}

/* ============================================================ *
 * Reduced-shear power spectrum correction (K10).		*
 * If q_mag_size!=0, magnification and size bias is added to    *
 * the reduced-shear result, where q = 2*(alpha+beta-1), and    *
 * alpha (beta) are the number density slopes with flux (size). *
 * ============================================================ */
#define EPS 1.0e-6
double Pg1(cosmo_lens *self, double s, int i_bin, int j_bin, error **err)
{
   int Nzbin, Nzcorr, Na, alpha, i, ii, jj;
   double logsmin, logsmax, dlogs, da, *fmn, *dfmn_dp[M_PAR], dpar[M_PAR], res, slog, *table;
   cosmo_lens *wmap7;   /* Fiducial model */


   testErrorRetVA(s<ELL_MIN_REDUCED || s>ELL_MAX_REDUCED, reduced_fourier_limit,
		  "Fourier scale %g out of range [%g;%g], setting Pg^(1)=0",
		  *err, __LINE__, 0.0, s, ELL_MIN_REDUCED, ELL_MAX_REDUCED);

   testErrorRetVA(self->tomo==tomo_auto_only && i_bin!=j_bin, lensing_tomoij,
		  "Cross-correlation (bins # %d,%d) not valid if tomo=tomo_auto_only",
		  *err, __LINE__, -1.0, i_bin, j_bin);
   testErrorRetVA(self->tomo==tomo_cross_only && i_bin==j_bin, lensing_tomoij,
		  "Cross-correlation (bin #%d) not valid if tomo=tomo_auto_only",
		  *err, __LINE__, -1.0, i_bin);
   testErrorRetVA(s<=0, ce_negative, "Negative or zero 2d Fourier scale ell=%g",* err, __LINE__, -1, s);
   testErrorRet(i_bin>j_bin, lensing_tomoij, "Pg_ij defined for i<=j", *err, __LINE__, -1);


   Nzbin  = self->redshift->Nzbin;
   Nzcorr = Nzbin*(Nzbin+1)/2;
   Na     = self->cosmo->N_a;
   if (self->Pg1==NULL) {

      logsmin   = log(ELL_MIN_REDUCED);
      logsmax   = log(ELL_MAX_REDUCED);
      dlogs     = (logsmax - logsmin)/((double)NELL_REDUCED - 1.0);
      da        = (1.0 - self->cosmo->a_min)/((double)Na); /* Not Na-1 to avoid a=1! */

      self->Pg1 = init_interTable_arr(Nzcorr, NELL_REDUCED, logsmin, logsmax, dlogs, 0.0, 0.0, err);

      fmn       = malloc_err(Na*sizeof(double), err);  forwardError(*err, __LINE__, 0.0);
      for (alpha=1; alpha<M_PAR; alpha++) {
         dfmn_dp[alpha] = calloc_err(Na, sizeof(double), err);  forwardError(*err, __LINE__, 0.0);
      }

      wmap7 = set_cosmological_parameters_lens_to_WMAP7(self->redshift, self->tomo, err);
      forwardError(*err, __LINE__, 0.0);
      fill_dpar(self->cosmo, wmap7->cosmo, (double *)dpar, err);
      forwardError(*err, __LINE__, 0.0);
      alpha = check_limits(dpar);
      testErrorRetVA(alpha!=0, reduced_limit, "Parameter dpar[%d]=%g out of range [%g;%g] for reduced-shear correction",
		     *err, __LINE__, 0.0, alpha, dpar[alpha], limits_lower[alpha], limits_upper[alpha]);

      for (ii=0; ii<Nzbin; ii++) {
	 for (jj=ii; jj<Nzbin; jj++) {

	    if (self->tomo==tomo_auto_only && ii!=jj) continue;
	    if (self->tomo==tomo_cross_only && ii==jj) continue;

	    table = self->Pg1[idx_zz(ii,jj,Nzbin, err)]->table;
	    forwardError(*err, __LINE__, -1.0);

	    /* F^{mn}, K10 eq. 10 */
	    fill_Fbar_array(wmap7, fmn, ii, jj, self->cosmo->a_min, Na, da, err);
	    forwardError(*err, __LINE__, 0.0);

	    /* Calculate dF^{mn}/dp for p!=p_0 (fiducial WMAP7 cosmology) */
	    for (alpha=1; alpha<M_PAR; alpha++) {
	       if (fabs(dpar[alpha])>EPS) {
		   fill_dFbar_dp_array(wmap7, parameter[alpha], dfmn_dp[alpha], ii, jj,
		  		      self->cosmo->a_min, Na, da, err);
		  forwardError(*err, __LINE__, 0.0);
	       }
	    }

	    for (i=0,slog=logsmin; i<NELL_REDUCED; i++,slog+=dlogs) {
	       res = sum_a_for_Pg1(slog, self->cosmo->a_min, Na, da, fmn, (const double **)dfmn_dp, dpar);
	       testErrorRetVA(!finite(res), math_negative, "Pg^1 (=%g) has to be positive",
			      *err, __LINE__, 0.0, res);
	       table[i] = res;
	    }

	 }
      }

      free_parameters_lens(&wmap7);
      free(fmn);
      for (alpha=1; alpha<M_PAR; alpha++) free(dfmn_dp[alpha]);
   }

   slog = log(s);
   i    = idx_zz(i_bin,j_bin,Nzbin, err);
   forwardError(*err, __LINE__, -1.0);
   res  = interpol_wr(self->Pg1[i], slog, err);
   forwardError(*err, __LINE__, -1.0);

   /* The factor 2 is already included in the fitting formulae */
   res = res + res/2.0*self->q_mag_size;

   return res;
}
#undef EPS


int change_xi(cosmo_lens* avant, cosmo_lens* apres)
{
   return change_Pshear(avant, apres);
}

/* ============================================================ *
 * Shear correlation function xi_+ (pm=+1) and xi_- (pm=-1).    *
 * ============================================================ */
double xi(cosmo_lens* self, int pm, double theta, int i_bin, int j_bin, error **err)
{
   double *table[2];
   double dlogtheta, logthetamin, logthetamax;
   double res;
   int Nzbin, Nzcorr, ii, jj, index;


   testErrorRetVA(theta<=0, ce_negative, "Negative angular scale theta=%g", *err, __LINE__, -1, theta);

   Nzbin  = self->redshift->Nzbin;
   Nzcorr = Nzbin*(Nzbin+1)/2;
   if (self->xiP==NULL || self->xiM==NULL) {
      self->xiP = init_interTable_arr(Nzcorr, N_thetaH, 0, 1, 1, 0.0, 0.0, err);
      forwardError(*err,__LINE__,0);
      self->xiM = init_interTable_arr(Nzcorr, N_thetaH, 0, 1, 1, 0.0, 0.0, err);
      forwardError(*err,__LINE__,0);

      for (ii=0; ii<Nzbin; ii++) {
	 for (jj=ii; jj<Nzbin; jj++) {

	    if (self->tomo==tomo_auto_only && ii!=jj) continue;
	    if (self->tomo==tomo_cross_only && ii==jj) continue;

	    index = idx_zz(ii,jj,Nzbin, err);
	    forwardError(*err, __LINE__, 0.0);
	    table[0] = self->xiP[index]->table;
	    table[1] = self->xiM[index]->table;
	    tpstat_via_hankel(self, table, &logthetamin, &logthetamax, tp_xipm, P_projected_kappa,
			      ii, jj, err);
	    forwardError(*err,__LINE__,0);

	    dlogtheta = (logthetamax-logthetamin)/((double)N_thetaH-1.0);
	    self->xiP[index]->a  = logthetamin;
	    self->xiM[index]->a  = logthetamin;
	    self->xiP[index]->b  = logthetamax;
	    self->xiM[index]->b  = logthetamax;
	    self->xiP[index]->dx = dlogtheta;
	    self->xiM[index]->dx = dlogtheta;

	 }
      }
   }

   testErrorRetVA(self->tomo==tomo_auto_only && i_bin!=j_bin, lensing_tomoij,
		  "Cross-correlation (bins # %d,%d) not valid if tomo=tomo_auto_only",
		  *err, __LINE__, -1.0, i_bin, j_bin);
   testErrorRetVA(self->tomo==tomo_cross_only && i_bin==j_bin, lensing_tomoij,
		  "Cross-correlation (bin #%d) not valid if tomo=tomo_auto_only",
		  *err, __LINE__, -1.0, i_bin);

   index = idx_zz(i_bin,j_bin,Nzbin, err);
      forwardError(*err,__LINE__,0);
   if (pm==1) {
      res = interpol_wr(self->xiP[index], log(theta), err);
      forwardError(*err,__LINE__,0);
   } else if (pm==-1) {
      res = interpol_wr(self->xiM[index], log(theta), err);
      forwardError(*err,__LINE__,0);
   } else {
      *err = addErrorVA(lensing_pm, "pm=%d not valid, has to be +1 or -1", *err, __LINE__, pm);
      res = -1.0;
      return res;
   }

   testErrorRetVA(self->reduced==reduced_K10 && pm==+1 && (theta<THETA_P_MIN_REDUCED || theta>THETA_MAX_REDUCED),
		  reduced_realsp_limit,
		  "Angular scale theta=%g' out of range [%g';%g'] for xi+ using reduced-shear correction",
		  *err, __LINE__, -1.0, theta/arcmin, THETA_P_MIN_REDUCED/arcmin, THETA_MAX_REDUCED/arcmin);
   testErrorRetVA(self->reduced==reduced_K10 && pm==-1 && (theta<THETA_M_MIN_REDUCED || theta>THETA_MAX_REDUCED),
		  reduced_realsp_limit,
		  "Angular scale theta=%g' out of range [%g';%g'] for xi- using reduced-shear correction",
		  *err, __LINE__, -1.0, theta/arcmin, THETA_M_MIN_REDUCED/arcmin, THETA_MAX_REDUCED/arcmin);

   return res;
}

int change_gamma2(cosmo_lens* avant, cosmo_lens* apres)
{
   return change_Pshear(avant,apres);
}

/* ============================================================ *
 * Tophat shear variance <|gamma|!2>.				*
 * ============================================================ */
double gamma2(cosmo_lens* self, double theta, int i_bin, int j_bin, error **err)
{
   double *table;
   double dlogtheta, logthetamin, logthetamax;
   double res;
   int Nzbin, Nzcorr, ii, jj, index;
   
   testErrorVA(theta<=0, ce_negative,"Negative scale theta=%g", *err, __LINE__, -1.0, theta);

   Nzbin  = self->redshift->Nzbin;
   Nzcorr = Nzbin*(Nzbin+1)/2;
   if (self->gamma==NULL) {

      self->gamma = init_interTable_arr(Nzcorr, N_thetaH, 0, 1, 1, 0.0, 0.0,err);
      forwardError(*err,__LINE__,0);

      for (ii=0; ii<Nzbin; ii++) {
	 for (jj=ii; jj<Nzbin; jj++) {

	    if (self->tomo==tomo_auto_only && ii!=jj) continue;
	    if (self->tomo==tomo_cross_only && ii==jj) continue;

	    index   = idx_zz(ii,jj,Nzbin,err);
	    forwardError(*err, __LINE__, 0.0);
	    table   = self->gamma[index]->table;
	    tpstat_via_hankel(self, &table, &logthetamin, &logthetamax, tp_gsqr, P_projected_kappa,
			      ii, jj, err);
	    forwardError(*err,__LINE__,0);
	    dlogtheta = (logthetamax-logthetamin)/((double)N_thetaH-1.0);
	    self->gamma[index]->a  = logthetamin;
	    self->gamma[index]->b  = logthetamax;
	    self->gamma[index]->dx = dlogtheta;

	 }
      }
   }

   testErrorRetVA(self->tomo==tomo_auto_only && i_bin!=j_bin, lensing_tomoij,
		  "Cross-correlation (bins # %d,%d) not valid if tomo=tomo_auto_only",
		  *err, __LINE__, -1.0, i_bin, j_bin);
   testErrorRetVA(self->tomo==tomo_cross_only && i_bin==j_bin, lensing_tomoij,
		  "Cross-correlation (bin #%d) not valid if tomo=tomo_auto_only",
		  *err, __LINE__, -1.0, i_bin);

   index = idx_zz(i_bin,j_bin,Nzbin,err);
   forwardError(*err,__LINE__,0);
   res   = interpol_wr(self->gamma[index], log(theta), err);
   forwardError(*err,__LINE__,0);
   return res;

}

int change_map2(cosmo_lens* avant, cosmo_lens* apres)
{
   return change_Pshear(avant,apres);
}

/* ============================================================ *
 * Aperture mass variance <M_ap^2> with polynomial filter.	*
 * ============================================================ */
double map2(cosmo_lens* self, double theta, int i_bin, int j_bin, error **err)
{
  double res;
  res = map2_poly(self, theta, i_bin, j_bin, err);
  forwardError(*err, __LINE__, 0);
  return res;
}


double map2_poly(cosmo_lens* self, double theta,  int i_bin, int j_bin, error **err)
{
   double *table;
   double dlogtheta, logthetamin, logthetamax;
   double res;
   int Nzbin, Nzcorr, ii, jj, index;


   testErrorRetVA(theta<=0, ce_negative, "Negative angular scale theta=%g", *err, __LINE__, -1, theta);

   Nzbin = self->redshift->Nzbin;
   Nzcorr = Nzbin*(Nzbin+1)/2;
   if (self->map_poly == NULL) {
      self->map_poly = init_interTable_arr(Nzcorr, N_thetaH, 0, 1, 1, 0.0, 0.0, err);
      forwardError(*err, __LINE__, 0);
      for (ii=0; ii<Nzbin; ii++) {
	 for (jj=ii; jj<Nzbin; jj++) {

	    if (self->tomo==tomo_auto_only && ii!=jj) continue;
	    if (self->tomo==tomo_cross_only && ii==jj) continue;

	    index   = idx_zz(ii,jj,Nzbin,err);
	    forwardError(*err,__LINE__,0);
	    table   = self->map_poly[index]->table;
	    tpstat_via_hankel(self, &table, &logthetamin, &logthetamax, tp_map2_poly, P_projected_kappa,
			      ii, jj, err);
	    forwardError(*err, __LINE__, 0);
	    dlogtheta = (logthetamax-logthetamin)/((double)N_thetaH-1.0);
	    self->map_poly[index]->a = logthetamin;
	    self->map_poly[index]->b = logthetamax;
	    self->map_poly[index]->dx = dlogtheta;

	 }
      }
   }

   testErrorRetVA(self->tomo==tomo_auto_only && i_bin!=j_bin, lensing_tomoij,
		  "Cross-correlation (bins # %d,%d) not valid if tomo=tomo_auto_only",
		  *err, __LINE__, -1.0, i_bin, j_bin);
   testErrorRetVA(self->tomo==tomo_cross_only && i_bin==j_bin, lensing_tomoij,
		  "Cross-correlation (bin #%d) not valid if tomo=tomo_auto_only",
		  *err, __LINE__, -1.0, i_bin);

   index = idx_zz(i_bin,j_bin,Nzbin, err);
   forwardError(*err,__LINE__,0);
   res   = interpol_wr(self->map_poly[index], log(theta), err);
   forwardError(*err, __LINE__, 0);

   return res;
}

double map2_gauss(cosmo_lens* self, double theta, int i_bin, int j_bin, error **err)
{
   double *table;
   double dlogtheta, logthetamin, logthetamax;
   double res;
   int Nzbin, Nzcorr, ii, jj, index;

   testErrorRetVA(theta<=0, ce_negative, "Negative angular scale theta=%g", *err, __LINE__, -1, theta);

   Nzbin = self->redshift->Nzbin;
   Nzcorr = Nzbin*(Nzbin+1)/2;
   if (self->map_gauss==NULL) {
      self->map_gauss = init_interTable_arr(Nzcorr, N_thetaH, 0, 1, 1, 0.0, 0.0, err);
      forwardError(*err, __LINE__, 0);
      for (ii=0; ii<Nzbin; ii++) {
	 for (jj=ii; jj<Nzbin; jj++) {

	    if (self->tomo==tomo_auto_only && ii!=jj) continue;
	    if (self->tomo==tomo_cross_only && ii==jj) continue;

	    index   = idx_zz(ii,jj,Nzbin,err);
	    forwardError(*err,__LINE__,0);
	    table   = self->map_gauss[index]->table;
	    tpstat_via_hankel(self, &table, &logthetamin, &logthetamax, tp_map2_gauss, P_projected_kappa,
			      ii, jj, err);
	    forwardError(*err, __LINE__, 0);
	    dlogtheta = (logthetamax-logthetamin)/((double)N_thetaH-1.0);
	    self->map_gauss[index]->a = logthetamin;
	    self->map_gauss[index]->b = logthetamax;
	    self->map_gauss[index]->dx = dlogtheta;
	 }
      }
   }

   testErrorRetVA(self->tomo==tomo_auto_only && i_bin!=j_bin, lensing_tomoij,
		  "Cross-correlation (bins # %d,%d) not valid if tomo=tomo_auto_only",
		  *err, __LINE__, -1.0, i_bin, j_bin);
   testErrorRetVA(self->tomo==tomo_cross_only && i_bin==j_bin, lensing_tomoij,
		  "Cross-correlation (bin #%d) not valid if tomo=tomo_auto_only",
		  *err, __LINE__, -1.0, i_bin);

   index = idx_zz(i_bin,j_bin,Nzbin,err);
   forwardError(*err, __LINE__, 0);
   res = interpol_wr(self->map_gauss[index], log(theta), err);
   forwardError(*err, __LINE__, 0);

   return res;
}

/* Returns the generalised ring statistic using a as the coefficients *
 * of the T_+ decomposition, see FK09 (4, 8, 11)		      *
 * If a==NULL, Z+ from SK07 is used.				      */
double RR(cosmo_lens *lens, double THETA_MIN, double THETA_MAX, const double *a, int N,
	  poly_t poly, int pm, error **err)
{
   double theta, dtheta, res[2]={0,0}, A, B, x, summand;

   testErrorRetVA(abs(pm)!=1, mr_incompatible, "pm=%d not valid, has to be +1 or -1",
		  *err, __LINE__, 0.0, pm);

   if (a==NULL) {
#ifdef __MRING_H
      /* Using the (non-optimised) Schneider&Kilbinger (2007) functions Z+,Z- */
      R_from_xi(lens, THETA_MIN/THETA_MAX, THETA_MAX, res, res+1, err);
      forwardError(*err, __LINE__, 0.0);
      return res[0];
#else
      *err = addError(mr_null, "Coefficients a=NULL", *err, __LINE__);
      return 0.0;
#endif
   }


   /* Decrease step size dtheta for more precision,     *
    * in particular, to make (numerical) B-mode smaller	*/
   dtheta  = 1.0*arcsec/10.0;
   if (THETA_MAX>30*arcmin) dtheta *= (THETA_MAX/(30*arcmin));


   /* FK09 (8) */
   A = (THETA_MAX - THETA_MIN)/2.0;
   B = (THETA_MAX + THETA_MIN)/2.0;

   for (theta=THETA_MIN,res[0]=0.0; theta<=THETA_MAX; theta+=dtheta) {

      summand  = xi(lens, pm, theta, 0, 0, err)*theta;
      forwardError(*err, __LINE__, -1.0);

     switch (poly) {
	case cheby : case cheby2:
	   x        = (theta-B)/A;
	   if (pm==+1) {
	      summand *= Tp(x, a, N, poly, err);
	      forwardError(*err, __LINE__, -1.0);
	   } else if (pm==-1) {
	      summand *= Tm(x, a, N, poly, B/A, err);
	      forwardError(*err, __LINE__, -1.0);
	   }
	   summand  /= dsqr(THETA_MAX);
	   break;
	default :
	   *err = addErrorVA(ce_unknown, "Wrong poly_t (%d)", *err, __LINE__, poly);
	   return -1;
     }

     res[0] += summand;

   }

   testErrorRet(!finite(res[0]), math_infnan, "R_E is not finite", *err, __LINE__, -1.0);

   return res[0]*dtheta;
}

double E_cosebi(cosmo_lens *lens, int n, double Psimin, double Psimax, int i_bin, int j_bin,
		const char *path, double *B_cosebi, error **err)
{
   double rp, rm;

   testErrorRetVA(n > NMAX_COSEBI, mr_range, "COSEBI mode n=%d cannot be larger than NMAX_COSEBI=%d",
		  *err, __LINE__, -1.0, n, NMAX_COSEBI);

   if (lens->c_cosebi == NULL) {
      /* Read COSEBI zeros and normalisation from file, calculate polynomial coefficients */
      lens->c_cosebi = read_zeros_norm_cosebi_auto_check(Psimin, Psimax, path, err);
      forwardError(*err, __LINE__, 0.0);
   }

   rp = RR_cosebi(lens, Psimin, Psimax, i_bin, j_bin, n, +1, err);
   forwardError(*err, __LINE__, 0.0);

   rm = RR_cosebi(lens, Psimin, Psimax, i_bin, j_bin, n, -1, err);
   forwardError(*err, __LINE__, 0.0);

   if (B_cosebi != NULL) {
      *B_cosebi = 0.5 * (rp - rm);
   }

   return 0.5 * (rp + rm);
}

/* Returns the COSEBI of order N using a as the coefficients *
 * of the T_+ decomposition, see SKE10 (36 - 38)	     */
double RR_cosebi(cosmo_lens *lens, double THETA_MIN, double THETA_MAX, int i_bin, int j_bin,
		 int n, int pm, error **err)
{
   double res, zmax, z, dz;
   int Nz;
   cosmo_lensANDextra intpar;

   testErrorRetVA(abs(pm)!=1, mr_incompatible, "pm=%d not valid, has to be +1 or -1",
		  *err, __LINE__, 0.0, pm);
   zmax = log(THETA_MAX / THETA_MIN);

   intpar.self  = lens;
   intpar.err   = err;
   intpar.pm    = pm;
   intpar.i_bin = i_bin;
   intpar.j_bin = j_bin;
   intpar.n     = n;
   intpar.thmin = THETA_MIN;


   // Romberg integration. Up to N=13, B-mode is < 6e-16. E_13 = 2e-14
   Nz = 10;
   dz = zmax/(double)Nz;
   for (z=0.0,res=0.0; z<=zmax-dz; z+=dz) {
      res += sm2_qromberg(dRR_cosebi_dz, (void*)&intpar, z, z+dz, 1.0e-6, err);
   }
   forwardError(*err, __LINE__, -1.0);


   // Adaptive Gaussian integration
   /*
   int n = 1000;
   size_t neval;
   double eps_abs, eps_rel, res_err, tmp;
   gsl_function F;
   gsl_integration_workspace *w = gsl_integration_workspace_alloc(n);
   F.function = &dRR_cosebi_dz;
   F.params   = &intpar;
   eps_abs    = 1.0e-5;
   eps_rel    = 1.0e-3;

   Nz = 1;
   dz = zmax / (double)Nz;
   for (z=0.0,res=0.0; z<=zmax-dz; z+=dz) {
      gsl_integration_qag(&F, z, z+dz, eps_abs, eps_rel, n, GSL_INTEG_GAUSS41, w, &tmp, &res_err);
      //gsl_integration_qng (&F, z, z+dz, eps_abs, eps_rel, &tmp, &res_err, &neval);
      res += tmp;
   }
   gsl_integration_workspace_free(w);
   */

   // Monte-Carlo integration
   /*
   gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(1);
   gsl_monte_vegas_init(s);

   gsl_rng *rng;
   gsl_monte_function F;
   double min[0], max[0], res_err;
   size_t calls;

   calls    = 10000*N;
   rng      = gsl_rng_alloc(gsl_rng_default);
   F.f      = &dRR_cosebi_dz_MC;
   F.dim    = 1;
   F.params = &intpar;
   min[0]   = 0.0;
   max[0]   = zmax;

   gsl_monte_vegas_integrate(&F, min, max, 1, calls, rng, s, &res, &res_err);

   gsl_monte_vegas_free(s);
   */

   testErrorRet(!finite(res), math_infnan, "R is not finite", *err, __LINE__, -1.0);

   res *= THETA_MIN*THETA_MIN;

   return res;
}

/*
double dRR_cosebi_dz_MC(double *z, int ndim, void *intpar)
{
}
*/

double dRR_cosebi_dz(double z, void *intpar, error **err)
{
   double theta, summand;
   int i_bin, j_bin, pm, n;
   cosmo_lensANDextra *extra;
   cosmo_lens *lens;
   //error **err;

   extra = (cosmo_lensANDextra*)intpar;
   lens  = extra->self;
   pm    = extra->pm;
   i_bin = extra->i_bin;
   j_bin = extra->j_bin;
   n     = extra->n;
   theta = extra->thmin * exp(z);

   summand = xi(lens, pm, theta, i_bin, j_bin, err);
   forwardError(*err, __LINE__, -1.0);

   //printf("%g %g %d %d\n", theta/arcmin, summand, n, pm);

   if (pm==+1) {
      summand *= Tplog_c(z, lens->c_cosebi, n, err);
      forwardError(*err, __LINE__, -1.0);
   } else if (pm==-1) {
      summand *= Tmlog(z, lens->c_cosebi, n, err);
      forwardError(*err, __LINE__, -1.0);
   }

   summand *= exp(2.0*z);

   return summand;
}


double int_for_map2_slow(double logell, void *intpar, error **err)
{
   cosmo_lensANDiiid *cplusplus;
   double theta, etasqr, res, ell;
   tpstat_t tpstat;
   int i_bin, j_bin;

   ell       = exp(logell);
   cplusplus = (cosmo_lensANDiiid*)intpar;
   theta     = cplusplus->r;
   etasqr    = dsqr(theta*ell);
   tpstat    = (tpstat_t)cplusplus->t;
   i_bin     = cplusplus->i;
   j_bin     = cplusplus->j;

   res = dsqr(ell)*Pshear(cplusplus->self, ell, i_bin, j_bin, err);
   forwardError(*err, __LINE__, 0);
   if (tpstat==tp_map2_poly) {
      /* TODO: implement Bessel function J_4 */
      *err = addError(ce_unknown, "wrong tpstat", *err, __LINE__);
      return -1;
   } else if (tpstat==tp_map2_gauss) {
      res *= dsqr(etasqr)/4.0*exp(-etasqr);
   } else {
      *err = addError(ce_unknown, "wrong tpstat", *err, __LINE__);
      return -1;
   }

   return res;
}

double map2_slow(cosmo_lens *self, double theta, tpstat_t tpstat, int i_bin, int j_bin, error **err)
{
   cosmo_lensANDiiid intpar;
   double res;

   intpar.self = self;
   intpar.t    = (int)tpstat;
   intpar.r    = theta;
   intpar.i    = i_bin;
   intpar.j    = j_bin;

   res = 1.0/(2.0*pi)*sm2_qromberg(int_for_map2_slow, (void*)&intpar, log(s_min), log(s_max),
				   1.0e-6, err);
   forwardError(*err, __LINE__, 0);

   return res;
}

/* ============================================================ *
 * Data input, likelihood stuff (used to be in likeli.c).	*
 * ============================================================ */


/* ============================================================ *
 * Useful if a datcov structure only contains the covariance.   *
 * ============================================================ */
datcov *init_datcov_for_cov_only(int Nzbin, int Ntheta, error **err)
{
   datcov *dc;

   dc = malloc_err(sizeof(datcov), err);         forwardError(*err, __LINE__, NULL);
   dc->Nzbin    = Nzbin;
   dc->Nzcorr   = Nzbin*(Nzbin+1)/2;
   dc->Ntheta   = Ntheta;
   dc->Ntheta2  = 0;
   dc->n        = dc->Nzcorr*dc->Ntheta;
   dc->usecov   = 1;
   dc->lndetC   = 0.0;
   dc->theta    = dc->data = dc->var = dc->cov[0] = dc->cov[1] = dc->cov[2] = NULL;
   dc->cov_scaling = cov_const;

   return dc;
}

datcov *init_data_cov_tomo(char* dataname, char *dataname2, char** covname_ptr, lensdata_t type,
			   decomp_eb_filter_t decomp_eb_filter, lensformat_t format, double corr_invcov, 
			   double a1, double a2, order_t order,
			   cov_scaling_t cov_scaling, error **err)
{
   datcov *res;
   gsl_matrix_view A;
   
   testErrorRetVA(type<0||type>=Nlensdata_t, ce_unknown, "Unknown lens data type %d",
		  *err, __LINE__, NULL, type);

   res = malloc_err(sizeof(datcov), err);           forwardError(*err,__LINE__,NULL);

   res->format = format;
   res->usecov = 1;
   res->type   = type;
   res->order  = order;
   res->decomp_eb_filter = decomp_eb_filter;

   /* a1 and a2 are zero if format!=angle_wquadr */
   res->a1 = a1;
   res->a2 = a2;

   switch (res->type) {

      case map3gauss:
         /* Three angular scales per bin */
         read_data_3rd(res, dataname, 0, err);
         forwardError(*err, __LINE__, NULL);
         break;

      case map2gauss_map3gauss_diag: case map2gauss_map3gauss:
      case decomp_eb_map3gauss_diag: case decomp_eb_map3gauss:
         /* Two files */
         read_data_2nd_3rd(res, dataname, dataname2, err);
         forwardError(*err, __LINE__, NULL);
         break;

      default:
         /* One angular scale per bin */
         read_data_tomo(res, dataname, 0, order, err);
         forwardError(*err,__LINE__,NULL);
         break;
   }      

   read_cov_tomo(res, covname_ptr[0], 0, err);                      forwardError(*err,__LINE__,NULL);

   if (cov_scaling == cov_ESH09) {

      read_cov_tomo(res, covname_ptr[1], 1, err);                   forwardError(*err,__LINE__,NULL);
      read_cov_tomo(res, covname_ptr[2], 2, err);                   forwardError(*err,__LINE__,NULL);

      res->lndetC = -1.0; /* dummy value; recalculated in chi2_lensing */

   } else if (cov_scaling == cov_const) {

      gsl_set_error_handler_off();
      A = gsl_matrix_view_array(res->cov[0], res->n, res->n);

      /* Replace res->cov with L, where C = L L^T */
      if (gsl_linalg_cholesky_decomp(&A.matrix) == GSL_EDOM) {
         del_data_cov(&res);
         *err = addError(mv_cholesky, "Cholesky decomposition failed", *err, __LINE__);
         return NULL;
      }


      /* Now cov is L, where L * L^T is the covariance matrix.           *
       * Therefore, multiply with the inverse square root Hartlap factor */
      multiply_all(res->cov[0], res->n*res->n, sqrt(1.0/corr_invcov));

      /* Determinant of the covariance.    *
       * A contains L with C = L^T L, so   *
       * ln|C| = ln|L|^2 = 2 ln|L| .       */
      res->lndetC = 2.0 * ln_determinant(res->cov[0], res->n, err);
      forwardError(*err,__LINE__,NULL);


   } else {

      *err = addErrorVA(ce_unknown, "Unknown or invalid cov_scaling type %d(%s)",
                              *err, __LINE__, cov_scaling, scov_scaling_t(cov_scaling));
      return NULL;

   }

   res->cov_scaling = cov_scaling;

   return res;
}

void del_data_cov(datcov** dc)
{
   datcov *d;
   d = *dc;

   if (d->theta) free(d->theta);
   if (d->data) free(d->data);
   if (d->usecov) {
     free(d->cov[0]);
     if (d->cov[1] != NULL) free(d->cov[1]);
     if (d->cov[2] != NULL) free(d->cov[2]);
   } else {
      free(d->var);
   }
   free(d);
   d = NULL;
}

/* ============================================================ *
 * Reads data file into structure dc (has to be initialised     *
 * before). The file has to have the structure                  *
 *   theta[arcmin] c_00 c_01 ... c_0{n-1} c_11 ... c_{n-1}{n-1} *
 * with c_ij some (correlation) quantity corr. to bins i, j.    *
 * If Nzbin>0 on input: performs a consistency check for second *
 *  order.							*
 * The file can contain xi+ and xi- data, in this case it       *
 * corresponds to 'cat xi+.dat xi-.dat > xi+-.dat'.		*
 * Used for 2nd-order and 3rd-order diagonal, also for combined *
 * 2nd and 3rd.							*
 * ============================================================ */
void read_data_tomo(datcov *dc, char data_name[], int Nzbin, order_t order, error **err)
{
   int n, i, Ncol, offset, nz;
   size_t Nrec, Ntheta, Nzcorr;
   double dNcol;
   double *ptr;

   /* Read file to double pointer */
   Nrec = 0;
   ptr = (double*)read_any_list_count(data_name, &Nrec, "%lg", sizeof(double), &Ntheta, err);
   forwardError(*err, __LINE__,);

   /* Number of lines and columns */
   dNcol = (double)Nrec/(double)Ntheta;
   Ncol = Nrec/Ntheta;
   testErrorRet(fabs(dNcol-(double)Ncol)>EPSILON, lensing_inconsistent, "Columns do not have equal length",
		*err, __LINE__,);

   if (dc->format==angle_mean || dc->format==angle_wlinear || dc->format==angle_wquadr) {
      Nzcorr = Ncol - 2;
      testErrorRetVA(Nzcorr<1, lensing_inconsistent, "Not enough columns in file for lensformat %s(%d)",
		     *err, __LINE__,, slensformat_t(dc->format), dc->format);
   } else {
      Nzcorr = Ncol - 1;
   }

   /* Fill datcov structure */
   dc->Ntheta  = (int)Ntheta;
   dc->Ntheta2 = 0;
   if (order == third_order) {
      dc->Nzbin   = Nperm_to_Ntheta(Nzcorr, err);                 forwardError(*err, __LINE__,);
   } else {
      dc->Nzbin   = get_and_check_Nzbin(Nzcorr, Nzbin, err);      forwardError(*err, __LINE__,);
   }
   dc->Nzcorr  = Nzcorr;
   dc->n       = Ntheta*Nzcorr;

   dc->data  = malloc_err(sizeof(double)*dc->n, err);          forwardError(*err, __LINE__,);
   dc->theta = malloc_err(sizeof(double)*Ntheta, err);         forwardError(*err, __LINE__,);
   if (dc->format==angle_mean || dc->format==angle_wlinear || dc->format==angle_wquadr) {
      dc->theta2 = malloc_err(sizeof(double)*Ntheta, err);     forwardError(*err, __LINE__,);
   } else {
      dc->theta2 = NULL;
   }
   dc->var   = NULL;

   /* First column: angular scale */
   for (i=0; i<Ntheta; i++) {
      dc->theta[i] = ptr[i*Ncol+0]*arcmin;
      if (dc->format==angle_mean || dc->format==angle_wlinear || dc->format==angle_wquadr) {
         /* Second column: upper bin limit */
         dc->theta2[i] = ptr[i*Ncol+1]*arcmin;
      }
   }

   if (dc->format==angle_mean || dc->format==angle_wlinear || dc->format==angle_wquadr) {
      offset = 2;
   } else {
      offset = 1;
   }

   /* In final data vector: angular scale is varying fast, redshift is varying slow */
   for (nz=0,n=0; nz<Nzcorr; nz++) {
      for (i=0; i<Ntheta; i++,n++) {
         testErrorRetVA(n>=dc->n, math_overflow, "Index overflow (%d>=%d)", *err, __LINE__,, n, dc->n);

         dc->data[n] = ptr[i*Ncol+offset+nz];
      }
   }
}

/* ============================================================ *
 * Reada the covariance matrix in block format.			*
 * ============================================================ */
void read_cov_tomo(datcov* dc, char cov_name[], int icov, error **err)
{
   int i, j, index, Nzcorr;
   size_t Nrec, Nrow, Nangular, Nangular2;
   double *ptr;

   testErrorRetVA(icov<0 || icov>2, lensing_range, "Cov matrix number %d out of range [0;2]",
		  *err, __LINE__,, icov);

   testErrorRet(dc->Ntheta<=0, lensing_initialised,
		"Data file has to be read before covariance, to set number of bins (Ntheta)", *err, __LINE__,);
   testErrorRet(dc->Nzcorr<=0, lensing_initialised,
		"Data file has to be read before covariance, to set number of bins (Nzcorr)", *err, __LINE__,);
   testErrorRet(dc->Nzbin<=0, lensing_initialised,
		"Data file has to be read before covariance, to set number of bins (Nzbin)", *err, __LINE__,);

   /* Read file to double pointer */
   Nrec = 0;
   ptr  = (double*)read_any_list_count(cov_name, &Nrec, "%lg", sizeof(double), &Nrow, err);
   forwardError(*err, __LINE__,);

   testErrorVA(Nrow!=dc->n, lensing_inconsistent,
	       "Covariance matrix (%d rows) inconsistent with data vector of length %d",
	       *err, __LINE__, Nrow, dc->n);
   testErrorRetVA(Nrec!=Nrow*Nrow, lensing_inconsistent, "Covariance matrix is not square, Nrow=%d, Nrec=%d\n",
		  *err, __LINE__,, Nrec, Nrow);

   /* The tests above ensure that different covariance files are consistent */

   dc->cov[icov] = malloc_err(dc->n*dc->n*sizeof(double), err);   forwardError(*err, __LINE__,);

   switch (dc->type) {
      case map3gauss :
         Nangular  = dc->Ntheta * (dc->Ntheta + 1) * (dc->Ntheta + 2) / 6;
         Nangular2 = 0;
         break;
      case map2gauss_map3gauss : case decomp_eb_map3gauss :
         Nangular  = dc->Ntheta;
         Nangular2 = dc->Ntheta2 * (dc->Ntheta2 + 1) * (dc->Ntheta2 + 2) / 6;
         break;
      case map2gauss_map3gauss_diag : case decomp_eb_map3gauss_diag :
         Nangular  = dc->Ntheta;
         Nangular2 = dc->Ntheta2;
         break;
      default :
         Nangular  = dc->Ntheta;
         Nangular2 = 0;
         break;
   }

   Nzcorr  = dc->Nzcorr;

   /* MKDEBUG: Second term had Nzcorr2 instead of Nzcorr */
   testErrorRetVA(dc->n != Nangular * Nzcorr + Nangular2 * Nzcorr,
		  lensing_range, "Inconsistent number of angular/redshift bins (%d != %d*%d + %d*%d",
		  *err, __LINE__,, dc->n, Nangular, Nzcorr, Nangular2, Nzcorr);


   /* MKDEBUG: New! Re-ordering of loop. Takes into account 2nd and 3rd. *
    * 16 Nov 2012. */

   for (i=0,index=0; i<dc->n; i++) {
      for (j=0; j<dc->n; j++,index++) {

         testErrorRetVA(index >= dc->n * dc->n, lensing_range, "Index overflow, %d >= %d",
			*err, __LINE__,, index, dc->n * dc->n);
         dc->cov[icov][index] = ptr[index];

      }
   }

}

/* ============================================================ *
 * Creates the vectors xip, xim, theta and theta2 (if required  *
 * by the lens format), all of size N,				*
 * and copies the content of datcov to those vectors.		*
 * ============================================================ */
void datcov2xipm(const datcov *dc, int i_bin, int j_bin, double **xip, double **xim, double **theta,
		 double **theta2, int *N, error **err)
{
   int i, index, offset;

   testErrorRetVA(dc->type!=xipm, lensing_type, "lenstype has to be %d('%s'), not %d",
   		  *err, __LINE__,, xipm, slensdata_t(xipm), dc->type);

   testErrorRetVA(i_bin>j_bin, lensing_tomoij, "i_bin(=%d) has to be smaller or equal j_bin(=%d)",
		  *err, __LINE__,, i_bin, j_bin);

   testErrorRetVA(j_bin>=dc->Nzbin, lensing_tomoij, "j_bin(=%d) has to be smaller than Nzbin=%d\n",
   		  *err, __LINE__,, j_bin, dc->Nzbin);

   *N = dc->Ntheta/2;

   *xip   = malloc_err(*N*sizeof(double), err);  forwardError(*err, __LINE__,);
   *xim   = malloc_err(*N*sizeof(double), err);  forwardError(*err, __LINE__,);
   *theta = malloc_err(*N*sizeof(double), err);  forwardError(*err, __LINE__,);

   if (dc->format==angle_mean || dc->format==angle_wlinear || dc->format==angle_wquadr) {
      *theta2 = malloc_err(*N*sizeof(double), err);  forwardError(*err, __LINE__,);
   } else {
      *theta2 = NULL;
   }


   offset = dc->Ntheta * idx_zz(i_bin, j_bin, dc->Nzbin, err);
   forwardError(*err, __LINE__,);
   for (i=0; i<*N; i++) {
      index = offset + i;
      (*xip)[i]   = dc->data[index];
      (*theta)[i] = dc->theta[i];
      if (dc->format==angle_mean || dc->format==angle_wlinear || dc->format==angle_wquadr) {
	 (*theta2)[i] = dc->theta2[i];
      }
   }
   
   for (i=*N; i<dc->Ntheta; i++) {
      index = offset + i;
      (*xim)[i-(*N)] = dc->data[index];
   }

}

/* Really only reads a covariance (in column format) for xi+ */
void read_cov_col(datcov *dc, char cov_name[], error **err)
{
   FILE *F;
   int i, j, res, nn;
   unsigned int ncomment;
   double dummy[5];

   nn          = numberoflines_comments(cov_name, &ncomment, err); forwardError(*err, __LINE__,);
   dc->n       = (int)sqrt(nn);
   dc->Ntheta  = dc->n;
   dc->Ntheta2 = 0;
   dc->Nzbin   = dc->Nzcorr = 1;
   F           = fopen_err(cov_name, "r", err);                 forwardError(*err, __LINE__,);
   dc->cov[0]  = malloc_err(sizeof(double)*dc->n*dc->n, err);   forwardError(*err, __LINE__,);
   dc->cov[1]  = dc->cov[2] = NULL;
   dc->theta   = malloc_err(sizeof(double)*dc->n, err);         forwardError(*err, __LINE__,);

   for (i=0; i<dc->n; i++) {
      for (j=0; j<dc->n; j++) {
	 res = fscanf(F, "%lf %lf %lf %lf %lf\n", dummy, dummy+1, dummy+2, dummy+3, dummy+4);
	 testErrorRet(res==EOF, io_eof, "Premature eof", *err, __LINE__,);
	 if (i==0) {
	    dc->theta[j] = dummy[1];
	 }
	 dc->cov[0][i*dc->n+j] = dummy[2];
      }
   }

}


/* ============================================================ *
 * The following functions implement the Eifler, Schneider &    *
 * Hartlap 2009 (ESH09) cosmology-dependent covariance.		*
 * ============================================================ */

/* ============================================================ *
 * Scales the cosmic variance xi+/xi- covariance term, see      *
 * ESH09, eq. (20, 21).						*
 * ============================================================ */
void scale_cosmic_variance_ESH09(cosmo_lens *model, gsl_matrix *cov, const datcov *dc, error **err)
{
   int pm1, pm2, i, j, i_bin, j_bin, k_bin, l_bin, ci, cj;
   double f, xi1, xi2, xi1_fid, xi2_fid, old;

   for (i_bin=ci=0; i_bin<dc->Nzbin; i_bin++) {
      for (j_bin=i_bin; j_bin<dc->Nzbin; j_bin++) {
	 for (i=0; i<dc->Ntheta; i++,ci++) {

	    pm1  = get_pm(dc->type, i, dc->Ntheta, err);
       	    forwardError(*err, __LINE__,);
	    xi1  = xi(model, pm1, dc->theta[i], i_bin, j_bin, err);
       	    forwardError(*err, __LINE__,);
	    xi1_fid = xi(dc->fiducial, pm1, dc->theta[i], i_bin, j_bin, err);
	    forwardError(*err, __LINE__,);

	    for (k_bin=cj=0; k_bin<dc->Nzbin; k_bin++) {
	       for (l_bin=k_bin; l_bin<dc->Nzbin; l_bin++) {
		  for (j=0;j<dc->Ntheta;j++,cj++) {

		     pm2  = get_pm(dc->type, j, dc->Ntheta, err);
                     forwardError(*err, __LINE__,);
		     xi2  = xi(model, pm2, dc->theta[j], k_bin, l_bin, err);
                     forwardError(*err, __LINE__,);
		     xi2_fid  = xi(dc->fiducial, pm2, dc->theta[j], k_bin, l_bin, err);
          	     forwardError(*err, __LINE__,);

		     f   = xi1 / xi1_fid * xi2 / xi2_fid;
		     old = gsl_matrix_get(cov, ci, cj);
		     gsl_matrix_set(cov, ci, cj, old * f);

		  }
	       }
	    }

	 }
      }
   }
}

/* ============================================================ *
 * Returns +-1, corresponding to type (xip, xim) and angular    *
 * index (for xipm).						*
 * ============================================================ */
int get_pm(lensdata_t type, int i, int Ntheta, error **err)
{
   /* First half of data vector = xi+, second half = xi- */
   if (type == xipm) {
      if (i < Ntheta/2) return +1;
      else return -1;
   }
   if (type == xip) return +1;
   if (type == xim) return -1;

   *err = addErrorVA(lensing_type, "Invalid lens type %d", *err, __LINE__, type);
   return 0;
}

/* ============================================================ *
 * Returns index c such that list[c] is the closest entry for   *
 * x. prev (if >0) is starting index for search.		*
 * ============================================================ */
int find_bin(double x, const double *list, int N, int prev, error **err)
{
   int c;

   if (prev>=0) {
      /* Start search at previous index */
      for (c=prev; c<N-1; c++) {
	 if (list[c] <= x && list[c+1] > x) {
	    /* Nearer to lower or upper bin corner? */
	    if (x-list[c] < (list[c+1]-list[c])/2.0) return c;
	    return c + 1;
	 }
      }
   }

   /* x outside of range? -> constant extrapolation */
   if (list[N-1] < x) return N-1;
   if (x < list[0]) return 0;

   /* Start search from beginning */
   for (c=0; c<N-1; c++) {
      if (list[c] <= x && list[c+1] > x) {
	 if (x-list[c] < list[c+1]-list[c]) return c;
	 return c + 1;
      }
   }

   *err = addErrorVA(lensing_range, "Invalid range for bin search, th=%g", *err, __LINE__, x);
   return -1;
}

/* ============================================================ *
 * Scales the mixed xi+/xi- covariance term, according to       *
 * ESH09, eq. (22).			*
 * ============================================================ */
#define NTH_M_ESH09 12
void scale_mixed_ESH09(const cosmo_lens *model, gsl_matrix *cov, const datcov *dc, error **err)
{
   const double theta[NTH_M_ESH09] = {1.0, 4.9, 10.3, 15.7, 33.0, 45.4, 69.3, 85.7, 106.6, 131.0, 162.0, 180.0};

   /* ESH09, http://www.astro.uni-bonn.de/~teifler/fit-parameters.pdf */
   const double alpha_pp[NTH_M_ESH09][NTH_M_ESH09] = {
     {1.1892,1.3888,1.4163,1.4212,1.4485,1.4726,1.5186,1.5485,1.5841,1.6262,1.6768,1.7063},
     {1.3888,1.2955,1.4126,1.4214,1.4483,1.4723,1.5183,1.5483,1.5838,1.6260,1.6767,1.7061},
     {1.4163,1.4126,1.3243,1.4190,1.4476,1.4712,1.5173,1.5474,1.5830,1.6253,1.6761,1.7055},
     {1.4212,1.4214,1.4190,1.3395,1.4465,1.4695,1.5156,1.5458,1.5816,1.6240,1.6750,1.7045},
     {1.4493,1.4483,1.4476,1.4465,1.3677,1.4607,1.5055,1.5365,1.5733,1.6166,1.6684,1.6982},
     {1.4726,1.4723,1.4712,1.4695,1.4607,1.3802,1.4944,1.5259,1.5637,1.6082,1.6609,1.6911},
     {1.5186,1.5183,1.5173,1.5156,1.5055,1.4944,1.396, 1.4976,1.5370,1.5845,1.6399,1.6712},
     {1.5485,1.5483,1.5474,1.5458,1.5365,1.5259,1.4976,1.4023,1.5135,1.5629,1.6208,1.6532},
     {1.5841,1.5838,1.5830,1.5816,1.5733,1.5637,1.5370,1.5135,1.4078,1.5308,1.5922,1.6264},
     {1.6262,1.6260,1.6253,1.6240,1.6166,1.6082,1.5845,1.5629,1.5308,1.4120,1.5496,1.5864},
     {1.6768,1.6767,1.6761,1.6750,1.6684,1.6609,1.6399,1.6208,1.5922,1.5496,1.4148,1.5279},
     {1.7063,1.7061,1.7055,1.7045,1.6982,1.6911,1.6712,1.6532,1.6264,1.5864,1.5279,1.4157},
   };
   const double beta_pp[NTH_M_ESH09][NTH_M_ESH09] = {
     {2.689, 2.4905,2.1437,2.0102,1.9220,1.9198,1.9322,1.9418,1.9522,1.9627, 1.9734,1.9781},
     {2.4905,2.5137,2.1993,2.0283,1.9232,1.9201,1.9322,1.9417,1.9520,1.9628, 1.9731,1.9782},
     {2.1437,2.1993,2.3725,2.1168,1.9283,1.9210,1.9319,1.9414,1.9518,1.9625, 1.9729,1.9781},
     {2.0102,2.0283,2.1168,2.3110,1.9401,1.9230,1.9314,1.9408,1.9513,1.9621, 1.9728,1.9778},
     {1.9222,1.9232,1.9283,1.9401,2.2460,1.9614,1.9301,1.9377,1.9484,1.9598, 1.9710,1.9765},
     {1.9198,1.9201,1.9210,1.9230,1.9614,2.2323,1.9340,1.9351,1.9451,1.9571, 1.9691,1.9749},
     {1.9322,1.9322,1.9319,1.9314,1.9301,1.9334,2.2234,1.9491,1.9381,1.9496, 1.9634,1.9701},
     {1.9418,1.9417,1.9414,1.9408,1.9377,1.9351,1.9491,2.2219,1.9433,1.943,  1.9580,1.9655},
     {1.9522,1.9520,1.9518,1.9513,1.9484,1.9451,1.9381,1.9433,2.2218,1.94082,1.9497,1.9583},
     {1.9627,1.9628,1.9625,1.9621,1.9598,1.9571,1.9496,1.9433,1.9408,2.2226, 1.9409,1.9473},
     {1.9734,1.9731,1.9729,1.9728,1.9710,1.9691,1.9634,1.9580,1.9497,1.9409, 2.2242,1.9495},
     {1.9781,1.9782,1.9781,1.9778,1.9765,1.9749,1.9701,1.9655,1.9583,1.9473, 1.9495,2.225},
   };

   const double alpha_mm[NTH_M_ESH09][NTH_M_ESH09] = {
     {0.8146, 1.1532, 1.331,  1.4036,0.6502,0.7446,0.6435,1.4427,0.8483,1.115, 0.900, 0.8539},
     {1.1532, 0.9987, 1.3058, 1.3935,1.3981,1.3473,0.6796,1.1596,1.1655,1.1080,1.1260,1.063},
     {1.3307, 1.3058, 1.0946, 1.3548,1.4015,1.375, 1.3997,1.4481,1.3982,0.6967,1.2972,1.1362},
     {1.4036, 1.3935, 1.3548, 1.1431,1.4052,1.3816,1.3711,1.3803,1.3884,1.3681,1.4430,1.4899},
     {0.6519, 1.3981, 1.4015, 1.4052,1.2103,1.3983,1.3711,1.3728,1.3832,1.4006,1.4284,1.4404},
     {0.7447, 1.3473, 1.3750, 1.3816,1.3983,1.2321,1.3761,1.3724,1.3813,1.40,  1.4244,1.4405},
     {0.64346,0.67964,0.13997,1.3711,1.3711,1.3761,1.2562,1.3808,1.3775,1.3931,1.4193,1.4355},
     {1.4427, 1.1596, 1.4481, 1.3803,1.3728,1.3724,1.3808,1.2668,1.3793,1.3878,1.4133,1.4298},
     {0.8482, 1.1655, 1.3982, 1.3884,1.3832,1.3813,1.3775,1.3793,1.2770,1.3833,1.4045,1.4211},
     {1.1150, 1.1080, 0.69672,1.3681,1.4006,1.3995,1.3931,1.3878,1.3833,1.2871,1.3931,1.4081},
     {0.90,   1.1260, 1.2972, 1.4430,1.4284,1.4244,1.4193,1.4133,1.4045,1.3931,1.2973,1.3944},
     {0.8539, 1.063,  1.1362, 1.4899,1.4404,1.4405,1.4354,1.4298,1.4211,1.4081,1.3944,1.3025},
   };
   const double beta_mm[NTH_M_ESH09][NTH_M_ESH09] = {
     {2.3102,2.9769,3.1407,3.0122,1.4087,1.5601,0.7488,4.3083,2.0481,5.0329,1.8159,2.6019},
     {2.9769,2.5973,3.1406,2.9523,2.2487,2.0178,3.0619,2.5425,3.0956,3.9419,3.3497,3.2596},
     {3.1407,3.1406,2.777, 3.0472,2.280, 2.053, 1.8821,1.708, 1.9703,6.7254,3.6954,7.0468},
     {3.0122,2.9523,3.0472,2.8242,2.3463,2.0730,1.9219,1.8869,1.8507,2.1066,1.8902,0.8581},
     {1.4122,2.2487,2.280, 2.3463,2.7801,2.2915,1.9407,1.8973,1.8868,1.898, 1.8862,1.9033},
     {1.5601,2.0178,2.0530,2.0730,2.2915,2.7192,1.9933,1.910, 1.8879,1.8913,1.9095,1.9213},
     {0.7488,3.0619,1.8821,1.9219,1.9407,1.9933,2.6167,2.0248,1.9049,1.8889,1.8982,1.9059},
     {4.3084,2.5425,1.7080,1.8869,1.8973,1.910, 2.0248,2.5626,1.961, 1.8931,1.8964,1.9045},
     {2.0481,3.0956,1.9703,1.8507,1.8868,1.8879,1.9050,1.961, 2.5101,1.9236,1.8944,1.9005},
     {5.0329,3.9419,6.7254,2.1067,1.898, 1.8913,1.8889,1.8930,1.9237,2.4610,1.9060,1.8975},
     {1.8159,3.3497,3.6954,1.8902,1.8862,1.9095,1.8982,1.8964,1.8944,1.9060,2.4166,1.9428},
     {2.6019,3.2596,7.0468,0.8581,1.9033,1.9213,1.906, 1.9045,1.9005,1.8975,1.943, 2.3965},
   };

   const double alpha_pm[NTH_M_ESH09][NTH_M_ESH09] = {
     {0.9198,1.0776, 1.2039, 1.2690,1.3444,1.3588,1.3680,1.3720,1.3776,1.3857,1.3968,1.4036},
     {0.7385,1.8641, 1.2002, 1.2616,1.3433,1.3585, 1.368,1.3721,1.3776,1.3857,1.3968,1.4036},
     {1.1172,1.0939, 0.3245, 1.2212,1.3388,1.3572,1.3679,1.3729,1.3777,1.3857,1.3968,1.4036},
     {1.0298,1.2851, 1.2017, 0.7253,1.3241,1.3544,1.3678,1.3721,1.3777,1.3857,1.3968,1.4036},
     {0.7624,1.5316, 1.4728, 1.4719,0.9712,1.3427,1.3660,1.3730,1.3785,1.3861,1.3968,1.4035},
     {0.7845,1.0214, 1.3619, 1.4555,1.4584,1.0303,1.3650,1.3824,1.3810,1.3871,1.3970,1.4034},
     //{0.5863,0.9360,-5.4805, 1.3951,1.3855,1.4222,1.0895,1.3595,1.3635,1.4203,1.4007,1.4049},
     {0.5863,0.9360, 1.128,  1.3951,1.3855,1.4222,1.0895,1.3595,1.3635,1.4203,1.4007,1.4049},
     {1.1916,1.0523, 0.8196, 1.4040,1.3075,1.3291,1.4146,1.1132,1.3654,1.3648,1.4387,1.4129},
     {0.8716,1.0998, 1.4396, 1.3243,1.2660,1.2710,1.3063,1.3641,1.1337,1.3707,1.3699,1.3430},
     {1.1017,1.1333, 1.1835, 1.1471,1.2451,1.2547,1.2624,1.2748,1.3179,1.1518,1.3768,1.3791},
     {0.8894,1.1091, 1.2446, 1.4713,1.2692,1.2635,1.2663,1.2652,1.2674,1.2897,1.1679,1.3794},
     {1.0537,1.0673, 1.3307, 1.2411,1.2369,1.2882,1.2774,1.2751,1.2713,1.2731,1.3207,1.1755},
   };
   const double beta_pm[NTH_M_ESH09][NTH_M_ESH09] = {
     {2.1759,2.7594,2.9587, 2.9510,2.6699,2.4892,2.2673,2.1796,2.1105,2.0591,2.0232,2.0098},
     {2.1639,4.6092,2.8784, 2.9351,2.6728,2.4916,2.2684,2.1802,2.1109,2.0594,2.0232,2.0099},
     {3.3921,3.2594,1.5224, 3.0572,2.6790,2.4998,2.2725,2.1827,2.1123,2.0601,2.0236,2.0102},
     {3.9545,3.7686,3.5351, 2.2997,2.6384,2.5113,2.2797,2.1871,2.1147,2.0614,2.0243,2.0106},
     {1.5130,3.5641,3.7076, 3.7696,2.8402,2.7344,2.321, 2.2166,2.1312,2.0699,2.0284,2.0135},
     {1.5519,2.2561,3.0393, 3.1508,3.6129,2.9527,2.4350,2.2705,2.1564,2.0818,2.0340,2.0173},
     //{0.6329,2.6217,5.1207, 2.3331,2.4026,2.6896,3.0243,2.4258,2.1906,2.1918,2.0579,2.0318},
     {0.6329,2.6217,2.858,  2.3331,2.4026,2.6896,3.0243,2.4258,2.1906,2.1918,2.0579,2.0318},
     {3.7856,2.9128,3.4400, 1.9722,2.0450,2.1320,2.7362,3.0278,2.3185,2.1055,2.1568,2.060},
     {2.0595,3.5045,2.8142, 0.8366,1.8996,1.9017,2.040, 2.3663,3.0103,2.2253,2.0447,1.9414},
     {4.6080,3.7127,3.7228, 3.1046,1.9377,1.8477,1.8426,1.8956,2.0997,2.9732,2.1490,2.0619},
     {1.8539,3.4068,4.2135, 1.9034,1.6319,1.9126,1.8088,1.8213,1.8395,1.9450,2.9192,2.2337},
     {2.5585,3.1909,5.2090, 5.6088,1.7140,1.9389,1.8113,1.8261,1.8255,1.8484,2.1093,2.8872},
   };

   int nz, mz, i, j, ii, jj, ci, cj;
   double f, alpha, beta, old;
   int pm1, pm2;

   ii = jj = -1;
   for (nz=ci=0; nz<dc->Nzcorr; nz++) {
      for (i=0; i<dc->Ntheta; i++,ci++) {

	 ii = find_bin(dc->theta[i]/arcmin, theta, NTH_M_ESH09, ii, err);
	 forwardError(*err, __LINE__,);

	 pm1 = get_pm(dc->type, i, dc->Ntheta, err);
	 forwardError(*err, __LINE__,);

	 for (mz=cj=0; mz<dc->Nzcorr; mz++) {
	    for (j=0;j<dc->Ntheta;j++,cj++) {

	       jj = find_bin(dc->theta[j]/arcmin, theta, NTH_M_ESH09, jj, err);
	       forwardError(*err, __LINE__,);

	       pm2 = get_pm(dc->type, j, dc->Ntheta, err);
	       forwardError(*err, __LINE__,);

	       switch (2*pm1 - pm2) {

		  case  1 : /* +2 -1: ++ */
		     alpha = alpha_pp[ii][jj];
		     beta  = beta_pp[ii][jj];
		     break;
		  case  3 : /* +2 +1: +- */
		     alpha = alpha_pm[ii][jj];
		     beta  = beta_pm[ii][jj];
		     break;
		  case -3 : /* -2 -1: -+ */
		     alpha = alpha_pm[jj][ii];
		     beta  = beta_pm[jj][ii];
		     break;
		  case -1 : /* -2 +1: -- */
		     alpha = alpha_mm[ii][jj];
		     beta  = beta_mm[ii][jj];
		     break;
		  default :
		     *err = addErrorVA(lensing_pm, "Invalid pm (%d, %d)", *err, __LINE__, pm1, pm2);
		     return;
	       }

	       f  = pow(model->cosmo->Omega_m / 0.25, alpha) * pow(model->cosmo->sigma_8 / 0.9, beta);
	       old = gsl_matrix_get(cov, ci, cj);
	       gsl_matrix_set(cov, ci, cj, f * old);

	    }
	 }
    
      }
   }
}
#undef NTH_M_ESH09

/* ============================================================ *
 * Lensing signal, shear second-order statistics for angular    *
 * scale theta [rad].						*
 * ============================================================ */

double lensing_signal(cosmo_lens *model, double theta, int i_bin, int j_bin, lensdata_t type,
		      decomp_eb_filter_t decomp_eb_filter, const cosebi_info_t *cosebi_info, error **err)
{
   double res, resp, resm, eta;
   const double *a;
   int N;

   testErrorRetVA(type!=decomp_eb && decomp_eb_filter!=decomp_eb_none, lensing_type,
		  "lensdata type (%d) and decomp_eb_filter type (%d) not compatible",
		  *err, __LINE__, 0.0, type, decomp_eb_filter);

   switch (type) {

      case xip :
         res = xi(model, +1, theta, i_bin, j_bin, err);
         forwardError(*err, __LINE__, 0);
         break;

      case xim :
         res = xi(model, -1, theta, i_bin, j_bin, err);
         forwardError(*err, __LINE__, 0);
         break;

      case map2poly :
         res = map2_poly(model, theta, i_bin, j_bin, err);
         forwardError(*err, __LINE__, 0);
         break;

      case map2gauss :
         res = map2_gauss(model, theta, i_bin, j_bin, err);
         forwardError(*err, __LINE__, 0);
         break;

      case gsqr :
         res = gamma2(model, theta, i_bin, j_bin, err);
         forwardError(*err, __LINE__, 0);
         break;

      case decomp_eb :

         if (decomp_eb_filter == COSEBIs_log) {

            int n = (int)round(theta/arcmin);
            testErrorRetVA(n > cosebi_info->n_max, lensing_cosebi_n_max,
                  "COSEBIs mode %d larger than maximum mode %d", *err, __LINE__, 0.0, n, cosebi_info->n_max);
            res = E_cosebi(model, n, cosebi_info->th_min, cosebi_info->th_max,
                  i_bin, j_bin, cosebi_info->path, NULL, err);
            forwardError(*err, __LINE__, 0);

         } else {

            switch (decomp_eb_filter) {
               case FK10_SN        : a = a_FK10_SN;        eta = eta_FK10_SN;        N = N_FK10; break;
               case FK10_FoM_eta10 : a = a_FK10_FoM_eta10; eta = eta_FK10_FoM_eta10; N = N_FK10; break;
               case FK10_FoM_eta50 : a = a_FK10_FoM_eta50; eta = eta_FK10_FoM_eta10; N = N_FK10; break;
               case decomp_eb_none : *err = addErrorVA(lensing_type, "decomp_eb_filter type cannot be 'decomp_eb_none'",
                                           *err, __LINE__, decomp_eb_filter);
                                     return 0.0;
               default             : *err = addErrorVA(lensing_type, "Unknown decomp_eb_filter type %d",
                                           *err, __LINE__, decomp_eb_filter);
                                     return 0.0;
            }
            resp = RR(model, theta*eta, theta, a, N, cheby2, +1, err);
            forwardError(*err, __LINE__, 0.0);
            resm = RR(model, theta*eta, theta, a, N, cheby2, -1, err);
            forwardError(*err, __LINE__, 0.0);
            res = 0.5*(resp + resm); /* E-mode */
         }
         break;

      case pkappa :
         /* Interpret theta as ell, undo arcminute transformation. */
	 res = Pshear(model, theta/arcmin, i_bin, j_bin, err);
	 forwardError(*err, __LINE__, 0.0);
	 break;

      case xipm : case nofz :
      default :
	 *err = addErrorVA(ce_unknown, "Unknown or invalid lensdata type %d(%s)",
			   *err, __LINE__, type, slensdata_t(type));
	 return 0.0;

   }

   return res;
}

/* ============================================================ *
 * Lensing log-likelihood function, for second-order real-space *
 * functions.							*
 * Returns chi^2 = -2 log L					*
 * ============================================================ */
#define NPERBIN 20
double chi2_lensing(cosmo_lens* csm, datcov* dc, int return_model, double **model_array, int *Nmodel,
		    const cosebi_info_t *cosebi_info, error **err)
{
   double *data_minus_model, model, th, dth, w, wtot;
   int i, j, in,i_bin, j_bin, Nzbin;
   double res, logL, lndetC;
   lensdata_t type;
   gsl_matrix_view A, B;
   gsl_matrix *tmp, *tmp2;
   gsl_vector_view x;


   /* Unphysically high baryon fraction, caused infinite chi^2 values. */
   testErrorRet(csm->cosmo->Omega_b/csm->cosmo->Omega_m>BARYON_FRAC, lensing_baryon_fraction,
                "Baryon fraction unphysically high", *err, __LINE__, 0.0);

   Nzbin = csm->redshift->Nzbin;
   testErrorRetVA(Nzbin!=dc->Nzbin, redshift_Nzbin,
		  "Number of redshift bins for model (%d) inconsistent with data (%d)",
		  *err, __LINE__, 0, Nzbin, dc->Nzbin);

   *model_array = malloc_err(sizeof(double)*dc->n, err);
   forwardError(*err, __LINE__, 0);

   data_minus_model = malloc_err(sizeof(double)*dc->n, err);
   forwardError(*err, __LINE__, 0);


   /* Fill model vector */

   for (i_bin=0,in=0; i_bin<Nzbin; i_bin++) {
      for (j_bin=i_bin; j_bin<Nzbin; j_bin++) {
         for (j=0; j<dc->Ntheta; j++,in++) {

            /* First half of data vector = xi+, second half = xi- */
            if (dc->type==xipm) {
               if (j<dc->Ntheta/2) type = xip;
               else type = xim;
            } else {
               type = dc->type;
            }

            if (dc->format==angle_mean) {

               /* Average model over bin width */
               dth = (dc->theta2[j]-dc->theta[j])/(double)(NPERBIN-1.0);
               for (i=0,model=0.0,th=dc->theta[j]; i<NPERBIN; i++,th+=dth) {
                  model += lensing_signal(csm, th, i_bin, j_bin, type, dc->decomp_eb_filter, NULL, err);
                  forwardError(*err, __LINE__, 0.0);
               }
               model /= (double)NPERBIN;

            }  else if (dc->format==angle_wlinear) {

               /* Linear weighting over angular bin */
               dth = (dc->theta2[j]-dc->theta[j])/(double)(NPERBIN-1.0);
               for (i=0,model=wtot=0.0,th=dc->theta[j]; i<NPERBIN; i++,th+=dth) {
                  w      = th/arcmin;
                  model += w*lensing_signal(csm, th, i_bin, j_bin, type, dc->decomp_eb_filter, NULL, err);
                  forwardError(*err, __LINE__, 0.0);
                  wtot  += w;
               }
               model /= wtot;

            }  else if (dc->format==angle_wquadr) {

               /* Quadratic weighting over angular bin */
               dth = (dc->theta2[j]-dc->theta[j])/(double)(NPERBIN-1.0);
               for (i=0,model=wtot=0.0,th=dc->theta[j]; i<NPERBIN; i++,th+=dth) {
                  w      = dc->a1*th/arcmin + dc->a2*dsqr(th/arcmin);
                  model += w*lensing_signal(csm, th, i_bin, j_bin, type, dc->decomp_eb_filter, NULL, err);
                  forwardError(*err, __LINE__, 0.0);
                  wtot  += w;
               }
               model /= wtot;

            } else {

               /* Model at bin center */
               model = lensing_signal(csm, dc->theta[j], i_bin, j_bin, type, dc->decomp_eb_filter, cosebi_info, err);
               forwardError(*err, __LINE__, -1);

            }

            testErrorRetVA(in>=dc->n, math_overflow, "Overflow, data index %d>=%d", *err, __LINE__, 0, in, dc->n);

            (*model_array)[in]   = model;
            data_minus_model[in] = dc->data[in] - model;

         }
      }
   }

   *Nmodel = in;

   /* Calculate the log-likelihood */

   res = 0.0;

   if (dc->usecov) {

      gsl_set_error_handler_off();
      x = gsl_vector_view_array(data_minus_model, dc->n);

      switch (dc->cov_scaling) {
         case cov_const :

            A = gsl_matrix_view_array(dc->cov[0], dc->n, dc->n);
            lndetC = dc->lndetC;
            break;

         case cov_ESH09 :

            /* Copy shot noise term D -> tmp */
            tmp = gsl_matrix_alloc(dc->n, dc->n);
            A   = gsl_matrix_view_array(dc->cov[2], dc->n, dc->n);
            i   = gsl_matrix_memcpy(tmp, &A.matrix);
            testErrorRetVA(i != 0, math_unknown, "Matrix copying failed, gsl return value = %d",
                  *err, __LINE__, -1, i);


            /* Copy mixed term M -> tmp2 */
            tmp2 = gsl_matrix_alloc(dc->n, dc->n);
            B    = gsl_matrix_view_array(dc->cov[1], dc->n, dc->n);
            i    = gsl_matrix_memcpy(tmp2, &B.matrix);
            testErrorRetVA(i != 0, math_unknown, "Matrix copying failed, gsl return value = %d",
                  *err, __LINE__, -1, i);

            /* Scale mixed term */
            scale_mixed_ESH09(csm, tmp2, dc, err);
            //write_matrix(tmp2->data, dc->n, "M_ESH09", err); forwardError(*err, __LINE__, -1.0);

            /* Add mixed to shot noise -> tmp */
            A = gsl_matrix_submatrix(tmp, 0, 0, dc->n, dc->n);
            B = gsl_matrix_submatrix(tmp2, 0, 0, dc->n, dc->n);
            i = gsl_matrix_add(&A.matrix, &B.matrix);
            testErrorRetVA(i != 0, math_unknown, "Matrix addition failed, gsl return value = %d",
                  *err, __LINE__, -1, i);


            /* Copy cosmic variance term V -> tmp2 */
            B = gsl_matrix_view_array(dc->cov[0], dc->n, dc->n);
            i = gsl_matrix_memcpy(tmp2, &B.matrix);
            testErrorRetVA(i != 0, math_unknown, "Matrix copying failed, gsl return value = %d",
                  *err, __LINE__, -1, i);

            /* Scale cosmic variance term */
            scale_cosmic_variance_ESH09(csm, tmp2, dc, err);

            /* Add cosmic variance to shot+mixed -> tmp */
            A = gsl_matrix_submatrix(tmp, 0, 0, dc->n, dc->n);
            B = gsl_matrix_submatrix(tmp2, 0, 0, dc->n, dc->n);
            i = gsl_matrix_add(&A.matrix, &B.matrix);
            testErrorRetVA(i != 0, math_unknown, "Matrix addition failed, gsl return value = %d",
                  *err, __LINE__, -1, i);
            //write_matrix(tmp->data, dc->n, "tot", err); forwardError(*err, __LINE__, -1.0);

            gsl_matrix_free(tmp2);

            /* Replace res->cov with L, where C = L L^T */
            testErrorRet(gsl_linalg_cholesky_decomp(&A.matrix) == GSL_EDOM, mv_cholesky,
                  "Cholesky decomposition of covariance matrix failed", *err, __LINE__, -1);

            /* Determinant of the covariance.  ln C = 2 ln L */
            //write_matrix(A.matrix.data, dc->n, "L", err); forwardError(*err, __LINE__, -1.0);
            lndetC = 2.0 * ln_determinant(A.matrix.data, dc->n, err);
            forwardError(*err, __LINE__, -1.0);

            break;

         default :

            *err = addErrorVA(ce_unknown, "Unknown or invalid cov_scaling type %d(%s)",
                  *err, __LINE__, dc->cov_scaling, scov_scaling_t(dc->cov_scaling));
            return -1.0;

      }


      /* A has to contain L, with C = L^T L    *
       * Calculates x = L^{-1} . (data-model). *
       * A remains unchanged.		       */
      gsl_blas_dtrsv(CblasLower, CblasNoTrans, CblasNonUnit, (const gsl_matrix*)(&A.matrix), &x.vector);

      /* x^T . x */
      for (i=0; i<dc->n; i++) {
         res += dsqr(gsl_vector_get(&x.vector, i));
      }

      if (dc->cov_scaling == cov_ESH09) {
         gsl_matrix_free(tmp);
      }

   } else {

      *err = addError(lensing_inconsistent, "usecov = 0 not valid", *err, __LINE__);
      return -1.0;

   }

   testErrorRetVA(res<0.0, math_negative, "Negative chi^2 %g. Maybe the covariance matrix is not positive",
         *err, __LINE__, -1.0, res);

   if (! return_model) free(*model_array);
   logL = -0.5 * (res + dc->n * ln2pi + lndetC);

   /* New v1.2: Problem with infinite weights solved */

   return logL;
}
#undef NPERBIN

/* Solve the equation Nperm = N(N+1)(N+2)/6 for N and returns N */
#define NMAX 50
int Nperm_to_Ntheta(int Nperm, error **err)
{
   int N;

   for (N=1; N<NMAX; N++) {
      if (N * (N + 1) * (N + 2) / 6 == Nperm) return N;
   }

   *err = addErrorVA(lensing_nperm, "Number of angular bins not found, maybe > %d?",
		     *err, __LINE__, NMAX);
   return 0;
}
#undef NMAX

/* ============================================================ *
 * Reads data file into structure dc (has to be initialised     *
 * before). The file has to have the structure                  *
 *   theta_1 theta_2 theta_3 [arcmin] <M_ap^3>.                 *
 * If Nzbin>0 on input performs a consistency check.		*
 * ============================================================ */
void read_data_3rd(datcov *dc, char data_name[], int Nzbin, error **err)
{
   int n, i, j, k, Ncol;
   size_t Nrec, Nperm, Nzcorr;
   double dNcol;
   double *ptr;

   /* Read file to double pointer */
   Nrec = 0;
   ptr = (double*)read_any_list_count(data_name, &Nrec, "%lg", sizeof(double), &Nperm, err);
   forwardError(*err, __LINE__,);

   /* ============================================================ *
    * Nperm is the number of permutations without repetition.      *
    * ============================================================ */

   /* Number of lines and columns */
   dNcol = (double)Nrec/(double)Nperm;
   Ncol = Nrec/Nperm;
   testErrorRet(fabs(dNcol-(double)Ncol)>EPSILON, lensing_inconsistent, "Columns do not have equal length",
		*err, __LINE__,);

   Nzcorr = Ncol - 3;

   /* Fill datcov structure */
   dc->Ntheta  = Nperm_to_Ntheta(Nperm, err);                  forwardError(*err, __LINE__,);
   dc->Ntheta2 = 0;
   dc->Nzbin   = Nperm_to_Ntheta(Nzcorr, err);                 forwardError(*err, __LINE__,);

   dc->Nzcorr  = Nzcorr;
   dc->n       = Nperm*Nzcorr;

   dc->data  = malloc_err(sizeof(double)*dc->n, err);          forwardError(*err, __LINE__,);
   dc->theta = malloc_err(sizeof(double)*dc->Ntheta, err);     forwardError(*err, __LINE__,);
   dc->theta2 = NULL;  /* Only valid format is angle_center */
   dc->var    = NULL;

   /* First three columns: angular scales (theta_1 theta_2 theta_3) */
   for (i=n=0; i<dc->Ntheta; i++) {
      for (j=i; j<dc->Ntheta; j++) {
         for (k=j; k<dc->Ntheta; k++,n++) {

            /* Store angular scale from theta_k for i,j=0 */
            if (i==0 && j==0) {
               dc->theta[k] = ptr[Ncol*k+2]*arcmin;
            }

            testErrorRetVA(n>=dc->n, math_overflow, "Index overflow (%d>=%d)",
                  *err, __LINE__,, n, dc->n);

            dc->data[n] = ptr[n*Ncol+3];

         }
      }
   }

}

/* ============================================================ *
 * Reads two files with 2nd- and 3rd-order data, respectively.  *
 * Copies the data to res. The 3rd-order file can be diagonal   *
 * or general (3theta) format					                      *
 * ============================================================ */
void read_data_2nd_3rd(datcov *res, char *dataname, char *dataname2, error **err)
{
   datcov *dc, *dc2;

   dc  = malloc_err(sizeof(datcov), err);              forwardError(*err,__LINE__,);
   dc2 = malloc_err(sizeof(datcov), err);              forwardError(*err,__LINE__,);

   /* Read second-order file */
   dc->format = res->format;
   read_data_tomo(dc, dataname, 0, second_order, err); forwardError(*err, __LINE__,);

   /* Read third-order file, verifying same number of redshift bins as 2nd-order */
   if (res->type==map2gauss_map3gauss_diag || res->type==decomp_eb_map3gauss_diag) {
      read_data_tomo(dc2, dataname2, dc->Nzbin, third_order, err);  forwardError(*err, __LINE__,);
   } else {
      read_data_3rd(dc2, dataname2, dc->Nzbin, err);   forwardError(*err, __LINE__,);
   }

   res->Ntheta  = dc->Ntheta;
   res->Ntheta2 = dc2->Ntheta;
   res->n       = dc->n + dc2->n;
   res->Nzbin   = dc->Nzbin;
   res->Nzcorr  = dc->Nzcorr;

   testErrorRetVA(dc->Nzbin != dc2->Nzbin, lensing_nzbin,
	       "Different number of redshift bins for 2nd (%d) and 3rd (%d) order",
	       *err, __LINE__,, dc->Nzbin, dc2->Nzbin);

   res->theta   = malloc_err(sizeof(double)*res->Ntheta, err);   forwardError(*err, __LINE__,);
   res->theta2  = malloc_err(sizeof(double)*res->Ntheta2, err);  forwardError(*err, __LINE__,);
   memcpy(res->theta, dc->theta, res->Ntheta*sizeof(double));
   memcpy(res->theta2, dc2->theta, res->Ntheta2*sizeof(double));

   res->data    = malloc_err(sizeof(double)*res->n, err);        forwardError(*err, __LINE__,);
   memcpy(res->data, dc->data, dc->n*sizeof(double));
   memcpy(res->data + dc->n, dc2->data, dc2->n*sizeof(double));

   del_data_cov(&dc);
   del_data_cov(&dc2);
}
