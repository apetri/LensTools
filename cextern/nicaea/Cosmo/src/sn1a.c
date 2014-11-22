/* ============================================================ *
 * sn1a.c							*
 * Pierre Astier, Martin Kilbinger 2006-2009 			*
 * ============================================================ */


#include "sn1a.h"

cosmo_SN *init_parameters_SN(double OMEGAM, double OMEGADE, double W0_DE, double W1_DE,
			     double *W_POLY_DE, int N_POLY_DE,
			     double H100, double OMEGAB, double OMEGANUMASS, 
			     double NEFFNUMASS, double NORM, double NSPEC,
			     nonlinear_t NONLINEAR, transfer_t TRANSFER, growth_t GROWTH,
			     de_param_t DEPARAM, norm_t normmode,
			     chi2mode_t CHI2MODE, double THETA1[], double THETA2[], double BETA_D,
			     double AMIN, error **err)
{
   cosmo_SN *res;
   int i;

   res = (cosmo_SN*)malloc_err(sizeof(cosmo_SN), err);
   forwardError(*err, __LINE__, NULL);

   res->cosmo = init_parameters(OMEGAM, OMEGADE, W0_DE, W1_DE, W_POLY_DE, N_POLY_DE,
				H100, OMEGAB, OMEGANUMASS, NEFFNUMASS, NORM, NSPEC,
				NONLINEAR, TRANSFER, GROWTH, DEPARAM, normmode, AMIN, err);
   forwardError(*err, __LINE__, 0);

   for (i=0; i<NDER; i++) {
      res->Theta1[i] = THETA1[i];
   }
   for (i=0; i<NLCP; i++) {
      res->Theta2[i]       = THETA2[i];
      res->Theta2_denom[i] = THETA2[i];
   }
   res->beta_d = BETA_D;

   res->chi2mode = CHI2MODE;

   return res;
}

cosmo_SN *copy_parameters_SN_only(cosmo_SN *source, error **err)
{
   cosmo_SN* res;
   res = init_parameters_SN(source->cosmo->Omega_m, source->cosmo->Omega_de, source->cosmo->w0_de,
			    source->cosmo->w1_de,
			    source->cosmo->w_poly_de, source->cosmo->N_poly_de,
			    source->cosmo->h_100, source->cosmo->Omega_b,
			    source->cosmo->Omega_nu_mass, source->cosmo->Neff_nu_mass, 
			    source->cosmo->normalization, source->cosmo->n_spec, 
			    source->cosmo->nonlinear, source->cosmo->transfer, source->cosmo->growth,
			    source->cosmo->de_param, source->cosmo->normmode, source->chi2mode,
			    source->Theta1, source->Theta2, source->beta_d, source->cosmo->a_min, err);
   forwardError(*err, __LINE__, NULL);
   return res;
}

void read_cosmological_parameters_SN(cosmo_SN **self, FILE *F, error **err)
{
   cosmo_SN *tmp;
   struct { char cosmo_file[128], nofz_file[128]; } tmp2;
   config_element c = {0, 0.0, ""};
   FILE *FD;
   int i;
   char stmp[128];

   tmp = set_cosmological_parameters_to_default_SN(err);
   forwardError(*err, __LINE__,);

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

   /* SNIa parameters */
   //CONFIG_READ_ARR(tmp, Theta1, d, i, NDER, s, F, c, err);  /* Photometric calibrations */
   CONFIG_READ_ARR(tmp, Theta2, d, i, NLCP, stmp, F, c, err);  /* -M, alpha, -beta         */

   tmp->beta_d = tmp->stretch = tmp->color = 0.0;

   *self = copy_parameters_SN_only(tmp, err);
   forwardError(*err, __LINE__,);

   free_parameters_SN(&tmp);
}

void updateFrom_SN(cosmo_SN* avant, cosmo_SN* apres, error **err)
{
   updateFrom(avant->cosmo, apres->cosmo, err);
   forwardError(*err, __LINE__,);
}

cosmo_SN* set_cosmological_parameters_to_default_SN(error **err)
{
   int i;
   double theta1[NDER], theta2[NLCP], h100, amin, beta_d;

   for (i=0; i<NDER; i++) theta1[i] = 0.0;

   /* Astier 2006 */
   h100 = 0.7;
   theta2[0] = 19.31 - 5.0*log10(h100/0.7);    /* - M */
   theta2[1] = 1.52;                           /* alpha */
   theta2[2] = -1.57;                          /* -beta */
   theta2[3] = 0.0;                            /* -beta_z */
   beta_d    = 0.0;

   amin = 1.0/(1.0+2.0);  /* zmax of 2 is enough for SNLS */

   return init_parameters_SN(0.26,0.74,-1.0,0.0,NULL,0,h100,0.04,0.0,0.0,0.85,
			     1.0,smith03,eisenhu,growth_de,linder,norm_s8,
			     chi2_simple,theta1,theta2,beta_d,amin,err);
}

cosmo_SN* set_cosmological_parameters_to_best_fit_SNLS_WMAP5(error **err)
{
   int i;
   double theta1[NDER], theta2[NLCP], h100, amin, beta_d;

   for (i=0; i<NDER; i++) theta1[i] = 0.0;

   /* Kilbinger et al. 2009 */
   h100 = 0.719;
   theta2[0] = 19.3137 - 5.0*log10(h100/0.7);  /* - M */
   theta2[1] = 1.62;                           /* alpha */
   theta2[2] = -1.80;                          /* -beta */
   theta2[3] = 0.0;                            /* -beta_z */
   beta_d    = 0.0;

   amin = 1.0/(1.0+2.0);  /* zmax of 2 is enough for SNLS */

   return init_parameters_SN(0.257,0.743,-1.025,0.0,NULL,0,h100,0.0433,0.0,0.0,0.807,
			     0.962,smith03,eisenhu,growth_de,linder,norm_s8,
			     chi2_simple,theta1,theta2,beta_d,amin,err);
}

cosmo_SN* set_cosmological_parameters_to_best_fit_SNLS(error **err)
{
   int i;
   double theta1[NDER], theta2[NLCP], h100, amin, beta_d;

   for (i=0; i<NDER; i++) theta1[i] = 0.0;

   h100 = 0.719;
   theta2[0] = 19.306 - 5.0*log10(h100/0.7);    /* - M */
   theta2[1] = 1.515;                           /* alpha */
   theta2[2] = -1.56;                           /* -beta */
   theta2[3] = 0.0;                             /* -beta_z */
   beta_d    = 0.0;

   amin = 1.0/(1.0+2.0);  /* zmax of 2 is enough for SNLS */

   return init_parameters_SN(0.266,0.734,-1.0,0.0,NULL,0,h100,0.0433,0.0,0.0,0.807,
			     0.962,smith03,eisenhu,growth_de,linder,norm_s8,
			     chi2_simple,theta1,theta2,beta_d,amin,err);
}

cosmo_SN* set_cosmological_parameters_to_best_fit_Union(error **err)
{
   int i;
   double theta1[NDER], theta2[NLCP], h100, amin, beta_d;

   for (i=0; i<NDER; i++) theta1[i] = 0.0;

   h100 = 0.719;
   theta2[0] = 19.31 - 5.0*log10(h100/0.7);    /* - M */
   theta2[1] = 1.37;                           /* alpha */
   theta2[2] = -2.45;                          /* -beta */
   theta2[3] = 0.0;                            /* -beta_z */
   beta_d    = 0.0;

   amin = 1.0/(1.0+2.0);  /* zmax of 2 is enough for Union */

   return init_parameters_SN(0.291,0.709,-1.0,0.0,NULL,0,h100,0.0433,0.0,0.0,0.807,
			     0.962,smith03,eisenhu,growth_de,linder,norm_s8,
			     chi2_simple,theta1,theta2,beta_d,amin,err);
}

cosmo_SN* set_cosmological_parameters_to_EdS_SN(error **err)
{
   int i;
   double theta1[NDER], theta2[NLCP], h100, amin, beta_d;

   for (i=0; i<NDER; i++) theta1[i] = 0.0;

   /* Astier 2006 */
   h100 = 0.7;
   theta2[0] = 19.31 - 5.0*log10(h100/0.7);    /* - M */
   theta2[1] = 1.52;
   theta2[2] = -1.57;
   theta2[3] = 0.0;
   beta_d    = 0.0;

   amin = 1.0/(1.0+2.0);  /* zmax of 2 is enough for SNLS */

   return init_parameters_SN(1.0,0.0,-1.0,0.0,NULL,0,h100,0.04,0.0,0.0,0.85,1.0,
			     smith03,eisenhu,growth_de,linder,norm_s8,
			     chi2_simple,theta1,theta2,beta_d,amin,err);
}

void free_parameters_SN(cosmo_SN **self)
{
   cosmo_SN *s;
   s = *self;
   free_parameters(&(s->cosmo));
   free(s);
   s = NULL;
}

/* ************    I/O stuff */

void dump_param_SN(cosmo_SN *self, FILE *F)
{
   int i;

   if (F==NULL) F = stderr;
   dump_param(self->cosmo, F);
   fprintf(F, "# ");
   for (i=0; i<NDER; i++) fprintf(F, "  Th1%d", i);
   for (i=0; i<NLCP; i++) fprintf(F, "   Th2%d", i);
   fprintf(F, " m\n");
   fprintf(F, "# ");
   for (i=0; i<NDER; i++) fprintf(F, "% .3f", self->Theta1[i]);
   for (i=0; i<NLCP; i++) fprintf(F, "% 7.3f", self->Theta2[i]);
   fprintf(F, " %d\n", self->chi2mode);
}

/* Reads a single line of supernovae data */
void readSnData(char *line, SnData *sndata, sndatformat_t sndatformat, error **err)
{
   char nm[100], *p=line;
   int ok, k;
   int nread=0;

  ok = sscanf(p,"%s %n",nm, &nread);
  testErrorRet(ok!=1, ce_file, "Bad line", *err, __LINE__,);
  /* on s'en fou du nom  */
  p += nread;
 
  ok = read_double(&p,&(sndata->z));
  testErrorRet(ok==0, ce_file, "Bad line", *err, __LINE__, );

  ok = read_double(&p,&(sndata->musb));
  testErrorRet(ok==0, ce_file, "Bad line", *err, __LINE__, );
  ok = read_double(&p,&(sndata->cov[0][0]));
  testErrorRet(ok==0, ce_file, "Bad line", *err, __LINE__, );

  ok = read_double(&p,&(sndata->s));
  testErrorRet(ok==0, ce_file, "Bad line", *err, __LINE__, );
  ok = read_double(&p,&(sndata->cov[1][1]));
  testErrorRet(ok==0, ce_file, "Bad line", *err, __LINE__, );

  ok = read_double (&p,&(sndata->c));
  testErrorRet(ok==0, ce_file, "Bad line", *err, __LINE__, );
  ok = read_double(&p,&(sndata->cov[2][2]));
  testErrorRet(ok==0, ce_file, "Bad line", *err, __LINE__, );

  ok = read_double(&p,&(sndata->cov[0][1]));
  testErrorRet(ok==0, ce_file, "Bad line", *err, __LINE__, );
  ok = read_double(&p,&(sndata->cov[0][2]));
  testErrorRet(ok==0, ce_file, "Bad line", *err, __LINE__, );  
  ok = read_double(&p,&(sndata->cov[1][2]));

  // convert diagonal r.m.s. to variance
  sndata->cov[0][0] *= sndata->cov[0][0];
  sndata->cov[1][1] *= sndata->cov[1][1];
  sndata->cov[2][2] *= sndata->cov[2][2];
  // symetrise
  sndata->cov[1][0] = sndata->cov[0][1];
  sndata->cov[2][0] = sndata->cov[0][2];
  sndata->cov[2][1] = sndata->cov[1][2];

  for (k=0; k<NDER; ++k) {
     if (sndatformat==SNLS_firstyear) {
	ok = read_double(&p,&(sndata->derivative[k]));
	testErrorRet(ok==0, ce_file, "Bad line", *err, __LINE__, );
     } else {
	sndata->derivative[k] = 0;
     }
  }
}

/* Reads a file with supernovae data */
SnSample *SnSample_read(const char *FileName, sndatformat_t sndatformat, error **err)
{
   int i, cur=-1, Nsample, dummy;
   double peculiar_vel, int_disp;
   SnSample *sn;
   char key[256], line[4096], *sval;
   FILE *f;

   //fprintf(stderr, "Reading %s\n", FileName);

   /* default values */
   int_disp = 0.13;
   peculiar_vel = 300.0;

   sn = malloc_err(sizeof(SnSample), err);     forwardError(*err, __LINE__, NULL);
   sn->W1       = NULL;
   sn->W1dim    = 0;
   sn->logdetW1 = 0.0;

   f = fopen(FileName,"r");
   if (f==NULL) {
      sprintf(line, "Could not open file %s", FileName);
      *err = addError(ce_file, line, *err, __LINE__);
      return NULL;
   }

   /* Read twice: 1st time for header and number of lines, 2nd time read data */
   for (i=0,Nsample=0; i<=1; i++) {

      if (i==1) {
	 sn->data = (SnData*)malloc_err(Nsample*sizeof(SnData), err);
	 forwardError(*err, __LINE__, NULL);
	 fseek(f, 0, SEEK_SET);
	 cur = 0;
      }

      while (fgets(line,4096,f))
      {
	 if (line[0] == '#') continue; // comment
	 if (line[strlen(line)-1] == '\n') line[strlen(line)-1] = '\0';
	 if (line[0] == '@') 
	 {
	    if (i==1) continue;

	    sval = key+128;
	    if (sscanf(line+1,"%s %s",key,sval) != 2)
	    {
	       *err = addError(ce_file, "Could not read line", *err, __LINE__);
	       /* Pierre: continue */
	       return NULL;
	    }
	    if (strcmp(key,"INSTRINSIC_DISPERSION")==0)
	      if (!read_double(&sval,&int_disp)) {
		 *err = addError(ce_file, "Could not read line", *err, __LINE__);
		 return NULL;
	      }
	    if (strcmp(key,"PECULIAR_VELOCITY")==0)
	      if (!read_double(&sval,&peculiar_vel)) {
		 *err = addError(ce_file, "Could not read line", *err, __LINE__);
		 return NULL;
	      }
	    if (strcmp(key,"THETA1_ERRORS")==0 && (NDER != 0))
	    {
	       if (i==0) {
		  sn->W1 = readASCII(sval, &(sn->W1dim), &dummy, err);
		  forwardError(*err, __LINE__, NULL);
		  testErrorRet(sn->W1dim!=dummy, ce_nonsquare, "Not a square matrix", *err, 
			       __LINE__, NULL);
		  sm2_inverse(sn->W1, sn->W1dim, err);
		  forwardError(*err, __LINE__, NULL);
	       }
	    }
	 }
	 else  // read data
	 {
	    if (i==0) {
	       Nsample++;
	    } else {
	       readSnData(line, sn->data+cur, sndatformat, err);
	       forwardError(*err, __LINE__, 0);
	       cur++;
	    }
	 }
      }   /* while line */

   }   /* for i */

   sn->int_disp = int_disp;
   sn->sig_mu_pec_vel = 5.0*peculiar_vel/3.0e5/log(10.0);
   sn->Nsample = Nsample;

   for (i=0,sn->logdetW1=0.0; i<sn->W1dim; i++) {
      sn->logdetW1 += log(sn->W1[i*sn->W1dim+i]);
   }

   return sn;
}


/* Prints SN Sample to file */
void out_SnSample(const SnSample *sn, FILE *F)
{
   int i, j, k;
   FILE *OUT;

   if (F==NULL) OUT = stdout;
   else OUT = F;

   fprintf(OUT, "# SnSample\n");
   fprintf(OUT, "# Nsample = %d int_disp = %f sig_mu_pec_vel = %g\n",
	   sn->Nsample, sn->int_disp, sn->sig_mu_pec_vel);

   fprintf(OUT, "# z musb s c dl mu_c derivative[8] covd[%d][%d]\n", NLCP, NLCP);
   for (i=0; i<sn->Nsample; i++) {
      fprintf(OUT, "%f %f %f %f    %7.2f %f", sn->data[i].z, sn->data[i].musb, sn->data[i].s,
	      sn->data[i].c, sn->data[i].dl, sn->data[i].mu_c);
      fprintf(OUT, "   ");
      for (j=0; j<NDER; j++) {
	 fprintf(OUT, " % f", sn->data[i].derivative[j]);
      }
      fprintf(OUT, "   ");
      for (j=0; j<NLCP; j++) {
	 for (k=0; k<NLCP; k++) {
	    fprintf(OUT, " % .3e", sn->data[i].cov[j][k]);
	 }
      }
      fprintf(OUT, "\n");
   }
}

void out_model(const cosmo_SN *cosmo, FILE *F, error **err)
{
   double z, dlum;

   fprintf(F, "# z Dlum [Mpc/h] mu_c\n");

   for (z=0.002; z<=1.5; z+=0.002) {
      dlum = D_lum(cosmo->cosmo, 1.0/(1.0+z), err);
      forwardError(*err, __LINE__,);
      fprintf(F, "%f %10.3f % f\n", z, dlum, distance_module(cosmo->cosmo, dlum, err));
      forwardError(*err, __LINE__,);
   }
}

void SetDl(cosmo_SN *self, SnSample *sn, error **err)
{
   int i;

   for (i=0; i<sn->Nsample; i++) {
      sn->data[i].dl   = D_lum(self->cosmo, 1.0/(sn->data[i].z+1.0), err);
      forwardError(*err, __LINE__,);
      testErrorRet(sn->data[i].dl<=0, ce_negative, "D_lum not positive", *err, __LINE__,);

      /* Distance module is 5 log10(dl/10pc). The dependence on the Hubble constant
	 is put into the (nuissance) parameter M (Theta2[0]) */
      sn->data[i].mu_c = distance_module(self->cosmo, sn->data[i].dl, err);
      forwardError(*err, __LINE__,);
   }
}

/* Returns Theta2^T * cov * Theta2 with Theta2 = (1, alpha, -beta) and cov = cov(m_B^*,s,c) *
 * For chi2mode=chi2_beta_z, Theta2[2] = -beta + beta_z*z.                 		    */
double DistModVariance(const SnData *snd, const double *Theta2)
{
   double res;

   res = snd->cov[0][0]+2.0*Theta2[1]*(snd->cov[0][1]+Theta2[2]*snd->cov[1][2])
     +2.0*Theta2[2]*snd->cov[0][2]+snd->cov[1][1]*Theta2[1]*Theta2[1]
     +snd->cov[2][2]*Theta2[2]*Theta2[2];

   return res;
}

/* Returns the scalar product x^T * A * x */
double vect_scalar_product(const double *x, const double *A, int N)
{
   int i, j;
   double res;

   for (j=0,res=0.0; j<N; j++) {
      for (i=0; i<N; i++) {
	 res += x[i]*A[i+j*N]*x[j];
      }
   }

   return res;
}

double int_for_Nhalo_z(double z, void *intpar, error **err)
{
   cosmo *self;
   double res;

   self = (cosmo*)intpar;

   res  = dsqr(1.0+z)/sqrt(Esqr(self, 1.0/(1.0+z), 0, err));
   forwardError(*err, __LINE__, 0.0);
   /* Correction for rest-frame shift of absorption wave-length */
   res *= pow(1.0+z, 1.2);

   return res;
}

double Nhalo_z(cosmo *self, double z, error **err)
{
   double res;
   double sigma, n;

   /* Menard et al. 2009 */
   n     = 0.037;         /* Mpc^-3 h^-1 */
   sigma = 0.01*pi;       /* Mpc^2       */
   res  = sm2_qromberg(int_for_Nhalo_z, (void*)self, 0.0, z, 1.0e-6, err);
   forwardError(*err, __LINE__, 0.0);

   return res*R_HUBBLE*sigma*n;
}

double chi2_SN_residual(const cosmo_SN *cosmo, const SnSample *sn, error **err)
{
   int i, j, k;
   double chi2, logdetC, res, x[3], tmp;
   mvdens *g;

   g = mvdens_alloc(3, err); forwardError(*err, __LINE__, 0.0);

   for (i=0,chi2=0.0,logdetC=0.0; i<sn->Nsample; i++) {

      /* === Data === */
      g->mean[0] = sn->data[i].musb;
      g->mean[1] = sn->data[i].s;
      g->mean[2] = sn->data[i].c;

      /* === Model === */
      x[0] = -cosmo->Theta2[0]
	- cosmo->Theta2[1]*(cosmo->stretch-1.0)
	- cosmo->Theta2[2]*cosmo->color + sn->data[i].mu_c;
      x[1] = cosmo->stretch;
      x[2] = cosmo->color;


      /* Covariance */
      for (j=0; j<3; j++) {
	 for (k=0; k<3; k++) {
	    g->std[j*3+k] = sn->data[i].cov[j][k];
	 }
      }
      g->std[0] += dsqr(sn->int_disp) + dsqr(sn->sig_mu_pec_vel/sn->data[i].z);


      res = (g->mean[1]-x[1]) * 1 * (g->mean[1]-x[1]);
      goto end;

      tmp = sm2_inverse(g->std, 3, err);  forwardError(*err, __LINE__, 0.0);
      logdetC += log(tmp);

      res = 0.0;
      for (j=0; j<3; j++) {
	 for (k=0; k<3; k++) {
	    res += (g->mean[j]-x[j]) * g->std[j*3+k] * (g->mean[k]-x[k]);
	 }
      }

   end:
      chi2 += res;

   }

   mvdens_free(&g);

   return -0.5*(chi2 + logdetC);

}

/* ============================================================ *
 * Returns the negative log-likelihood.				*
 * Variables: Theta2[0] = -M, Theta2[1] = alpha,		*
 * Theta2[2] = -beta.						*
 * ============================================================ */
double chi2_SN(const cosmo_SN *cosmo, const SnSample *sn,
	       mvdens *data_beta_d, int wTheta1, int add_logdetCov, error **err)
{
   double chi2, res, res_all, wi, mu_meas, sigintsqr, kTheta1, Theta2[NLCP], Theta2_betaz[NLCP],
     logdetC, dmv;
   double c_d;
   int i, j;


   if (cosmo->chi2mode==chi2_residual) {
      res = chi2_SN_residual(cosmo, sn, err);
      forwardError(*err, __LINE__, 0.0);
      return res;
   }

   for (i=0; i<NLCP; i++) Theta2[i] = cosmo->Theta2[i];
   if (cosmo->chi2mode==chi2_Theta2_denom_fixed) {
      /* Copy the fixed denominator values for alpha and beta */
      Theta2[1] = cosmo->Theta2_denom[1];
      Theta2[2] = cosmo->Theta2_denom[2];
   }

   sigintsqr = dsqr(sn->int_disp);

   for (i=0,chi2=0.0,logdetC=0.0; i<sn->Nsample; i++) {

      /* m - M */
      mu_meas = sn->data[i].musb + cosmo->Theta2[0];

      if (cosmo->chi2mode!=chi2_no_sc) {

	 if (cosmo->chi2mode!=chi2_betaz) {

	    /* Add stretch and color terms alpha*(s-1) - beta*c */
	    mu_meas += cosmo->Theta2[1]*(sn->data[i].s-1.0) + cosmo->Theta2[2]*sn->data[i].c;

	    if (cosmo->chi2mode==chi2_dust) {

	       /* Correction for intergalactic dust along line of sight: *
		* add (beta - beta_d)*c_d				 */

	       c_d = sn->data[i].dust;

	       mu_meas += -cosmo->Theta2[2]*c_d - cosmo->beta_d*c_d;

	       /* Approximation: linear law */
	       //mu_meas += -cosmo->Theta2[2]*2.0*1e-2*sn->data[i].z
	       //  -3.0*2.0*1e-2*sn->data[i].z;

	    }

	    dmv = DistModVariance(&(sn->data[i]), Theta2);

	 } else {

	    /* beta = beta_0 + beta_1*z */
	    mu_meas += cosmo->Theta2[1]*(sn->data[i].s-1.0) +
	      (cosmo->Theta2[2] + cosmo->Theta2[3]*sn->data[i].z)*sn->data[i].c;
	    Theta2_betaz[0] = Theta2[0];
	    Theta2_betaz[1] = Theta2[1];
	    Theta2_betaz[2] = Theta2[2] + Theta2[3]*sn->data[i].z;
	    /* Include beta_z in the error */
	    dmv             = DistModVariance(&(sn->data[i]), Theta2_betaz);
	 }

      } else {

	 /* No stretch and color */
	 dmv = sn->data[i].cov[0][0];

      }

      /* Weight = inverse variance */
      wi  = 1.0/(sigintsqr + dmv + dsqr(sn->sig_mu_pec_vel/sn->data[i].z));
		 //+ dsqr(0.093*sn->data[i].z)); // Lensing correction (Kowalski08)

      if (wTheta1) {
	 /* Take into account Theta1 (zero-points) */
	 for (j=0,kTheta1=0.0; j<NDER; j++) kTheta1 += cosmo->Theta1[j]*sn->data[i].derivative[j];
      } else {
	 kTheta1 = 0.0;
      }

      chi2    += wi*dsqr(mu_meas + kTheta1 - sn->data[i].mu_c);
      /* |C| = Product_i [ w_i^{-1} ]; log|C| = - Sum_i [ w_i ] */
      logdetC -= log(wi);

   }


   /* Add up everything */
   res_all = -0.5*(sn->Nsample*ln2pi + chi2);

   /* Adding the covariance-term is correct for the log-likelihood in a Bayesian *
    * context. However, it leads to a biased estimation of parameters, in        *
    * particular for alpha and beta (J. Guy, P. Astier, priv. comm.)		 */
   if (add_logdetCov==1) {
      res      = -0.5*logdetC;
      res_all += res;
   }

   /* Zero-points (Kilbinger et al. 2009) */
   if (wTheta1==1) {
      /* TODO: use mvdens_log_pdf */
      testErrorRet(sn->W1==NULL, io_null, "Matrix for zero-points W1 not initialised",
		   *err, __LINE__, 0.0);
      chi2     = vect_scalar_product(cosmo->Theta1, sn->W1, sn->W1dim);
      res      = -0.5*(sn->W1dim*ln2pi + chi2 + sn->logdetW1);
      res_all += res;
   }

   /* Dust (Menard, Kilbinger & Scranton 2009) */
   /* In case of beta_d=const=data_beta_d->mean[0] this term is not added */
   if (cosmo->chi2mode==chi2_dust && data_beta_d!=NULL) {
      res = mvdens_log_pdf(data_beta_d, &(cosmo->beta_d), err);
      forwardError(*err, __LINE__, 0.0);
      res_all += res;
   }

   return res_all;
}

/* SNIa distance modulus. To make it independent of the Hubble constant, the *
 * H_0-dependence from w(z) is divided out -> 5 log(d_L/h) + 25.	     */
double distance_module(cosmo *self, double dlum, error **err)
{
   testErrorRet(dlum<0, ce_negative, "Luminosity distance is negative", *err, __LINE__, -1);
   return 5.0*log10(dlum/self->h_100) + 25;
}
