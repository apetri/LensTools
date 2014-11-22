/* ============================================================ *
 * cmb_bao.c							*
 * Martin Kilbinger 2008					*
 * Contains the CMB 'distance priors' (Komatsu et al.2008, K08) *
 * and BAO distance ratios (e.g. Eisenstein et al. 2005, E05,   *
 * Percival et al. 2007, P07).					*
 * ============================================================ */

#include "cmb_bao.h"

/* ============================================================ *
 * Hu & Sugiyama (1996) fitting formula for the decoupling red- *
 * shift, see also K08 (66-68).					*
 * ============================================================ */
double z_star(cosmo *self)
{
   double g1, g2, ob, om, h2, zs;

   h2 = self->h_100*self->h_100;
   ob = self->Omega_b*h2;
   om = self->Omega_m*h2;

   g1 = 0.0783*pow(ob, -0.238)/(1.0 + 39.5*pow(ob, 0.763));
   g2 = 0.560*(1.0 + 21.1*pow(ob, 1.81));
   zs = 1048.0*(1.0 + 0.00124*pow(ob, -0.738))*(1.0 + g1*pow(om, g2));

   return zs;
}

/* K08 (65) */
double acoustic_scale(cosmo *self, error **err)
{
   double zs, l_A, ws, as, r_s;

   zs   = z_star(self);
   as   = 1.0/(1.0+zs);
   ws   = w(self, as, 1, err);     forwardError(*err, __LINE__, -1.0);

   /* f_K is comoving angular diameter distance, therefore no (1+z_*) prefactor */
   l_A  = pi*f_K(self, ws, err);   forwardError(*err, __LINE__, -1.0);
   r_s  = r_sound_integral(self, as, err);  forwardError(*err, __LINE__, -1.0);
   l_A /= r_s;

   return l_A;
}

/* K08 (69) */
double shift_parameter(cosmo *self, error **err)
{
   double R, zs, ws, as;

   zs = z_star(self);
   as = 1.0/(1.0+zs);
   ws = w(self, as, 0, err);                      forwardError(*err, __LINE__, -1.0);
   R  = sqrt(self->Omega_m)*f_K(self, ws, err);   forwardError(*err, __LINE__, -1.0);
   R /= R_HUBBLE;

   return R;
}

/* Spherically average distance, E05, eq. (2), in units of Mpc/h */
double D_V(cosmo *self, double a, error **err)
{
   double DV, ww, fK, EE;

   ww  = w(self, a, 0, err);                     forwardError(*err, __LINE__, -1.0);
   fK  = f_K(self, ww, err);			 forwardError(*err, __LINE__, -1.0);
   EE  = Esqr(self, a, 0, err);			 forwardError(*err, __LINE__, -1.0);

   DV  = fK*fK/sqrt(EE)*R_HUBBLE*(1.0/a - 1.0);

   return cbrt(DV);
}

/* ============================================================ *
 * Log-likelihood for the CMB distance priors.			*
 * ============================================================ */
double chi2_cmbDP(cosmo *model, mvdens *g, error **err)
{
   double *x, c2;
   size_t ndim;

   ndim = g->ndim;

   testErrorRet(ndim!=3 && ndim!=4, ce_npar,
		"Number of dimensions in likelihood for CMB distance prior has to be 3 or 4",
		*err, __LINE__, 0.0);

   x = malloc_err(sizeof(double)*ndim, err);          forwardError(*err, __LINE__, 0.0);

   //dump_param(model, stderr);

   x[0] = acoustic_scale(model, err);                 forwardError(*err, __LINE__, 0.0);
   x[1] = shift_parameter(model, err);                forwardError(*err, __LINE__, 0.0);
   x[2] = z_star(model);

   if (ndim==4) {
      /* 100*omega_b */
      x[3] = 100.0*model->Omega_b*model->h_100*model->h_100;
   }

   c2 = mvdens_log_pdf(g, x, err);
   forwardError(*err, __LINE__, 0.0);

   free(x);

   return c2;
}

double chi2_bao_A(cosmo *model, mvdens *g, const double *z_BAO, error **err)
{
   double *x, c2;
   int i, ndim;

   ndim = g->ndim;

   x = malloc_err(sizeof(double)*ndim, err);    forwardError(*err, __LINE__, 0.0);

   /* BAO data ("A" parameter, Eisenstein 2005 eq 4): *
    * x = D_V(z)/cz * sqrt(Omega_m) * H_0             */

   for (i=0; i<ndim; i++) {
      x[i]  = D_V(model, 1.0/(z_BAO[i]+1.0), err);  forwardError(*err, __LINE__, 0.0);
      x[i] *= sqrt(model->Omega_m)/z_BAO[i]/R_HUBBLE;
   }

   c2 = mvdens_log_pdf(g, x, err);      forwardError(*err, __LINE__, 0.0);
   free(x);

   testErrorRet(!finite(c2), ce_infnan, "inf or nan in logl", *err, __LINE__, 0.0);

   return c2;
}

double chi2_bao_d_z(cosmo *model, mvdens *g, const double *z_BAO, error **err)
{
   double *d_z, c2, r_s;
   int i, ndim;

   ndim = g->ndim;

   d_z = malloc_err(sizeof(double)*ndim, err);    forwardError(*err, __LINE__, 0.0);

   //dump_param(model, stderr);

   /* BAO data (Percival et al. 2009): d_z = r_sound(r_d)/D_V(z) */

   /* Corrections of r_s to match Percival 09 r_s = 111.426 Mpc/h.   *
    * With fnu!=0 in z_drag, r_sound_integral() gives r_s = 112.091. *
    * With fnu=0, we get 111.284, within 0.13% of Percival 09.       */

   double z_d, a_d;
   z_d   = z_drag(model);
   a_d   = 1.0/(1.0+z_d);
   r_s   = r_sound_integral(model, a_d, err); forwardError(*err, __LINE__, -1.0);

   //r_s   = r_sound_drag_analytical(model, err);

   //r_s = r_sound_drag_fit(model, err);

   //r_s = 153.5 * model->h_100;
   //double h2;
   //h2   = model->h_100 * model->h_100;
   //r_s *= pow(model->Omega_m*h2/0.1326, -0.255); //*pow(model->Omega_b*h2/0.02273, -0.134);

   for (i=0; i<ndim; i++) {
      d_z[i] = r_s/D_V(model, 1.0/(z_BAO[i]+1.0), err);
      forwardError(*err, __LINE__, -1.0);
   }

   c2 = mvdens_log_pdf(g, d_z, err);     forwardError(*err, __LINE__, 0.0);
   free(d_z);

   testErrorRet(!finite(c2), ce_infnan, "inf or nan in logl", *err, __LINE__, 0.0);

   return c2;
}

double chi2_bao_D_V_ratio(cosmo *model, mvdens *g, const double *z_BAO, error **err)
{
   double *R, c2;
   int i, ndim;

   ndim = g->ndim;

   R = malloc_err(sizeof(double)*ndim, err);    forwardError(*err, __LINE__, 0.0);

   //dump_param(model, stderr);

   /* BAO data (Percival et al. 2007, 2009): R = D_V(z_1)/D_V(z_2) */

   for (i=0; i<ndim; i++) {
      R[i] = D_V(model, 1.0/(z_BAO[2*i]+1.0), err);
      forwardError(*err, __LINE__, -1.0);
      R[i] = R[i]/D_V(model, 1.0/(z_BAO[2*i+1]+1.0), err);
      forwardError(*err, __LINE__, -1.0);

      //printf("%f %f %g\n", z_BAO[2*i], z_BAO[2*i+1], R[i]);
   }

   c2 = mvdens_log_pdf(g, R, err);     forwardError(*err, __LINE__, 0.0);
   free(R);

   testErrorRet(!finite(c2), ce_infnan, "inf or nan in logl", *err, __LINE__, 0.0);

   return c2;
}
