/* ============================================================ *
 * decomp_eb.c							*
 * Martin Kilbinger, Liping Fu 2008, 2009			*
 * ============================================================ */
#include "decomp_eb.h"
#include <gsl/gsl_sf_gamma.h>

/* === CFHTLS Wide 3rd data release, Fu&Kilbinger (2010) === */

/* S/N, Psi=19' eta=1/50 */
const double a_FK10_SN[N_FK10]        = {0.1197730890, -0.3881211865, 0.5212557875, -0.3440507036, 0.2761305382, -0.07286690971};

/* FoM, Psi=222', eta=1/10 */
const double a_FK10_FoM_eta10[N_FK10] = {0.009877788826, 0.1061397843, -0.4300211814, 0.5451016406, -0.3372272549, 0.1716983151};

/* FoM, Psi=222', eta=1/50 */
const double a_FK10_FoM_eta50[N_FK10] = {0.1239456383, -0.3881431858, 0.5579593467, -0.3679282338, 0.1540941993, 0.01293361618};


/* First-kind Chebyshev polynomial of order n */
#define EPS 1.0e-10
double Cheby(double x, int n, error **err)
{
   double Cn;

   testErrorRetVA(x<-1-EPS || x>1+EPS, mr_range, "x = %g out of range", *err, __LINE__, -1, x);

   if (x<-1) x = -1.0;
   if (x>1) x = 1.0;

   Cn = cos (n * acos (x));
   return Cn;
}
#undef EPS

/* Second-kind Chebyshev polynomial of order n */
#define EPS 1.0e-10
double Cheby2(double x, int n, error **err)
{
   double Un;

   testErrorRetVA(x<-1-EPS || x>1+EPS, mr_range, "x = %g out of range", *err, __LINE__, -1, x);

   if (x<-1) x = -1.0;
   if (x>1) x = 1.0;

   if (x == 1) {
      Un = n + 1;
   } else if (x == -1) {
      Un = pow(-1.0, n) * (n + 1.0);
   } else {
      Un = sin((n+1.0)*acos(x))/sin(acos (x));
   }

   return Un;
}
#undef EPS

/* Legendre polynomial of order n */
double Legen(double x, int n)
{

   switch (n) {
      case 0:
	 return 1.0;
      case 1:
	 return x;
      default:
	 return (2.0*n-1.0)/(double)n*x*Legen(x,n-1)-(n-1.0)/(double)n*Legen(x,n-2);  
   }
}

/* General basis function of order n */
#define EPS 1.0e-10
double C(double x, int n, poly_t poly, error **err)
{
   double c;

   testErrorRetVA(x<-1-EPS || x>1+EPS, mr_range, "x = %g out of range", *err, __LINE__, -1, x);

   if (x<-1) x = -1.0;
   if (x>1) x = 1.0;

   switch (poly) {
      case cheby :
	 c = Cheby(x, n, err);
	 forwardError(*err, __LINE__, -1.0);
	 break;
      case cheby2 :
	 c = Cheby2(x, n, err);
	 forwardError(*err, __LINE__, -1.0);
	 break;
      case legen :
	 c = Legen(x, n);
	 break;
      default :
	 *err = addErrorVA(mr_poly, "Unknown polynomial type %d", *err, __LINE__, poly);
	 return -1.0;
   }

   return c;
}
#undef EPS


/* ============================================================ *
 * The filter function for the generalised ring statistics.	*
 * FK09 (11).							*
 * ============================================================ */
#define	 EPS 1.0e-6
double Tp(double x, const double *a, int N, poly_t poly, error **err)
{
   int n;
   double res, Cn, Cnm1=0.0, Cnm2; 

   if (x<-1-EPS || x>+1+EPS) return 0;

   /* NEW! */
   //if (x<-1) x = -1;
   //if (x>+1) x = +1;


   testErrorRetVA(N<=0, mr_range, "N has to be larger than zero but is %d", *err, __LINE__, 0.0, N);

   testErrorRetVA(x>+1+EPS, mr_range, "x=%.10f out of range", *err, __LINE__, 0.0, x);
   testErrorRetVA(x<-1-EPS, mr_range, "x=%.10f out of range", *err, __LINE__, 0.0, x);

   if (poly==cheby || poly==cheby2) {

      Cnm2 = 1.0;
      if (poly==cheby)       Cnm1 = x;
      else if (poly==cheby2) Cnm1 = 2.0*x;

      res  = 0.0;
      res += a[0] * Cnm2;
      if (N==1) return res;
      res += a[1] * Cnm1;
      if (N==2) return res;
	
      for (n=2; n<N; n++) {
	 /* Chebyshev (both T,U) recurrence relation *
	  * C_n(x) = 2 x C_{n-1}(x) - C_{n-2}(x)     */
	 Cn   = 2.0*x*Cnm1 - Cnm2;
	 res += a[n] * Cn;
	 Cnm2 = Cnm1;
	 Cnm1 = Cn;
      }

   } else {

      for (n=0,res=0.0; n<N; n++) {
	 Cn   = C(x, n, poly, err);
	 forwardError(*err, __LINE__, -1);
	 res += a[n] * Cn;
      }

   }

   return res;
}
#undef EPS

/* ============================================================ *
 * FK09 (24).							*
 * ============================================================ */
double Fn0(double x, int n, poly_t poly, error **err)
{
   double fn0, denom=0.0;
   int nn;


   if (n<0) nn = -n;
   else nn = n;

   if (poly==cheby) {
      if (nn==0) return x + 1.0;
      if (nn==1) return 0.5*(x*x - 1.0);
   } else if (poly==cheby2) {
      if (n==-1) return 0.0;  /* n, not nn: no symmetry for cheby2! */
   } else {
      *err = addErrorVA(mr_poly, "Tm for polynomial type %d not defined", *err, __LINE__, poly);
      return 0.0;
   }

   fn0 = 0.0;

   if (nn%2==0) fn0 += 1.0;   /* n even */
   else fn0 -= 1.0;           /* n odd  */

   //printf("# ** 1 %g\n", fn0);

   fn0 += x*Cheby(x, nn, err);
   forwardError(*err, __LINE__, 0.0);

   //printf("# ** 2 %g\n", fn0);

   if (poly==cheby) fn0 += n*(1 - x*x)*Cheby2(x, n-1, err);
   else if (poly==cheby2) fn0 -= (1 - x*x)*Cheby2(x, n-1, err);
   forwardError(*err, __LINE__, 0.0);

   //printf("# ** 3 %g\n", fn0);

   if (poly==cheby) denom = 1.0 - (double)(n*n); 
   else if (poly==cheby2) denom = 1.0+(double)n;         /* n, not nn! */

   fn0 /= denom;

   return fn0;
}

/* ============================================================ *
 * FK09 (26)							*
 * ============================================================ */
void Fnnu(double x, int n, poly_t poly, double Fn[], error **err)
{
   double f0, fp1, fp2, fp3, fm1, fm2, fm3;

   f0  = Fn0(x, n, poly, err);      forwardError(*err, __LINE__,);
   fp1 = Fn0(x, n+1, poly, err);    forwardError(*err, __LINE__,);
   fp2 = Fn0(x, n+2, poly, err);    forwardError(*err, __LINE__,);
   fp3 = Fn0(x, n+3, poly, err);    forwardError(*err, __LINE__,);
   fm1 = Fn0(x, n-1, poly, err);    forwardError(*err, __LINE__,);
   fm2 = Fn0(x, n-2, poly, err);    forwardError(*err, __LINE__,);
   fm3 = Fn0(x, n-3, poly, err);    forwardError(*err, __LINE__,);

   Fn[0] = f0;
   Fn[1] = 0.5*(fp1 + fm1);
   Fn[2] = 0.25*(fp2 + 2.0*f0 + fm2);
   Fn[3] = 0.125*(fp3 + 3.0*fp1 + 3.0*fm1 + fm3);
}

/* ============================================================ *
 * FK09 (23)							*
 * ============================================================ */
#define EPS 1.0e-10
double alpha_explicit(double x, int n, double R, poly_t poly, error **err)
{
   double r, Fn[4], alpha;

   testErrorRet(poly!=cheby && poly!=cheby2, mr_poly,
		"Explicit expression only obtained so far for cheby1 and cheby2",
		*err, __LINE__, 0.0);

   if (x<-1.0+EPS) return 0.0;
   assert(x>-1.0);

   r = 3.0/dsqr(x + R);
   Fnnu(x, n, poly, Fn, err);
   forwardError(*err, __LINE__, 0.0);

   alpha  = R*(1.0-r*R*R)*Fn[0];
   alpha += (1.0 - 3.0*r*R*R)*Fn[1];
   alpha += -3.0*r*R*Fn[2];
   alpha += -r*Fn[3];
   alpha *= 4.0/3.0*r;

   return alpha;
}

/* ============================================================ *
 * The filter function T_- obtained from T_+, FK09 (20).	*
 * ============================================================ */
#define EPS 1.0e-10
double Tm(double x, const double *a, int N, poly_t poly, double R, error **err)
{
   double tm, c, alph;
   int n;

   testErrorRetVA(x<-1-EPS, mr_range, "x=%.20f smaller than -1", *err, __LINE__, 0.0, x);

   for (n=0,tm=0.0; n<N; n++) {
      c    = C(x, n, poly, err);                forwardError(*err, __LINE__, 0.0);
      alph = alpha_explicit(x, n, R, poly, err);
      forwardError(*err, __LINE__, 0.0);
      
      tm  += a[n]*(c + alph);
   }
   return tm;
}
#undef EPS

/* ================================================================== *
 * Returns the generalised ring statistic using a as the coefficients *
 * of the T_+ decomposition, see FK09 (4, 8, 11).		      *
 * The correlation functions xi+ and xi- are read from the            *
 * Nxi-dimensional arrays xip and xim, on angular scales th.          *
 * N is the size of the coefficient-array a.                          *
 * ================================================================== */
double RR_data(const double *xip, const double *xim, const double *th, const int Nxi,
	       double THETA_MIN, double THETA_MAX, const double *a, int N,
	       poly_t poly,  int pm, error **err)
{
   double theta, dtheta, res, A, B, x, summand;
   double dt1, dt2;
   int i;

   testErrorRetVA(abs(pm)!=1, mr_incompatible, "pm=%d not valid, has to be +1 or -1",
		  *err, __LINE__, 0.0, pm);
   
   testErrorRetVA(THETA_MIN<th[0], mr_range,
		  "THETA_MIN=%g' is smaller than minimum angular scale for xi+-, %g'",
		  *err, __LINE__, 0.0, THETA_MIN/arcmin, th[0]/arcmin);

   testErrorRetVA(THETA_MAX>th[Nxi-1], mr_range,
		  "THETA_MAX=%g' is larger than maximum angular scale for xi+-, %g'",
		  *err, __LINE__, 0.0, THETA_MAX/arcmin, th[Nxi-1]/arcmin);

   if (a==NULL) {
      *err = addError(mr_null, "Coefficients a=NULL", *err, __LINE__);
      return 0.0;
   }


   /* FK09 (8) */
   A = (THETA_MAX - THETA_MIN)/2.0;
   B = (THETA_MAX + THETA_MIN)/2.0;

   for (i=0,res=0.0; i<Nxi; i++) {
      theta  = th[i];
      if (theta < THETA_MIN) continue;
      else if (theta > THETA_MAX) break;

      /* theta[i] is the bin center */

      if (i==0) {
	 dt1 = (th[i+1] - th[i])/2.;
	 dtheta = 2*dt1;
      } else if (i==Nxi-1) {
	 dt2 = (th[i] - th[i-1])/2.;
	 dtheta = 2*dt2;
      } else { 
	 dt1 = (th[i+1] - th[i])/2.;
	 dt2 = (th[i] - th[i-1])/2.;
	 dtheta = dt2 + dt1;
      }
      x = (theta-B)/A;
      
      if (pm==+1)  {  
	 summand  = xip[i]*theta/dsqr(THETA_MAX);
	 summand *= Tp(x, a, N, poly, err);
	 forwardError(*err, __LINE__, -1.0);
      } else if (pm==-1) {
	 summand  = xim[i]*theta/dsqr(THETA_MAX);
	 summand *= Tm(x, a, N, poly, B/A, err);
	 forwardError(*err, __LINE__, -1.0);
      } else {
	 summand = 0.0;
      }	     

      res += summand*dtheta;
   }

   testErrorRet(!finite(res), math_infnan, "R is not finite", *err, __LINE__,
		  0.0);

   return res;
}


/* If cov_mode=outer, n, m are not used. For cov_mode=inner, a, N are not used.  *
 * If cov_mode=fromZ, the covariance using Z+ is returned and a, N are not used. */
#define EPS 1.0e-5
double cov_RR(const double *THETA_MIN, const double *THETA_MAX, const double *a, int N, poly_t poly,
	      const double *theta, const double *cov_xi, int Ntheta, 
	      cov_mode_t cov_mode, int n, int m, double fac, error **err)
{
   double sum, dlogtheta, tmp, A[2], B[2], xi, xj, yi, yj, thetai, thetaj;
   int i, j;
   static double  **Cov_xi = NULL;


   /* Testing for a==NULL not possible because cov_RR is also called from get_cov_RR_inner_array */
   if (cov_mode==fromZ) {
#ifdef __MRING_H
      sum = cov_RR_Z(THETA_MIN, THETA_MAX, theta, cov_xi, Ntheta, err);
      forwardError(*err, __LINE__, 0.0);
      return sum;
#else
      *err = addError(mr_type, "cov_mode 'fromZ' not supported here", *err, __LINE__);
      return 0.0;
#endif
   }

   for (i=0; i<2; i++) {
      A[i] = (THETA_MAX[i] - THETA_MIN[i])/2.0;
      B[i] = (THETA_MAX[i] + THETA_MIN[i])/2.0;
   }

   dlogtheta = log(theta[1]) - log(theta[0]);

   for (i=0; i<2; i++) {
      testErrorRetVA(theta[0]>THETA_MIN[i], mr_range,
		     "Minumum of xi-covariance (%g=%g') larger than THETA_MIN (%g=%g')",
		     *err, __LINE__, -1, theta[0], theta[0]/arcmin, THETA_MIN[i], THETA_MIN[i]/arcmin);
      testErrorRetVA(theta[Ntheta-1]<THETA_MAX[i], mr_range,
		     "Maximum of xi-covariance scale (%g=%g') smaller than THETA_MAX (%g=%g')",
		     *err, __LINE__, -1, theta[Ntheta-1], theta[Ntheta-1]/arcmin, THETA_MAX[i],
		     THETA_MAX[i]/arcmin);
   }

   if (Cov_xi==NULL) {
      Cov_xi = sm2_matrix(0, Ntheta-1, 0, Ntheta-1, err);
      forwardError(*err, __LINE__, -1.0);
      for (i=0; i<Ntheta; i++) {
	 for (j=0; j<Ntheta; j++) {

	    Cov_xi[i][j] = cov_xi[i*Ntheta+j];

	 }

      }

      //if (fac<1.12)
	//fprintf(stderr, "WARNING: Cov-xi-factor = %g, too small for stable inversion of covariance!!!\n", fac);
   }

   sum = 0.0;
   //for (i=0; i<Ntheta; i++) {
   //thetai = theta[i];
   for (thetai=theta[0]; thetai<=theta[Ntheta-1]; thetai*=fac) {
      if (thetai<THETA_MIN[0] || thetai>THETA_MAX[0]) continue;
      xi = (thetai-B[0])/A[0];
      yi = thetai/THETA_MAX[0];
   
      //for (j=0; j<Ntheta; j++) {
      //thetaj = theta[j];
      for (thetaj=theta[0]; thetaj<=theta[Ntheta-1]; thetaj*=fac) {
	 if (thetaj<THETA_MIN[1] || thetaj>THETA_MAX[1]) continue;
	 xj = (thetaj-B[1])/A[1];
	 yj = thetaj/THETA_MAX[1];

	 tmp  = dsqr(dlogtheta);

	 if (cov_mode==outer) {
	    tmp *= Tp(xi, a, N, poly, err);
	    forwardError(*err, __LINE__, -1);
	    tmp *= Tp(xj, a, N, poly, err);
	    forwardError(*err, __LINE__, -1);
	 } else {
	    tmp *= C(xi, n, poly, err);  	    forwardError(*err, __LINE__, -1);
	    tmp *= C(xj, m, poly, err);  	    forwardError(*err, __LINE__, -1);
	 }

	 tmp *= dsqr(yi*yj);

	 //tmp *= cov_xi[i*Ntheta+j];
	 tmp *= sm2_interpol2d(Cov_xi, Ntheta, log(theta[0]), log(theta[Ntheta-1]), dlogtheta, log(thetai),
	 		       Ntheta, log(theta[0]), log(theta[Ntheta-1]), dlogtheta, log(thetaj), 0.0, 0.0, err);
	 forwardError(*err, __LINE__, -1.0);

	 sum += tmp;

      }
 
  }

   /* Comment the following two lines for faster code, but be aware that no check is done
    * whether the input covariance changes! */
   sm2_free_matrix(Cov_xi, 0, Ntheta-1, 0, Ntheta-1);
   Cov_xi = NULL;

   return sum;
}
#undef EPS

/* RR-Covariance using a diagonal xi-covariance.				 *
 * If cov_mode=outer, n, m are not used. For cov_mode=inner, a, N are not used.  *
 * If cov_mode=fromZ, the covariance using Z+ is returned and a, N are not used. */
#define EPS 1.0e-5
double cov_RR_diag_xi(const double *THETA_MIN, const double *THETA_MAX, const double *a, int N, poly_t poly,
		      const double *theta, const double *var_xi, int Ntheta, 
		      int islog, error **err)
{
   double sum, dlogtheta, dtheta, tmp, A[2], B[2], xi, xj, yi, yj, thetai, thetaj;
   int i;

   for (i=0; i<2; i++) {
      A[i] = (THETA_MAX[i] - THETA_MIN[i])/2.0;
      B[i] = (THETA_MAX[i] + THETA_MIN[i])/2.0;
   }

   if (islog==0) {
      dtheta = theta[1] - theta[0];
      dlogtheta = 0.0;
   } else if (islog==1) {
      dlogtheta = log(theta[1]) - log(theta[0]);
      dtheta = 0.0;
   } else {
      *err = addErrorVA(mr_type, "Invalid flag islog = %d\n", *err, __LINE__, islog);
      return 0.0;
   }

   for (i=0; i<2; i++) {
      testErrorRetVA(theta[0]>THETA_MIN[i], mr_range,
		     "Minumum of xi-covariance (%g=%g') larger than THETA_MIN (%g=%g')",
		     *err, __LINE__, -1, theta[0], theta[0]/arcmin, THETA_MIN[i], THETA_MIN[i]/arcmin);
      testErrorRetVA(theta[Ntheta-1]<THETA_MAX[i], mr_range,
		     "Maximum of xi-covariance scale (%g=%g') smaller than THETA_MAX (%g=%g')",
		     *err, __LINE__, -1, theta[Ntheta-1], theta[Ntheta-1]/arcmin, THETA_MAX[i],
		     THETA_MAX[i]/arcmin);
   }

   sum = 0.0;
   for (i=0; i<Ntheta; i++) {
      thetai = theta[i];
      thetaj = thetai;

      if (thetai<THETA_MIN[0] || thetai>THETA_MAX[0]) continue;
      xi = (thetai-B[0])/A[0];
      yi = thetai/THETA_MAX[0];
   
      xj = (thetaj-B[1])/A[1];
      yj = thetaj/THETA_MAX[1];

      if (islog==0) {
	 //	 tmp = dsqr(dtheta) /thetai/thetai;
	 tmp = dsqr(dtheta)*yi*yj/THETA_MAX[0]/THETA_MAX[1];

      } else {
	
	 tmp  = dsqr(dlogtheta);
	 tmp *= dsqr(yi*yj);
      }

      tmp *= Tp(xi, a, N, poly, err);
      forwardError(*err, __LINE__, -1);
      tmp *= Tp(xj, a, N, poly, err);
      forwardError(*err, __LINE__, -1);


      tmp *= var_xi[i];
      //      tmp *= dsqr(yi*yj);

      sum += tmp; 
  }

   return sum;
}
#undef EPS

/* ============================================================ *
 * Returns the chi^2 for a null test (RR=0).			*
 * ============================================================ */
double chi2_RB_null(const double *RB, const double *covRB, int NRB)
{
   int i, j;
   double c2;

   for (i=0,c2=0.0; i<NRB; i++) {
      for (j=0; j<NRB; j++) {
	 c2 += RB[i] * covRB[i*NRB+j] * RB[j];
      }
   }

   return c2;
}

/* ============================================================ *
 * Reads COSEBIs zeros and normalisation coefficients from file *
 * and returns polynomial coefficients.	File name is automatic, *
 * scale checks are done.					*
 * ============================================================ */
double *read_zeros_norm_cosebi_auto_check(double Psimin, double Psimax, const char *path, error **err)
{
   double *c_cosebi, psimin, psimax;
   char rname[1024];

   sprintf(rname, "%s/cosebi_tplog_rN_%d_%.3f_%.3f",
	   path == NULL ? "." : path, NMAX_COSEBI, Psimin/arcmin, Psimax/arcmin);
   c_cosebi = read_zeros_norm_cosebi(rname, &psimin, &psimax, err);
   forwardError(*err, __LINE__, NULL);

   testErrorRetVA(fabs(psimin-Psimin/arcmin)>EPSILON, mr_file,
		  "Inconsistent file %s, psimin=%g should be %g according to the file name",
		  *err, __LINE__, NULL, rname, psimin, Psimin/arcmin);
   testErrorRetVA(fabs(psimax-Psimax/arcmin)>EPSILON, mr_file,
		  "Inconsistent file %s, psimax=%g should be %g according to the file name",
		  *err, __LINE__, NULL, rname, psimax, Psimax/arcmin);

   return c_cosebi;
}

/* ============================================================ *
 * Reads COSEBIs zeros and normalisation coefficients from file *
 * and sets corresponding min and max scales.			*
 * Calculates polynomial coefficients c and returns them.       *
 * ============================================================ */
double *read_zeros_norm_cosebi(const char *rname, double *psimin, double *psimax, error **err)
{
   FILE *F;
   int Nzeros, Ncoeff, k, n, off_c, off_R, nmax;
   ssize_t nread;
   double *Rn, *Norm, *c;
   char *str, *line=NULL;


   F = fopen_err(rname, "r", err);                  forwardError(*err, __LINE__, NULL);

   /* Read two header lines */
   line = malloc_err(1024*sizeof(char), err);       forwardError(*err, __LINE__, NULL);
   str = fgets(line, 1024, F);
   str = fgets(line, 1024, F);
   free(line);

   nread = fscanf(F, "%d  %lg %lg\n", &nmax, psimin, psimax);
   testErrorRet(nread != 3, mr_file, "File has wrong format.", *err, __LINE__, NULL);

   testErrorRetVA(nmax > NMAX_COSEBI, mr_range,
		  "COSEBI number of modes n=%d read from file %s cannot be larger than NMAX_COSEBI=%d",
		  *err, __LINE__, NULL, nmax, rname, NMAX_COSEBI);


   /* Number of zeros for nmax polynomials,               *
    * sum_{i=1}^{nmax+1} (i+1) = sum_{i=0}^{nmax+1} i - 1 */
   Nzeros = (nmax+1) * (nmax+2) / 2 - 1;

   //fprintf(stderr, "nmax = %d , Psi = [%g, %g]\n", nmax, *psimin, *psimax);
   //fprintf(stderr, "Reading %d Rn, %d Norm, from %s\n", Nzeros, nmax, rname);

   Rn   = malloc_err(sizeof(double)*Nzeros, err);   forwardError(*err, __LINE__, NULL);
   Norm = malloc_err(sizeof(double)*nmax, err);     forwardError(*err, __LINE__, NULL);

   /* Read zeros */
   for (k=0; k<Nzeros; k++) {
      fscanf(F, "%lg ", Rn+k);
   }
   /* Read normalisations */
   for (k=0; k<nmax; k++) {
      fscanf(F, "%lg ", Norm+k);
   }

   fclose(F);


   /* Calculate polynomial coefficients */
   Ncoeff = NMAX_COSEBI * (NMAX_COSEBI + 5) / 2;
   c =  malloc_err(sizeof(double)*Ncoeff, err);        forwardError(*err, __LINE__, NULL);

   /* Polynomial T_i of order i+1 */
   for (n=1,off_R=0; n<=nmax; n++) {

      off_c = n*(n + 3)/2 - 2;   // Offset for coefficients for polynomial #n

      for (k=0; k<n+2; k++) {
	 c[k + off_c] = sum_combinations(n+1-k, n+1, Rn+off_R, err);
	 forwardError(*err, __LINE__, NULL);
	 c[k + off_c] *= Norm[n-1];
	 if (k % 2 == 1) c[k + off_c] = -c[k + off_c];
	 if (n % 2 == 0) c[k + off_c] = -c[k + off_c];
      }

      off_R = (n+1)*(n+2)/2 - 1; // Ncoeff for index i, degree i+1
   }

   free(Rn);
   free(Norm);

   return c;
}

/* ============================================================ *
 * COSEBI logarithmis filter function, SEK10 (28). The          *
 * coefficients have to be precalculated from the zeros, output *
 * by Mathematica.						*
 * Not normalised. 						*
 * ============================================================ */
double Tplog_c(double z, const double *c, int n, error **err)
{
   int i, off_c;
   double res, zpowerm;
   const double *cn;

   off_c = n*(n + 3)/2 - 2; // Offset for coefficients for polynomial #n
   cn    = c + off_c;

   for (i=0,res=0.0,zpowerm=1.0; i<=n+1; i++) {
      res     += cn[i] * zpowerm;
      zpowerm *= z;
   }

   testErrorRet(!finite(res), math_infnan, "Tplog is not finite", *err, __LINE__, 0.0);

   return res;
}

/* ============================================================ *
 * Filter function, SEK10 (38, 39). Normalised.                 *
 * ============================================================ */
double Tmlog(double z, const double *c, int n, error **err)
{
   int m, off_c;
   double res, zpowerm;
   const double *cn;

   off_c = n*(n + 3)/2 - 2; // Offset for coefficients for polynomial #n
   cn    = c + off_c;

   for (m=0,res=0.0,zpowerm=1; m<=n; m++) { // n+1 MKDEBUG ??
      res += dnm(n, m ,cn) * zpowerm;
      zpowerm *= z;
   }

   res += an2(n, cn) * exp(-2.0*z) - an4(n, cn) * exp(-4.0*z);
   testErrorRet(!finite(res), math_infnan, "Tmlog is not finite", *err, __LINE__, 0.0);

   return res;
}

double an2(int n, const double *c)
{
   int j;
   double res, jfac, m2powerjp1;
   
   for (j=0,res=0.0,jfac=1.0,m2powerjp1=-2.0; j<=n+1; j++) {
      res        += c[j]*jfac / m2powerjp1;
      jfac       *= j + 1.0;      /* factorial(j) */
      m2powerjp1 *= -2.0;        /* (-2)^(j+1)   */
   }

   res *= 4.0;

   return res;
}

double an4(int n, const double *c)
{
   int j;
   double res, jfac, m4powerjp1;
   
   for (j=0,res=0.0,jfac=1.0,m4powerjp1=-4.0; j<=n+1; j++) {
      res        += c[j]*jfac / m4powerjp1;
      jfac       *= j + 1.0;
      m4powerjp1 *= -4.0;        /* (-4)^(j+1) */
   }

   res *= 12.0;

   return res;
}

/* dnm(n+1, n, c) = 0 */
double dnm(int n, const int m, const double *c)
{
   int j;
   double res, jfac, m2powermjm1, p2powermjm1, mfac;

   p2powermjm1 = 0.5;
   m2powermjm1 = -0.5;
   mfac = gsl_sf_fact(m);
   for (j=m,res=0.0,jfac=mfac; j<=n+1; j++) {
      res         += c[j] * jfac * m2powermjm1 * (3.0 * p2powermjm1 - 1.0);
      jfac        *= j + 1.0;        /* j! */
      m2powermjm1 /= -2.0;    /* (-2)^(m-j-1) */
      p2powermjm1 /= 2.0;     /* 2^(m-k-1)    */
   }

   res *= 4.0/mfac;
   res += c[m];

   return res;
}

double sum_combinations(int j, int n, const double *r, error **err)
{
   gsl_combination *c;
   size_t i, c_i;
   double sum, prod;


   if (n==0) return 1.0;

   c = gsl_combination_calloc(n, j);
     
   sum = 0.0;
   do {
      for (i=0,prod=1.0; i<j; i++) {
	 c_i =  gsl_combination_get(c, i);
	 prod *= r[c_i];
	 //printf(" %zu", c_i);
      }
      sum += prod;
      //printf("{"); gsl_combination_fprintf(stdout, c, " %u"); printf (" }\n");
   } while (gsl_combination_next(c) == GSL_SUCCESS);

   gsl_combination_free(c);

   return sum;
}


/* ============================================================ *
 * E-/B-mode correlation functions from COSEBIs, SEK10 (40).    *
 * The results are stored in xi_pm_EB = (E+, E-, B+, B-).       *
 * ============================================================ */
void xipmEB(double theta, double THETA_MIN, double THETA_MAX, const double *c,
	    const double *E, const double *B, int N, double xi_pm_EB[4], error **err) 
{
   double Tp, Tm, z;
   int n;

   for (n=0; n<4; n++) {
      xi_pm_EB[n] = 0.0;
   }

   for (n=1; n<=N; n++) {
      z  = log(theta/THETA_MIN);
      Tp = Tplog_c(z, c, n, err); 
      Tm = Tmlog(z, c, n, err); 
	 
      xi_pm_EB[0] += E[n] * Tp;   // E+
      xi_pm_EB[1] += E[n] * Tm;   // E-
      xi_pm_EB[2] += B[n] * Tp;   // B+
      xi_pm_EB[3] += B[n] * Tm;   // B-
   }

   for (n=0; n<4; n++) {
      xi_pm_EB[n] *= 2.0 / theta / (THETA_MAX - THETA_MIN);
   }
}

