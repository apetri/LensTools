/* ============================================================ *
 * math.h							*
 * Martin Kilbinger 06/2008					*
 * ============================================================ */

#ifndef __MATHS_H
#define __MATHS_H

#ifndef _NOFFTW_
#include <fftw3.h>
#endif

#ifdef __PLANCK__
#include "HL2_likely/tools/errorlist.h"
#include "HL2_likely/tools/io.h"
#include "HL2_likely/tools/maths_base.h"
#else
#include "maths_base.h"
#include "errorlist.h"
#include "io.h"
#endif

#define NR_END 1
#define FREE_ARG char*

/* See cosmo.h */
typedef enum {comp_c=0, comp_b=1, comp_nu=2} comp_t;
#define NCOMP 3

/* The following scales are Used for Hankel transforms */
#define theta_min 3.0e-7
#define theta_max 0.12

#define l_min     0.0001
#define l_max     5.0e6

/* Size of fast Hankel transform arrays (from lensing.h) */
#define N_thetaH  1024
typedef enum {tp_xipm, tp_gsqr, tp_map2_poly, tp_map2_gauss, tp_w, tp_xir} tpstat_t;


#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define MIN(a,b) ((a)<(b) ? (a) : (b))
#define MAX(a,b) ((a)>(b) ? (a) : (b))

/* Acceleration epoch (conservative, z=0.5) */
#define a_acc (2.0/3.0)


typedef struct {
  double* table;
  double a, b, dx, lower, upper;
  int n;
} interTable;

typedef struct  {
  double** table;
  double a1,b1,dx1,a2,b2,dx2,lower,upper;
  int n1,n2;
} interTable2D;

typedef struct {
  /* These tables start counting with 1 ! */
  double *x, *y, **z[NCOMP], **z2[NCOMP];
  int m, n;
} interTable2Dspline;

typedef struct {
  double *x, *y, **z[NCOMP];
  int m, n;
  unsigned long kx, ky;
} interTable2Dneq;

typedef struct {
  double *x, *y, *y2;
  double yp1, ypn;
  int n;
} splineTable;

/* Growth factor */
#define KMAXX 8
#define IMAXX (KMAXX+1)
typedef struct {
  int first, kmax, kopt;
  double epsold;
  double a[IMAXX + 1], alf[KMAXX + 1][KMAXX + 1];
} rkdrive_var_t;


/* ============================================================ *
 * Function types                                               *
 * ============================================================ */

typedef double (*funcwithpars)(double, void*, error **err);
typedef void (*derivative)(double a, double *y, double *yp, void* extra, error **err);
typedef double (*lgpdf)(void*,double*,error**);
typedef double (*lgpdf_c)(void*, const double*, error**);
typedef void (*rkdrive)(double [], double [], int, double *, double,
			double, double [], double *, double *, void *,
			derivative, rkdrive_var_t *, error **);

double nd_dfridr(lgpdf func, int idim, double *pars, double h,void* self,double *errn,  error **err);
double d2_func_dx_dy(lgpdf_c func, int idim, int jdim, double hh,double hx,double hy, double *x, double x0, double y0, void *self, error **err);
double nd_dfridr2(lgpdf_c func, int idim, int jdim, double *pars, double hx,double hy, void* self,
		  double *errn,  error **err);

double sm2_inverse(double *C, int N, error **err);
double *copy_matrix(double *C, int N, error **err);
double *multiply_matrix(const double *A, const double *B, int N, error **err);
void write_matrix(const double *A, int N, const char name[], error **err);
double *corr_coeff(double *C, int N, error **err);
void sm2_ludcmp(double *a, int n, int *indx, double *d, error **err);
void sm2_lubksb(double *a, int n, int *indx, double b[]);
double scalar_product_mk(const double *mean, const double *L, const double *x, int ndim, error **err);
void multiply_all(double *m, int n, double c);


#ifndef __APPLE__
double fmin(double maxarg1,double maxarg2);
double fmax(double maxarg1,double maxarg2);
#endif
double dsqr(double a);
double DCUB(double a);


double sinc(double x);
double gammln(double xx);
double bessj0(double x);


double **sm2_matrix(long nrl, long nrh, long ncl, long nch,error **err);
void sm2_free_matrix(double **m, long nrl, long nrh, long ncl, long nch);
double *sm2_vector(long nl, long nh,error **err);
void sm2_free_vector(double *v, long nl, long nh);
void sm2_polint(double xa[], double ya[], int n, double x, double *y, double *dy, error **err);
void sm2_spline(double x[], double y[], int n, double yp1, double ypn, double y2[], error **err);
void sm2_splint(double xa[], double ya[], double y2a[], int n, double x, double *y, error **err);
void sm2_splin2(double x1a[], double x2a[], double **ya, double **y2a, int m, int n,
                double x1, double x2, double *y, error **err);
splineTable *init_splineTable(int n, error **err);
splineTable *copy_splineTable(const splineTable *self, error **err);
void del_splineTable(splineTable** self);

void sm2_hunt(double xx[], unsigned long n, double x, unsigned long *jlo);

/* 1d integration */
double sm2_trapzdberg(funcwithpars func, void *intpar,
                      double a, double b, int n, double *s, error **err);
double sm2_qromberg(funcwithpars, double *intpar,
                    double a, double b, double EPS, error **err);
double sm2_qrombergo(funcwithpars  func, void *intpar,
                     double a, double b,
                     double (*choose)(double(*)(double,void *,error**), void *, double, double, 
				      int, double *,error**), double EPS,error **err);
double sm2_midpntberg(funcwithpars func, void *intpar,
                      double a, double b, int n, double *s,error **err);

/* ODEs */
void init_rkdrive_var(rkdrive_var_t *rkvar);
void bsstep(double y[], double dydx[], int nv, double *xx, double htry, double eps,
            double yscal[], double *hdid, double *hnext, void *extra,
            derivative, rkdrive_var_t *, error **);
void mmid(double y[], double dydx[], int nvar, double xs, double htot,
          int nstep, double yout[], void *ce,
          derivative,error **err);
void pzextr(int iest, double xest, double yest[], double yz[], double dy[], int nv,double*,
	    double**, error **err);
void odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1,
            void *ce, double hmin, int *nok, int *nbad,
            derivative, rkdrive, rkdrive_var_t *, error **);

/* Some complex functions */
void sm2_Cmul(my_complex a, my_complex b, my_complex *c);
void sm2_Cdiv(my_complex a, my_complex b, my_complex *c);
void cdgamma(my_complex x, my_complex *res);

/* ============================================================ *
 * Interpolation.						*
 * ============================================================ */

interTable* init_interTable(int n, double a, double b, double dx, double lower, double upper,
			    error **err);
interTable* copy_interTable(interTable* self, error **err);
void del_interTable(interTable** self);

interTable** init_interTable_arr(int N, int n, double a, double b, double dx, double lower, double upper, error **err);
interTable** copy_interTable_arr(interTable** self, int N, error **err);
void del_interTable_arr(interTable*** self, int N);


interTable2D* init_interTable2D(int n1, double a1, double b1, double dx1, int n2, double a2, double b2,
				double dx2, double lower, double upper, error **err);
interTable2D* copy_interTable2D(interTable2D* self, error **err);
void del_interTable2D(interTable2D** self);

interTable2D** init_interTable2D_arr(int N, int n1, double a1, double b1, double dx1, int n2, double a2, double b2, 
				     double dx2, double lower, double upper, error **err);
interTable2D** copy_interTable2D_arr(interTable2D **self, int N, error **err);
void del_interTable2D_arr(interTable2D*** self, int N);

interTable2Dspline *init_interTable2Dspline(int m, int n, error **err);
interTable2Dspline *copy_interTable2Dspline(interTable2Dspline *source, error **err);
void del_interTable2Dspline(interTable2Dspline **self);

interTable2Dneq *init_interTable2Dneq(int m, int n, error **err);
interTable2Dneq *copy_interTable2Dneq(interTable2Dneq *source, error **err);
void del_interTable2Dneq(interTable2Dneq **self);

double interpol_wr(interTable* self, double x, error **err);
double interpol2D(interTable2D* self, double x, double y, error **err);
double sm2_interpol(double *f, int n, double a, double b, double dx, double x,
		    double lower, double upper, error **err);
double sm2_interpol2d(double **f, int nx, double ax, double bx, double dx, double x,
		      int ny, double ay, double by, double dy, double y,
		      double lower, double upper, error **err);
double sm2_interpol2Dneq(interTable2Dneq *self, comp_t comp, double x0, double y0, error **err);


void jacobi_transform(const double *ainput, int n, double d[], double **v, int *nrot, error **err);
void indexx(unsigned long n, double arr[], unsigned long indx[], error **err);

#ifndef _NOFFTW_
/* Fast Hankel transform */
void tpstat_via_hankel(void*, double **xi, double *logthetamin, double *logthetamax,
		       tpstat_t tpstat, double (*P_projected)(void *, double, int, int, error**),
		       int i_bin, int j_bin, error **err);
void hankel_kernel_mu(double x, fftw_complex *res, double q, double mu, error **err);
void hankel_kernel_mumu(double x, fftw_complex *res, double q, double mu, error **err);
void hankel_kernel_exp(double k, fftw_complex *res, error **err);
void hankel_kernel_tophat(double k, fftw_complex *res, error **err);
#endif

#endif
