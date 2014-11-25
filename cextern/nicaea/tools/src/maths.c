/* ============================================================ *
 * maths.c							*
 * Martin Kilbinger 06/2008					*
 * ============================================================ */

#ifdef __PLANCK__
#include "HL2_likely/tools/maths.h"
#else
#include "maths.h"
#endif


#define CON 1.4
#define CON2 (CON*CON)
#define BIG 1.0e30
#define NTAB 10
#define SAFE 2.0
double nd_dfridr(lgpdf func,int idim, double *pars, double h, void* self, double *errn, error **err) {
  int i,j;
  double errt,fac,hh,*a,ans,*x,x0;
  double dval,qval;
  
  ans = 1e30; /* dummy initialization */
  testErrorRetVA(h==0.0, math_wrongValue, "h has to be larger than zero (got %g)", *err, __LINE__, 0.0,h);
  

  x = pars;
  x0 = x[idim];
  
  a = malloc_err(sizeof(double)*NTAB*NTAB,err);
  forwardError(*err,__LINE__,0);

  hh = h;
  
  x[idim] = x0+hh;
  
  dval = (*func)(self,x,err);
  forwardError(*err,__LINE__,0);
  //_DEBUGHERE_("F(%d-> %g + %g) = %15.10f",idim,x0,hh,dval);
  
  x[idim] = x0-hh;
  qval = (*func)(self,x,err);
  forwardError(*err,__LINE__,0);
  //_DEBUGHERE_("F(%d-> %g - %g) = %15.10f",idim,x0,hh,qval);
  
  a[0] = (dval-qval)/(2.0*hh);
  //_DEBUGHERE_("DF(%d-> %g : %d) = %15.10f",idim,x0,0,a[0]);
  
  *errn=BIG;
  for (i=1;i<NTAB;i++) {
    hh /= CON;
    
    x[idim] = x0+hh;
    dval = (*func)(self,x,err);
    forwardError(*err,__LINE__,0);
    //_DEBUGHERE_("F(%d-> %g + %g) = %15.10f",idim,x0,hh,dval);
    
    x[idim] = x0-hh;
    qval = (*func)(self,x,err);
    forwardError(*err,__LINE__,0);
    //_DEBUGHERE_("F(%d-> %g - %g) = %15.10f",idim,x0,hh,qval);

    a[i] = (dval-qval)/(2.0*hh);
    //_DEBUGHERE_("DF(%d-> %g : %d) = %15.10f",idim,x0,i,a[i]);

    fac = CON2;
    for (j=1;j<=i;j++) {
      a[j*NTAB+i] = (a[(j-1)*NTAB+i]*fac - a[(j-1)*NTAB+i-1]) / (fac-1.0);
      fac = CON2*fac;
      
      errt = fmax(fabs(a[j*NTAB+i] - a[(j-1)*NTAB+i]),fabs(a[j*NTAB+i] - a[(j-1)*NTAB+i-1]));
      if (errt <= *errn) {
        *errn = errt;
        ans = a[j*NTAB+i];
      }
    }
    
    if (fabs(a[i*NTAB+i]-a[(i-1)*NTAB+i-1]) >= SAFE*(*errn)) 
      break;
  }
  pars[idim]=x0;
  free(a);
  _DEBUGHERE_("DF(%d-> %g) = %15.10f +- %g (%d)",idim,x0,ans,*errn,i);
  //_DEBUGHERE_("deriv %g %g %g",pars[i],ans,errn);
  return ans;
}
#undef CON
#undef CON2
#undef BIG
#undef NTAB
#undef SAFE

double d2_func_dx_dy(lgpdf_c func, int idim, int jdim, double hh, double hx, double hy, double *x, double x0, double y0, void *self, error **err)
{
   double val[4], a;

   if (idim!=jdim) {
      /* Off-diagonal */

      x[idim] = x0+hh*hx;
      x[jdim] = y0+hh*hy;
      val[0]  = (*func)(self,x,err);
      forwardError(*err,__LINE__,0);
  
      x[idim] = x0+hh*hx;
      x[jdim] = y0-hh*hy;
      val[1]  = (*func)(self,x,err);
      forwardError(*err,__LINE__,0);

      x[idim] = x0-hh*hx;
      x[jdim] = y0+hh*hy;
      val[2]  = (*func)(self,x,err);
      forwardError(*err,__LINE__,0);

      x[idim] = x0-hh*hx;
      x[jdim] = y0-hh*hy;
      val[3]  = (*func)(self,x,err);
      forwardError(*err,__LINE__,0);
  
      a = (val[0]-val[1]-val[2]+val[3])/(4.0*hh*hh*hx*hy);
   } else {

      /* Diagonal */

      x[jdim] = y0;

      x[idim] = x0+2.0*hh*hx;
      val[0]  = (*func)(self, x, err);
      forwardError(*err, __LINE__, 0.0);

      x[idim] = x0;
      val[1]  = (*func)(self, x, err);
      forwardError(*err, __LINE__, 0.0);

      x[idim] = x0-2.0*hh*hx;
      val[2]  = (*func)(self, x, err);
      forwardError(*err, __LINE__, 0.0);

      a = (val[0] - 2.0*val[1] + val[2])/(4.0*hh*hh*hx*hx);
   }

   return a;
}

#define CON 1.4
#define CON2 (CON*CON)
#define BIG 1.0e30
#define NTAB 10
#define SAFE 2.0
double nd_dfridr2(lgpdf_c func, int idim, int jdim, double *pars, double hx, double hy, void* self,
		  double *errn,  error **err)
{
  int i,j;
  double errt,fac,hh,*a,ans,*x,x0,y0;
  
  ans = 1e30; /* dummy initialization */
  testErrorRetVA(hx==0.0 || hy==0.0, math_wrongValue, "h has to be larger than zero (got %g %g)", *err, __LINE__, 0.0, hx,hy);


  x = pars;
  x0 = x[idim];
  y0 = x[jdim];

  a = malloc_err(sizeof(double)*NTAB*NTAB,err);
  forwardError(*err,__LINE__,0);

  hh = 1;

  a[0] = d2_func_dx_dy(func, idim, jdim, hh, hx, hy, x, x0, y0, self, err);
  forwardError(*err, __LINE__, 0.0);

  *errn=BIG;
  for (i=1;i<NTAB;i++) {
    hh /= CON;

    a[i] = d2_func_dx_dy(func, idim, jdim, hh, hx, hy, x, x0, y0, self, err);

    fac = CON2;
    for (j=1;j<=i;j++) {
      a[j*NTAB+i] = (a[(j-1)*NTAB+i]*fac - a[(j-1)*NTAB+i-1]) / (fac-1.0);
      fac = CON2*fac;
      
      errt = fmax(fabs(a[j*NTAB+i] - a[(j-1)*NTAB+i]),fabs(a[j*NTAB+i] - a[(j-1)*NTAB+i-1]));
      if (errt <= *errn) {
        *errn = errt;
        ans = a[j*NTAB+i];
      }
    }
    
    if (fabs(a[i*NTAB+i]-a[(i-1)*NTAB+i-1]) >= SAFE*(*errn)) 
      break;
  }
  pars[idim]=x0;
  pars[jdim]=y0;
  free(a);

  _DEBUGHERE_("DF2(%d,%d) = %15.10f +- %g (i=%d, hh=%g)",idim,jdim,ans,*errn,i,hh);

  return ans;
}
#undef CON
#undef CON2
#undef BIG
#undef NTAB
#undef SAFE

/* ============================================================ *
 * Inverts the NxN matrix C.					*
 * ============================================================ */
double sm2_inverse(double *C, int N, error **err) {
   int i, j;
   double *col, det, *Cinv;
   int *p;

   Cinv = (double*)malloc_err((N*N)*sizeof(double), err);   forwardError(*err, __LINE__, 0);
   col = (double*)malloc_err((N)*sizeof(double), err);      forwardError(*err, __LINE__, 0);
   p = (int*)malloc_err((N)*sizeof(int), err);              forwardError(*err, __LINE__, 0);
   sm2_ludcmp(C, N, p, &det, err);                          forwardError(*err, __LINE__, 0);
   for (j=0; j<N; j++) {
      det *= C[j+j*N];
   }
   for (j=0; j<N; j++) {
      for (i=0; i<N; i++) {
         col[i] = 0.0;
      }
      col[j] = 1.0;
      sm2_lubksb(C, N, p, col);
      for (i=0; i<N; i++) {
         Cinv[i*N+j] = col[i];
      }
   }

   memcpy(C,Cinv,sizeof(double)*N*N);

   free(Cinv);
   free(col);
   free(p);
   return det;
}

/* ============================================================ *
 * Returns a copy of the NxN matrix C.				*
 * ============================================================ */
double *copy_matrix(double *C, int N, error **err)
{
   double *copy;

   copy = malloc_err(sizeof(double)*N*N, err);
   forwardError(*err, __LINE__, NULL);
   memcpy(copy, C, sizeof(double)*N*N);

   return copy;
}

/* ============================================================ *
 * Returns A*B where A and B are NxN matrices.			*
 * ============================================================ */
double *multiply_matrix(const double *A, const double *B, int N, error **err)
{
   int i, j, k;
   double *C, sum;

   C = malloc_err(sizeof(double)*N*N, err);
   forwardError(*err, __LINE__, NULL);

   for (i=0; i<N; i++) {
      for (j=0; j<N; j++) {
	 for (k=0,sum=0.0; k<N; k++) {
	    sum += A[N*i+k] * B[N*k+j];
	 }
	 C[N*i+j] = sum;
      }
   }
   return C;
}

/* ============================================================ *
 * Writes the NxN matrix A to file "name".			*
 * ============================================================ */
void write_matrix(const double *A, int N, const char name[], error **err)
{
   FILE *F;
   int i, j;

   F = fopen_err(name, "w", err);       forwardError(*err, __LINE__,);
   for (i=0; i<N; i++) {
      for (j=0; j<N; j++) {
	 fprintf(F, "%.10e ", A[N*i+j]);
      }
      fprintf(F, "\n");
   }
   fclose(F);
}

double *corr_coeff(double *C, int N, error **err)
{
   double *r;
   int j, k;

   r = malloc_err(sizeof(double)*N*N, err);
   forwardError(*err, __LINE__, NULL);

   for (j=0; j<N; j++) {
      for (k=0; k<j; k++) {
         testErrorRet(C[j*N+j]<=0 || C[k*N+k]<=0, math_negative, "Diagonal not positive", *err, __LINE__, NULL);
         r[j*N+k] = C[j*N+k]/sqrt(C[j*N+j]*C[k*N+k]);
      }
      r[j*N+j] = 1.0;
   }

   return r;
}

#define TINY 1.0e-20;
void sm2_ludcmp(double *a, int n, int *indx, double *d, error **err)
{
   int i,imax=-1,j,k;
   double big,dum,sum,temp;
   double *vv;

   vv = malloc_err(sizeof(double)*n, err);
   forwardError(*err, __LINE__,);
   *d=1.0;
   for (i=0;i<n;i++) {
      big=0.0;
      for (j=0;j<n;j++)
        if ((temp=fabs(a[i*n+j])) > big) big=temp;
      testErrorRet(big==0.0, math_singularValue, "Singular matrix", *err,__LINE__,);
      vv[i]=1.0/big;
   }
   for (j=0;j<n;j++) {
      for (i=0;i<j;i++) {
         sum=a[i*n+j];
         for (k=0;k<i;k++)
           sum -= a[i*n+k]*a[k*n+j];
         a[i*n+j]=sum;
      }
      big=0.0;
      for (i=j;i<n;i++) {
         sum=a[i*n+j];
         for (k=0;k<j;k++)
           sum -= a[i*n+k]*a[k*n+j];
         a[i*n+j]=sum;
         if ( (dum=vv[i]*fabs(sum)) >= big) {
            big=dum;
            imax=i;
         }
      }
      if (j != imax) {
         for (k=0;k<n;k++) {
            dum=a[imax*n+k];
            a[imax*n+k]=a[j*n+k];
            a[j*n+k]=dum;
         }
         *d = -(*d);
         vv[imax]=vv[j];
      }
      indx[j]=imax;
      if (a[j*n+j] == 0.0)
        a[j*n+j]=TINY;
      if (j != n-1) {
         dum=1.0/(a[j*n+j]);
         for (i=j+1;i<n;i++)
           a[i*n+j] *= dum;
      }
   }
   free(vv);
}
#undef TINY

void sm2_lubksb(double *a, int n, int *indx, double b[])
{
   int i,ii=-1,ip,j;
   double sum;

   for (i=0;i<n;i++) {
      ip=indx[i];
      sum=b[ip];
      b[ip]=b[i];
      if (ii>-1)
        for (j=ii;j<i;j++)
          sum -= a[i*n+j]*b[j];
      else if (sum)
        ii=i;
      b[i]=sum;
   }
   for (i=n-1;i>=0;i--) {
      sum=b[i];
      for (j=i+1;j<n;j++)
        sum -= a[i*n+j]*b[j];
      b[i]=sum/a[i*n+i];
   }
}

/* Cf. mvdens.c:scalar_product. L has to be a lower triangle matrix */
double scalar_product_mk(const double *mean, const double *L, const double *x, int ndim, error **err)
{
  double *x_tmp;
  double logl;
  int i, j;
  
  /* Copy vector */
  x_tmp = calloc_err(ndim, sizeof(double), err);
  forwardError(*err, __LINE__, 0);
  
  for (i=0; i<ndim; i++) {
     for (j=i; j<ndim; j++) {
	x_tmp[i] += (mean[i]-x[i])*L[i*ndim+j];
     }
  }

  /* Compute squared norm */
  for (i=0,logl=0.0; i<ndim; i++) {
    logl += x_tmp[i]*x_tmp[i];
  }
  
  free(x_tmp);

  return logl;
}

/* ============================================================ *
 * Multiplies the n-dimensional vector (or nxn-dim. matrix)     *
 * with the scalar c.						*
 * ============================================================ */
void multiply_all(double *m, int n, double c)
{
   int i;

   for (i=0; i<n; i++) {
      m[i] *= c;
   }
}

#ifndef __APPLE__
double fmin(double maxarg1,double maxarg2)
{
   return (maxarg1) < (maxarg2) ? maxarg1 : maxarg2;
}

double fmax(double maxarg1,double maxarg2)
{
   return (maxarg1) > (maxarg2) ? maxarg1 : maxarg2;
}
#endif

#ifndef __dsqr__
#define __dsqr__
double dsqr(double a)
{
   return a*a;
}
#endif

double DCUB(double a)
{
   return a*a*a;
}

/* ============================================================ *
 * Special functions.						*
 * ============================================================ */

double sinc(double x)
{
   if (x<0.001 && x>-0.001) return 1. - x*x/6.;
   else return sin(x)/x;
}

double gammln(double xx)
{
	double x,y,tmp,ser;
	const double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return (double)(-tmp+log(2.5066282746310005*ser/x));
}


double bessj0(double x)
{
	double ax,z;
	double xx,y,ans,ans1,ans2;

	if ((ax=fabs(x)) < 8.0) {
		y=x*x;
		ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7
			+y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
		ans2=57568490411.0+y*(1029532985.0+y*(9494680.718
			+y*(59272.64853+y*(267.8532712+y*1.0))));
		ans=ans1/ans2;
	} else {
		z=8.0/ax;
		y=z*z;
		xx=ax-0.785398164;
		ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
			+y*(-0.2073370639e-5+y*0.2093887211e-6)));
		ans2 = -0.1562499995e-1+y*(0.1430488765e-3
			+y*(-0.6911147651e-5+y*(0.7621095161e-6
			-y*0.934935152e-7)));
		ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
	}
	return ans;
}

double **sm2_matrix(long nrl, long nrh, long ncl, long nch, error **err)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{  
   long i, nrow=nrh-nrl+1, ncol=nch-ncl+1; 
   double **m;

   /* allocate pointers to rows */
   m=(double **) malloc_err((size_t)((nrow+NR_END)*sizeof(double*)), err);
   forwardError(*err, __LINE__, NULL);

   m += NR_END;
   m -= nrl;

   /* allocate rows and set pointers to them */
   m[nrl]=(double *) malloc_err((size_t)((nrow*ncol+NR_END)*sizeof(double)), err);
   forwardError(*err, __LINE__, NULL);

   m[nrl] += NR_END;
   m[nrl] -= ncl;

   for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

   /* return pointer to array of pointers to rows */
   return m;
}  

void sm2_free_matrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by sm2_matrix() */
{  
   free((FREE_ARG) (m[nrl]+ncl-NR_END));
   free((FREE_ARG) (m+nrl-NR_END));
}

double *sm2_vector(long nl, long nh, error **err)
/* allocate a double vector with subscript range v[nl..nh] */
{  
   double *v;

   v=(double *)malloc_err((size_t) ((nh-nl+1+NR_END)*sizeof(double)), err);
   forwardError(*err, __LINE__, NULL);

   return v-nl+NR_END;
}

void sm2_free_vector(double *v, long nl, long nh)
/* free a double vector allocated with vector() */
{
   free((FREE_ARG) (v+nl-NR_END));
}

void sm2_polint(double xa[], double ya[], int n, double x, double *y, double *dy, error **err)
{
   int i,m,ns=1;
   double den,dif,dift,ho,hp,w;
   double *c,*d;

   dif=fabs(x-xa[1]);
   c=sm2_vector(1,n,err);
   forwardError(*err,__LINE__,);
   d=sm2_vector(1,n,err);
   forwardError(*err,__LINE__,);
   for (i=1;i<=n;i++) {
      if ( (dift=fabs(x-xa[i])) < dif) {
         ns=i;
         dif=dift;
      }
      c[i]=ya[i];
      d[i]=ya[i];
   }
   *y=ya[ns--];
   for (m=1;m<n;m++) {
      for (i=1;i<=n-m;i++) {
         ho=xa[i]-x;
         hp=xa[i+m]-x;
         w=c[i+1]-d[i];
         den=ho-hp;
         testErrorRet(den==0.0, math_wrongValue, "den cannot be 0", *err, __LINE__,);
         den=w/den;
         d[i]=hp*den;
         c[i]=ho*den;
      }
      *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
   }
   sm2_free_vector(d,1,n);
   sm2_free_vector(c,1,n);
}

void sm2_spline(double x[], double y[], int n, double yp1, double ypn, double y2[], error **err)
{
   int i,k;
   double p,qn,sig,un,*u;

   u=sm2_vector(1,n-1,err);
   forwardError(*err,__LINE__,);
   if (yp1 > 0.99e30)
     y2[1]=u[1]=0.0;
   else {
      y2[1] = -0.5;
      u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
   }
   for (i=2;i<=n-1;i++) {
      sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
      p=sig*y2[i-1]+2.0;
      y2[i]=(sig-1.0)/p;
      u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
      u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
   }
   if (ypn > 0.99e30)
     qn=un=0.0;
   else {
      qn=0.5;
      un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
   }
   y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
   for (k=n-1;k>=1;k--) {
     y2[k]=y2[k]*y2[k+1]+u[k];
   }
   sm2_free_vector(u,1,n-1);
}

void sm2_splint(double xa[], double ya[], double y2a[], int n, double x, double *y, error **err)
{
   int klo,khi,k;
   double h,b,a;

   klo=1;
   khi=n;
   while (khi-klo > 1) {
      k=(khi+klo) >> 1;
      if (xa[k] > x) khi=k;
      else klo=k;
   }
   h=xa[khi]-xa[klo];
   testErrorRet(h==0.0, math_wrongValue, "h cannot be 0", *err, __LINE__,);
   a=(xa[khi]-x)/h;
   b=(x-xa[klo])/h;
   *y=a*ya[klo]+b*ya[khi]+
     ((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

void sm2_splin2(double x1a[], double x2a[], double **ya, double **y2a, int m, int n,
                double x1, double x2, double *y, error **err)
{
   int j;
   double *ytmp,*yytmp;

   ytmp  = sm2_vector(1,m,err); forwardError(*err, __LINE__,);
   yytmp = sm2_vector(1,m,err); forwardError(*err, __LINE__,);
   for (j=1;j<=m;j++) {
      sm2_splint(x2a,ya[j],y2a[j],n,x2,&yytmp[j],err);
      forwardError(*err, __LINE__,);
   }
   sm2_spline(x1a,yytmp,m,1.0e30,1.0e30,ytmp,err);  forwardError(*err, __LINE__,);
   sm2_splint(x1a,yytmp,ytmp,m,x1,y,err);           forwardError(*err, __LINE__,);
   sm2_free_vector(yytmp,1,m);
   sm2_free_vector(ytmp,1,m);
}

splineTable *init_splineTable(int n, error **err)
{
   splineTable *self;

   self     = malloc_err(sizeof(splineTable), err);            forwardError(*err, __LINE__, NULL);
   self->x  = (double*)malloc_err(sizeof(double)*(n+1), err);  forwardError(*err, __LINE__, NULL);
   self->y  = (double*)malloc_err(sizeof(double)*(n+1), err);  forwardError(*err, __LINE__, NULL);
   self->y2 = (double*)malloc_err(sizeof(double)*(n+1), err);  forwardError(*err, __LINE__, NULL);
   self->n  = n;

   return self;
}

splineTable *copy_splineTable(const splineTable *self, error **err)
{
   int n;
   splineTable *res;

   if (self==NULL) return NULL;

   n   = self->n;
   res = init_splineTable(n, err);   forwardError(*err, __LINE__, NULL);

   res->yp1 = self->yp1;
   res->ypn = self->ypn;
   memcpy(res->x, self->x, n*sizeof(double));
   testErrorRet(res->x==NULL, math_alloc, "Could not allocate memory", *err, __LINE__, NULL);
   memcpy(res->y, self->y, n*sizeof(double));
   testErrorRet(res->y==NULL, math_alloc, "Could not allocate memory", *err, __LINE__, NULL);
   memcpy(res->y2, self->y2, n*sizeof(double));
   testErrorRet(res->y2==NULL, math_alloc, "Could not allocate memory", *err, __LINE__, NULL);

   return res;
}

void del_splineTable(splineTable** self)
{
   if (*self==NULL) return;
   free((*self)->x);
   free((*self)->y);
   free((*self)->y2);
   free(*self);
   *self=NULL;
}

void sm2_hunt(double xx[], unsigned long n, double x, unsigned long *jlo)
{
   unsigned long jm,jhi,inc;
   int ascnd;

   ascnd=(xx[n] >= xx[1]);
   if (*jlo <= 0 || *jlo > n) {
      *jlo=0;
      jhi=n+1;
   } else {
      inc=1;
      if ((x >= xx[*jlo]) == ascnd) {
         if (*jlo == n) return;
         jhi=(*jlo)+1;
         while ((x >= xx[jhi]) == ascnd) {
            *jlo=jhi;
            inc += inc;
            jhi=(*jlo)+inc;
            if (jhi > n) {
               jhi=n+1;
               break;
            }
         }
      } else {
         if (*jlo == 1) {
            *jlo=0;
            return;
         }
         jhi=(*jlo)--;
         while ((x < xx[*jlo]) == ascnd) {
            jhi=(*jlo);
            inc <<= 1;
            if (inc >= jhi) {
               *jlo=0;
               break;
            }
            else *jlo=jhi-inc;
         }
      }
   }
   while (jhi-(*jlo) != 1) {
      jm=(jhi+(*jlo)) >> 1;
      if ((x >= xx[jm]) == ascnd)
        *jlo=jm;
      else
        jhi=jm;
   }
   if (x == xx[n]) *jlo=n-1;
   if (x == xx[1]) *jlo=1;
}

#define FUNC(x,y,err) ((*func)(x,y,err))
double sm2_trapzdberg(funcwithpars func, void *intpar,
                      double a, double b, int n, double *s, error **err)
{
   double x,tnm,sum,del;
   int it,j;
   double res,b1,b2;

   if (n == 1) {
      b1 = (*func)(a,intpar,err);
      forwardError(*err,__LINE__,0);
      b2 = (*func)(b,intpar,err);
      forwardError(*err,__LINE__,0);

      res = (*s=0.5*(b-a)*(b1+b2));
      return res;
   } else {
      for (it=1,j=1;j<n-1;j++) it <<= 1;
      tnm=it;
      del=(b-a)/tnm;
      x=a+0.5*del;
      for (sum=0.0,j=1;j<=it;j++,x+=del) {
         sum += FUNC(x,intpar,err);
         forwardError(*err,__LINE__,0);
      }
      *s = 0.5*(*s+(b-a)*sum/tnm);
      return *s;
   }
}

/* ============================================================ *
 * Romberg Integration, see NR p. 140                           *
 * ============================================================ */
   
#define JMAX 35
#define JMAXP (JMAX+1)
#define K 5
double sm2_qrombergo(funcwithpars  func, void *intpar, 
                     double a, double b, double (*choose)(double(*)(double,
                     void *, error **), void *, double, double, int, double *, error**),
                     double EPS, error **err)
{
   int j;
   double ss,dss,h[JMAXP+1],s[JMAXP], schoose;

   h[1]=1.0;
   for (j=1;j<=JMAX;j++) {
      s[j]=(*choose)(func,intpar,a,b,j,&schoose,err);
      forwardError(*err,__LINE__,0);
      if (j >= K) {
         sm2_polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss,err);
         forwardError(*err,__LINE__,0);
         testErrorRet(!finite(dss), math_infnan, "Inf or nan encountered", *err, __LINE__, 0.0);
         if (fabs(dss) <= EPS*fabs(ss)) return ss;
      }
      h[j+1]=h[j]/9.0;
   }
   *err=addError(math_tooManySteps,"Too many steps",*err,__LINE__);
   return 0.0;
}
#undef JMAX
#undef JMAXP

#define JMAX 35
#define JMAXP (JMAX+1)
#define K 5
double sm2_qromberg(funcwithpars func, double *intpar,
                    double a, double b, double EPS, error **err)
{
   double ss,dss;
   double s[JMAXP],h[JMAXP+1], strap;
   int j;

   h[1]=1.0;
   for (j=1;j<=JMAX;j++) {
      s[j] = sm2_trapzdberg(func,intpar,a,b,j,&strap,err);
      forwardError(*err,__LINE__,0);

      if (j >= K) {
         sm2_polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss,err);
         forwardError(*err,__LINE__,0);
         testErrorRet(!finite(dss), math_infnan, "Inf or nan encountered", *err, __LINE__, 0.0);
         if (fabs(dss) <= EPS*fabs(ss)) return ss;
      }
      h[j+1]=0.25*h[j];
   }
   *err=addError(math_tooManySteps,"Too many steps",*err,__LINE__);
   return 0.0;
}
#undef JMAX
#undef JMAXP
#undef K

double sm2_midpntberg(funcwithpars  func, void *intpar,
                      double a, double b, int n, double *s, error **err)
{
   double x,tnm,sum,del,ddel;
   int it,j;
   if (n == 1) {
      *s=(b-a)*FUNC(0.5*(a+b), intpar,err);
      forwardError(*err,__LINE__,0);
      return (*s);
   } else {
      for(it=1,j=1;j<n-1;j++) it *= 3;
      tnm=it;
      del=(b-a)/(3.0*tnm);
      ddel=del+del;
      x=a+0.5*del;
      sum=0.0;
      for (j=1;j<=it;j++) {
         sum += FUNC(x, intpar,err);
         forwardError(*err,__LINE__,0);
         x += ddel;
         sum += FUNC(x, intpar,err);
         forwardError(*err,__LINE__,0);
         x += del;
      }
      *s = (*s+(b-a)*sum/tnm)/3.0;
      return *s;
   }
}
#undef FUNC

/* Initialises an rkdrive structure, used e.g. in odeint. *
 * Necessary (for bsstep) to set initial values are only  *
 * first and epsold. The other variables are set for      *
 * consistency.						  */
void init_rkdrive_var(rkdrive_var_t *rkvar)
{
   int i, j;
   rkvar->first = 1;
   rkvar->kmax = rkvar->kopt = 0;
   rkvar->epsold = -1.0;

   for (i=0; i<=IMAXX; i++) {
      rkvar->a[i] = 0.0;
   }
   for (i=0; i<=KMAXX; i++) {
      for (j=0; j<KMAXX; j++) {
	 rkvar->alf[i][j] = 0.0;
      }
   }
}


/* IMAXX, KMAXX now defined in maths.h */
#define SAFE1 0.25
#define SAFE2 0.7
#define REDMAX 1.0e-5
#define REDMIN 0.7
#define TINY 1.0e-30
#define SCALMX 0.1
void bsstep(double y[], double dydx[], int nv, double *xx, double htry, double eps,
            double yscal[], double *hdid, double *hnext, void* extra,
            derivative derivs, rkdrive_var_t *rkvar, error **erro)
{
   int i,iq,k,kk,km=0;
   double xnew=0.0;
   double eps1,errmax,fact,h,red=0.0,scale=0.0,work,wrkmin,xest;
   double *err,*yerr,*ysav,*yseq;
   const int nseq[IMAXX+1]={0,2,4,6,8,10,12,14,16,18};
   int reduct,exitflag=0;
   double **d,*x;

   d    = sm2_matrix(1,nv,1,KMAXX,erro); forwardError(*erro,__LINE__,);
   err  = sm2_vector(1,KMAXX,erro);      forwardError(*erro,__LINE__,);
   x    = sm2_vector(1,KMAXX,erro);      forwardError(*erro,__LINE__,);
   yerr = sm2_vector(1,nv,erro);         forwardError(*erro,__LINE__,);
   ysav = sm2_vector(1,nv,erro);         forwardError(*erro,__LINE__,);
   yseq = sm2_vector(1,nv,erro);         forwardError(*erro,__LINE__,);

   if (eps != rkvar->epsold) {
      *hnext = xnew = -1.0e29;
      eps1   = SAFE1 * eps;
      rkvar->a[1] = nseq[1] + 1;
      for (k=1; k<=KMAXX; k++) rkvar->a[k+1] = rkvar->a[k] + nseq[k+1];
      for (iq=2; iq<=KMAXX; iq++) {
         for (k=1; k<iq; k++)
           rkvar->alf[k][iq] = pow(eps1, (rkvar->a[k+1] - rkvar->a[iq+1])/((rkvar->a[iq+1] - rkvar->a[1]+1.0) * (2*k+1)));
      }
      rkvar->epsold = eps;
      for (rkvar->kopt=2; rkvar->kopt<KMAXX; rkvar->kopt++)
        if (rkvar->a[rkvar->kopt+1] > rkvar->a[rkvar->kopt] * rkvar->alf[rkvar->kopt-1][rkvar->kopt]) break;
      rkvar->kmax = rkvar->kopt;
   }
   h=htry;
   for (i=1;i<=nv;i++) ysav[i]=y[i];
   if (*xx != xnew || h != (*hnext)) {
      rkvar->first = 1;
      rkvar->kopt  = rkvar->kmax;
   }
   reduct=0;
   for (;;) {
      for (k=1; k<=rkvar->kmax; k++) {
         xnew=(*xx)+h;
         testErrorRet(xnew==(*xx), math_underflow, "Step size underflow", *erro, __LINE__,);
         mmid(ysav,dydx,nv,*xx,h,nseq[k],yseq,extra,derivs,erro);
         forwardError(*erro,__LINE__,);
         xest=dsqr(h/nseq[k]);
         pzextr(k,xest,yseq,y,yerr,nv,x,d,erro);
         forwardError(*erro,__LINE__,);
         if (k != 1) {
            errmax=TINY;
            for (i=1;i<=nv;i++) errmax=fmax(errmax,fabs(yerr[i]/yscal[i]));
            errmax /= eps;
            km=k-1;
            err[km]=pow(errmax/SAFE1,1.0/(2*km+1));
         }
         if (k != 1 && (k >= rkvar->kopt-1 || rkvar->first)) {
            if (errmax < 1.0) {
               exitflag=1;
               break;
            }
            if (k == rkvar->kmax || k == rkvar->kopt+1) {
               red=SAFE2/err[km];
               break;
            }
            else if (k == rkvar->kopt && rkvar->alf[rkvar->kopt-1][rkvar->kopt] < err[km]) {
               red=1.0/err[km];
               break;
            }
            else if (rkvar->kopt == rkvar->kmax && rkvar->alf[km][rkvar->kmax-1] < err[km]) {
               red = rkvar->alf[km][rkvar->kmax-1]*SAFE2/err[km];
               break;
            }
            else if (rkvar->alf[km][rkvar->kopt] < err[km]) {
               red = rkvar->alf[km][rkvar->kopt-1]/err[km];
               break;
            }
         }
      }
      if (exitflag) break;
      red=fmin(red,REDMIN);
      red=fmax(red,REDMAX);
      h *= red;
      reduct=1;
   }
   *xx=xnew;
   *hdid=h;
   rkvar->first = 0;
   wrkmin=1.0e35;
   for (kk=1; kk<=km; kk++) {
      fact = fmax(err[kk],SCALMX);
      work = fact * rkvar->a[kk+1];
      if (work < wrkmin) {
         scale=fact;
         wrkmin=work;
         rkvar->kopt = kk+1;
      }
   }
   *hnext=h/scale;
   if (rkvar->kopt >= k && rkvar->kopt != rkvar->kmax && !reduct) {
      fact=fmax(scale/rkvar->alf[rkvar->kopt-1][rkvar->kopt],SCALMX);
      if (rkvar->a[rkvar->kopt+1] * fact <= wrkmin) {
         *hnext=h/fact;
         rkvar->kopt++;
      }
   }
   sm2_free_vector(yseq,1,nv);
   sm2_free_vector(ysav,1,nv);
   sm2_free_vector(yerr,1,nv);
   sm2_free_vector(x,1,KMAXX);
   sm2_free_vector(err,1,KMAXX);
   sm2_free_matrix(d,1,nv,1,KMAXX);
}
//#undef KMAXX
//#undef IMAXX
#undef SAFE1
#undef SAFE2
#undef REDMAX
#undef REDMIN
#undef TINY
#undef SCALMX

void mmid(double y[], double dydx[], int nvar, double xs, double htot, int nstep,
          double yout[], void* extra,
          derivative derivs, error **err)
{
   int n,i;
   double x,swap,h2,h,*ym,*yn;

   ym=sm2_vector(1,nvar,err);
   forwardError(*err,__LINE__,);
   yn=sm2_vector(1,nvar,err);
   forwardError(*err,__LINE__,);
   h=htot/nstep;
   for (i=1;i<=nvar;i++) {
      ym[i]=y[i];
      yn[i]=y[i]+h*dydx[i];
   }
   x=xs+h;
   (*derivs)(x,yn,yout,extra,err);
   forwardError(*err,__LINE__,);
   h2=2.0*h;
   for (n=2;n<=nstep;n++) {
      for (i=1;i<=nvar;i++) {
         swap=ym[i]+h2*yout[i];
         ym[i]=yn[i];
         yn[i]=swap;
      }
      x += h;
      (*derivs)(x,yn,yout,extra,err);
      forwardError(*err,__LINE__,);
   }
   for (i=1;i<=nvar;i++)
     yout[i]=0.5*(ym[i]+yn[i]+h*yout[i]);
   sm2_free_vector(yn,1,nvar);
   sm2_free_vector(ym,1,nvar);
}

void pzextr(int iest, double xest, double yest[], double yz[], double dy[], int nv,double *x, double **d, error **err)
{
   int k1,j;
   double q,f2,f1,delta,*c;

   c=sm2_vector(1,nv,err);
   forwardError(*err,__LINE__,);
   x[iest]=xest;
   for (j=1;j<=nv;j++) dy[j]=yz[j]=yest[j];
   if (iest == 1) {
      for (j=1;j<=nv;j++) d[j][1]=yest[j];
   } else {
      for (j=1;j<=nv;j++) c[j]=yest[j];
      for (k1=1;k1<iest;k1++) {
         delta=1.0/(x[iest-k1]-xest);
         f1=xest*delta;
         f2=x[iest-k1]*delta;
         for (j=1;j<=nv;j++) {
            q=d[j][k1];
            d[j][k1]=dy[j];
            delta=c[j]-q;
            dy[j]=f1*delta;
            c[j]=f2*delta;
            yz[j] += dy[j];
         }
      }
      for (j=1;j<=nv;j++) d[j][iest]=dy[j];
   }
   sm2_free_vector(c,1,nv);
}

#define MAXSTP 10000
#define TINY 1.0e-30
void odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1,
            void *extra, double hmin, int *nok, int *nbad,
            derivative derivs, rkdrive rkqs, rkdrive_var_t *rkvar, error **err)
{
   int nstp,i;
   double x,hnext,hdid,h;
   double *yscal,*y,*dydx;

   yscal=sm2_vector(1,nvar,err);   forwardError(*err,__LINE__,);
   y=sm2_vector(1,nvar,err);       forwardError(*err,__LINE__,);
   dydx=sm2_vector(1,nvar,err);    forwardError(*err,__LINE__,);
   x=x1;
   h=SIGN(h1,x2-x1);
   *nok = (*nbad) =  0;

   for (i=1;i<=nvar;i++) y[i]=ystart[i];

   for (nstp=1;nstp<=MAXSTP;nstp++) {

      (*derivs)(x,y,dydx,extra,err);
      forwardError(*err,__LINE__,);
      for (i=1;i<=nvar;i++)
        yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
      if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
      (*rkqs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,extra,derivs,rkvar,err);
      forwardError(*err,__LINE__,);
      if (hdid == h) ++(*nok); else ++(*nbad);
      if ((x-x2)*(x2-x1) >= 0.0) {

         for (i=1;i<=nvar;i++) ystart[i]=y[i];
         sm2_free_vector(dydx,1,nvar);
         sm2_free_vector(y,1,nvar);
         sm2_free_vector(yscal,1,nvar);
         return;

      }
      testErrorRet(fabs(hnext)<=hmin, math_underflow, "Step size too small", *err, __LINE__,);
      h=hnext;

   }

   *err = addError(math_tooManySteps, "Too many steps", *err, __LINE__);
   return;

}
#undef MAXSTP
#undef TINY

void sm2_Cmul(my_complex a, my_complex b, my_complex *c)
{
   (*c)[0] = a[0]*b[0] - a[1]*b[1];
   (*c)[1] = a[1]*b[0] + a[0]*b[1];
}

void sm2_Cdiv(my_complex a, my_complex b, my_complex *c)
{
   double r,den;
   if (fabs(b[0]) >= fabs(b[1])) {
      r=b[1]/b[0];
      den=b[0]+r*b[1];
      (*c)[0]=(a[0]+r*a[1])/den;
      (*c)[1]=(a[1]-r*a[0])/den;
   } else {
      r=b[0]/b[1];
      den=b[1]+r*b[0];
      (*c)[0]=(a[0]*r+a[1])/den;
      (*c)[1]=(a[1]*r-a[0])/den;
   }
}

/* ============================================================ *
 * Complex gamma function.					*
 * ============================================================ */
void cdgamma(my_complex x, my_complex *res)
{
   double          xr, xi, wr, wi, ur, ui, vr, vi, yr, yi, t;

   xr = (double) x[0];
   xi = (double) x[1];

   if (xr<0) {
      wr = 1 - xr;
      wi = -xi;
   } else {
      wr = xr;
      wi = xi;
   }

   ur = wr + 6.00009857740312429;
   vr = ur * (wr + 4.99999857982434025) - wi * wi;
   vi = wi * (wr + 4.99999857982434025) + ur * wi;
   yr = ur * 13.2280130755055088 + vr * 66.2756400966213521 + 
     0.293729529320536228;
   yi = wi * 13.2280130755055088 + vi * 66.2756400966213521;
   ur = vr * (wr + 4.00000003016801681) - vi * wi;
   ui = vi * (wr + 4.00000003016801681) + vr * wi;
   vr = ur * (wr + 2.99999999944915534) - ui * wi;
   vi = ui * (wr + 2.99999999944915534) + ur * wi;
   yr += ur * 91.1395751189899762 + vr * 47.3821439163096063;
   yi += ui * 91.1395751189899762 + vi * 47.3821439163096063;
   ur = vr * (wr + 2.00000000000603851) - vi * wi;
   ui = vi * (wr + 2.00000000000603851) + vr * wi;
   vr = ur * (wr + 0.999999999999975753) - ui * wi;
   vi = ui * (wr + 0.999999999999975753) + ur * wi;
   yr += ur * 10.5400280458730808 + vr;
   yi += ui * 10.5400280458730808 + vi;
   ur = vr * wr - vi * wi;
   ui = vi * wr + vr * wi;
   t = ur * ur + ui * ui;
   vr = yr * ur + yi * ui + t * 0.0327673720261526849;
   vi = yi * ur - yr * ui;
   yr = wr + 7.31790632447016203;
   ur = log(yr * yr + wi * wi) * 0.5 - 1;
   ui = atan2(wi, yr);
   yr = exp(ur * (wr - 0.5) - ui * wi - 3.48064577727581257) / t;
   yi = ui * (wr - 0.5) + ur * wi;
   ur = yr * cos(yi);
   ui = yr * sin(yi);
   yr = ur * vr - ui * vi;
   yi = ui * vr + ur * vi;
   if (xr<0) {
      wr = xr * 3.14159265358979324;
      wi = exp(xi * 3.14159265358979324);
      vi = 1 / wi;
      ur = (vi + wi) * sin(wr);
      ui = (vi - wi) * cos(wr);
      vr = ur * yr + ui * yi;
      vi = ui * yr - ur * yi;
      ur = 6.2831853071795862 / (vr * vr + vi * vi);
      yr = ur * vr;
      yi = ur * vi;
   }

   (*res)[0]=yr; (*res)[1]=yi;
}


/* ============================================================ *
 * Interpolation functions.					*
 * ============================================================ */

interTable* init_interTable(int n, double a, double b, double dx, double lower, double upper, error **err)
{
   interTable* self;
  
   self = malloc_err(sizeof(interTable), err);                forwardError(*err, __LINE__, NULL);
   self->table = (double*) malloc_err(sizeof(double)*n, err); forwardError(*err, __LINE__, NULL);
   self->n = n;
   self->a = a;
   self->b = b;
   self->dx = dx;
   self->lower = lower;
   self->upper = upper;

   return self;
}

interTable* copy_interTable(interTable* self, error **err)
{
   interTable* res;

   if (self == NULL) return NULL;

   res = init_interTable(self->n, self->a, self->b, self->dx, self->lower, self->upper, err);
   forwardError(*err,__LINE__, NULL);
   memcpy(res->table, self->table, self->n*sizeof(double));
   testErrorRet(res->table==NULL, math_alloc, "memcpy failed", *err, __LINE__, NULL);

   return res;
}

void del_interTable(interTable** self)
{
   if (*self==NULL) return;

   free((*self)->table);
   free(*self);
   *self = NULL;
}

interTable** init_interTable_arr(int N, int n, double a, double b, double dx, double lower, double upper, error **err)
{
   interTable** self;
   int i;

   self = malloc_err(sizeof(interTable*)*N, err);  forwardError(*err, __LINE__, NULL);

   for (i=0; i<N; i++) {
      self[i] = init_interTable(n, a, b, dx, lower, upper, err);
      forwardError(*err, __LINE__, NULL);
   }

   return self;
}

interTable** copy_interTable_arr(interTable** self, int N, error **err)
{
   interTable** res;
   int n;

   if (self == NULL) return NULL;

   res = malloc_err(sizeof(interTable*)*N, err);
   forwardError(*err, __LINE__, NULL);
   for (n=0; n<N; n++) {
      res[n] = init_interTable(self[n]->n, self[n]->a, self[n]->b, self[n]->dx, self[n]->lower, self[n]->upper, err);
      forwardError(*err,__LINE__, NULL);
      memcpy(res[n]->table, self[n]->table, self[n]->n*sizeof(double));
      testErrorRet(res[n]->table==NULL, math_alloc, "memcpy failed", *err, __LINE__, NULL);
   }

   return res;
}

void del_interTable_arr(interTable*** self, int N)
{
   int i;

   if (*self==NULL) return;

   for (i=0; i<N; i++) {
      free((*self)[i]->table);
   }
   free(*self);
   *self = NULL;
}

interTable2D* init_interTable2D(int n1, double a1, double b1, double dx1, int n2, double a2, double b2, 
				double dx2, double lower, double upper, error **err)
{
   interTable2D* self;
  
   self = malloc_err(sizeof(interTable2D), err);  forwardError(*err, __LINE__, NULL);
   self->table = sm2_matrix(0,n1-1,0,n2-1,err);
   forwardError(*err,__LINE__,NULL);
   self->n1    = n1;
   self->a1    = a1;
   self->b1    = b1;
   self->dx1   = dx1;
   self->n2    = n2;
   self->a2    = a2;
   self->b2    = b2;
   self->dx2   = dx2;
   self->lower = lower;
   self->upper = upper;

   return self;
}

interTable2D* copy_interTable2D(interTable2D* self, error **err)
{
   interTable2D* res;
   int i;
  
   if (self == NULL)
     return NULL;
   res = init_interTable2D(self->n1,self->a1,self->b1,self->dx1,self->n2,self->a2,self->b2,self->dx2,
			   self->lower,self->upper,err);
   forwardError(*err,__LINE__,NULL);
   for(i=0;i<self->n1;i++)
     memcpy(res->table[i],self->table[i],self->n2*sizeof(double));
   return res;
}

void del_interTable2D(interTable2D** self)
{
   if (*self==NULL) return;
   sm2_free_matrix((*self)->table,0,(*self)->n1-1,0,(*self)->n2-1);
   free(*self);
   *self = NULL;
}

interTable2D** init_interTable2D_arr(int N, int n1, double a1, double b1, double dx1, int n2, double a2, double b2, 
				    double dx2, double lower, double upper, error **err)
{
   interTable2D **self;
   int i;

   self = malloc_err(sizeof(interTable2D*)*N, err); forwardError(*err, __LINE__, NULL);

   for (i=0; i<N; i++) {
      self[i] = init_interTable2D(n1, a1, b1, dx1, n2, a2, b2, dx2, lower, upper, err);
      forwardError(*err, __LINE__, NULL);
   }

   return self;
}

interTable2D** copy_interTable2D_arr(interTable2D **self, int N, error **err)
{
   interTable2D **res;
   int n, i;

   if (self==NULL) return NULL;

   res = malloc_err(sizeof(interTable2D*)*N, err);
   forwardError(*err, __LINE__, NULL);
   for (n=0; n<N; n++) {
      res[n] = init_interTable2D(self[n]->n1, self[n]->a1, self[n]->b1, self[n]->dx1,
				 self[n]->n2, self[n]->a2, self[n]->b2, self[n]->dx2,
				 self[n]->lower, self[n]->upper, err);
      forwardError(*err, __LINE__, NULL);
      for (i=0; i<self[n]->n1; i++) {
	 memcpy(res[n]->table[i], self[n]->table[i], self[n]->n2*sizeof(double));
	 testErrorRet(res[n]->table==NULL, math_alloc, "memcpy failed", *err, __LINE__, NULL);
      }      
   }

   return res;
}

void del_interTable2D_arr(interTable2D*** self, int N)
{
   int i;

   if (*self==NULL) return;

   for (i=0; i<N; i++) {
      sm2_free_matrix((*self)[i]->table, 0, (*self)[i]->n1-1, 0, (*self)[i]->n2-1);
   }
   free(*self);
   *self = NULL;
}

interTable2Dspline *init_interTable2Dspline(int m, int n, error **err)
{
   interTable2Dspline *self;
   int comp;

   self = (interTable2Dspline*)malloc_err(sizeof(interTable2Dspline), err);
   forwardError(*err, __LINE__, NULL);
   self->x  = sm2_vector(1, m, err);         forwardError(*err, __LINE__, NULL);
   self->y  = sm2_vector(1, n, err);         forwardError(*err, __LINE__, NULL);
   for (comp=0; comp<NCOMP; comp++) {
      self->z[comp]  = sm2_matrix(1, m, 1, n, err);    forwardError(*err, __LINE__, NULL);
      self->z2[comp] = sm2_matrix(1, m, 1, n, err);    forwardError(*err, __LINE__, NULL);
   }
   self->m  = m;
   self->n  = n;

   return self;
}

interTable2Dspline *copy_interTable2Dspline(interTable2Dspline *source, error **err)
{
   interTable2Dspline *new;
   int i, comp;

   if (source==NULL) return NULL;
   new = init_interTable2Dspline(source->m, source->n, err);
   forwardError(*err, __LINE__, NULL);
   memcpy(new->x+1, source->x+1, source->m*sizeof(double));
   memcpy(new->y+1, source->y+1, source->n*sizeof(double));
   for (comp=0; comp<NCOMP; comp++) {
      for (i=1; i<=source->m; i++) {
	 memcpy(new->z[comp][i]+1, source->z[comp][i]+1, source->n*sizeof(double));
	 memcpy(new->z2[comp][i]+1, source->z[comp][i]+1, source->n*sizeof(double));
      }
   }
   new->m = source->m;
   new->n = source->n;

   return new;
}

void del_interTable2Dspline(interTable2Dspline **self)
{
   int comp;
   if (*self==NULL) return;
   sm2_free_vector((*self)->x, 1, (*self)->m);
   sm2_free_vector((*self)->y, 1, (*self)->n);
   for (comp=0; comp<=NCOMP; comp++) {
      sm2_free_matrix((*self)->z[comp], 1, (*self)->m, 1, (*self)->n);
      sm2_free_matrix((*self)->z2[comp], 1, (*self)->m, 1, (*self)->n);
   }
   free(*self);
   *self = NULL;
}

interTable2Dneq *init_interTable2Dneq(int m, int n, error **err)
{
   interTable2Dneq *self;
   int comp;

   self = (interTable2Dneq*)malloc_err(sizeof(interTable2Dneq), err);
   forwardError(*err, __LINE__, NULL);
   self->x  = sm2_vector(1, m, err);          forwardError(*err, __LINE__, NULL);
   self->y  = sm2_vector(1, n, err);          forwardError(*err, __LINE__, NULL);
   for (comp=0; comp<NCOMP; comp++) {
      self->z[comp]  = sm2_matrix(1, m, 1, n, err);    forwardError(*err, __LINE__, NULL);
   }
   self->m  = m;
   self->n  = n;
   self->kx = self->ky = 1;   /* Default value */

   return self;
}

interTable2Dneq *copy_interTable2Dneq(interTable2Dneq *source, error **err)
{
   interTable2Dneq *new;
   int i, comp;

   if (source==NULL) return NULL;
   new = init_interTable2Dneq(source->m, source->n, err);
   forwardError(*err, __LINE__, NULL);
   memcpy(new->x+1, source->x+1, source->m*sizeof(double));
   memcpy(new->y+1, source->y+1, source->n*sizeof(double));
   for (comp=0; comp<NCOMP; comp++) {
      for (i=1; i<=source->m; i++) {
	 memcpy(new->z[comp][i]+1, source->z[comp][i]+1, source->n*sizeof(double));
      }
   }
   new->m  = source->m;
   new->n  = source->n;
   new->kx = source->kx;
   new->kx = source->ky;

   return new;
}

void del_interTable2Dneq(interTable2Dneq **self)
{
   int comp;
   if (*self==NULL) return;
   sm2_free_vector((*self)->x, 1, (*self)->m);
   sm2_free_vector((*self)->y, 1, (*self)->n);
   for (comp=0; comp<=NCOMP; comp++) {
      sm2_free_matrix((*self)->z[comp], 1, (*self)->m, 1, (*self)->n);
   }
   free(*self);
   *self = NULL;
}

double sm2_interpol2Dneq(interTable2Dneq *self, comp_t comp, double x0, double y0, error **err)
{
   double rx, ry, z0, logk1, logk2, logT1, logT2;

   sm2_hunt(self->x, self->m, x0, &self->kx);
   sm2_hunt(self->y, self->n, y0, &self->ky);

   ry = (y0 - self->y[self->ky])/(self->y[self->ky+1] - self->y[self->ky]);

   if (self->kx==self->m) {
      if (self->z[comp][self->m][self->ky]<EPSILON) return 0.0;

      logk1 = log(self->x[self->m-1]);
      logk2 = log(self->x[self->m]);
      logT1 = (1.0-ry)*log(self->z[comp][self->m-1][self->ky])
	+ ry*log(self->z[comp][self->m-1][self->ky+1]);
      logT2 = (1.0-ry)*log(self->z[comp][self->m][self->ky])
	+ ry*log(self->z[comp][self->m][self->ky+1]);

      z0 = ((log(x0)-logk2)/(logk2-logk1))*(logT2-logT1) + logT2;
      return exp(z0);
   } else if (self->kx==0) {
      if (self->z[comp][1][self->ky]<EPSILON) return 0.0;

      logk1 = log(self->x[1]);
      logk2 = log(self->x[2]);
      logT1 = (1.0-ry)*log(self->z[comp][1][self->ky])
	+ ry*log(self->z[comp][1][self->ky+1]);
      logT2 = (1.0-ry)*log(self->z[comp][2][self->ky])
	+ ry*log(self->z[comp][2][self->ky+1]);

      z0 = ((log(x0)-logk1)/(logk2-logk1))*(logT2-logT1) + logT1;
      return exp(z0);
   }


   testErrorRet(self->kx<1 || self->kx>self->m-1 || self->ky<1 || self->ky>self->n-1,
		math_interpoloutofrange, "out of range", *err, __LINE__, 0);

   rx = (x0 - self->x[self->kx])/(self->x[self->kx+1] - self->x[self->kx]);

   z0 = (1.0-rx)*(1.0-ry)*self->z[comp][self->kx][self->ky] + rx*(1.0-ry)*self->z[comp][self->kx+1][self->ky]
     + (1.0-rx)*ry*self->z[comp][self->kx][self->ky+1] + rx*ry*self->z[comp][self->kx+1][self->ky+1];

   return z0;
}

/* ============================================================ *
 * Interpolates f at the value x, where f is a double[n] array,	*
 * representing a function between a and b, stepwidth dx.	*
 * 'lower' and 'upper' are powers of a logarithmic power law	*
 * extrapolation. If no	extrapolation desired, set these to 0	*
 * ============================================================ */
double interpol_wr(interTable* self, double x, error **err) {
   double res;
   res = sm2_interpol(self->table, self->n, self->a, self->b, self->dx, x,
		      self->lower, self->upper, err);
   forwardError(*err, __LINE__, 0.0);
   return res;
}

double interpol2D(interTable2D* self, double x,double y, error **err) {
   double res;
   res = sm2_interpol2d(self->table, self->n1, self->a1, self->b1, self->dx1, x, self->n2, 
			self->a2, self->b2, self->dx2, y, self->lower, self->upper, err);
   forwardError(*err, __LINE__, 0.0);
   return res;
}

#define EPS 1.0-5
double sm2_interpol(double *f, int n, double a, double b, double dx, double x,
		    double lower, double upper, error **err)
{
   double r;
   int  i;
   char str[1024];

   if (x < a) {
      if (fabs(lower)<EPS) {
	 sprintf(str, "Query point for 1d interpolation x=%g smaller than minimum a=%g.\n"
		 "For linear extrapolation set lower!=0", x, a);
	 *err = addError(math_interpol2small, str, *err, __LINE__);
	 return 0.0;
      }
      return f[0] + lower*(x - a);
   }
   r = (x - a)/dx;
   i = (int)(floor(r));
   if (i+1 >= n) {
      if (fabs(upper)<EPS) {
	 if (i+1==n) {
	    return f[i];  /* constant extrapolation */
	 } else {
	    sprintf(str, "Query point for 1d interpolation x=%g larger than maximum b=%g.\n"
		    "For linear extrapolation set upper!=0", x, a);
	    *err = addError(math_interpol2big, str, *err, __LINE__);
	    return 0.0;
	 }				
      } else {
	 return f[n-1] + upper*(x-b); /* linear extrapolation */
      }
   } else {
      return (r - i)*(f[i+1] - f[i]) + f[i]; /* interpolation */
   }
}

/* ============================================================ *
 * like interpol, but f beeing a 2d-function			*
 * 'lower' and 'upper' are the powers of a power law extra-	*
 * polation in the second argument				*
 * ============================================================ */
double sm2_interpol2d(double **f, int nx, double ax, double bx, double dx, double x,
		      int ny, double ay, double by, double dy, double y,
		      double lower, double upper, error **err)
{
   double t, dt, s, ds;
   int i, j;
   char str[1024];

   if (x<ax) {
      sprintf(str, "Query point for 2d interpolation x=%g smaller than minimum a=%g (first coordinate).\n",
	      x, ax);
      *err = addError(math_interpol2small, str, *err, __LINE__);
      return 0.0;
   }
   if (x>bx) {
      sprintf(str, "Query point for 2d interpolation x=%g larger than maximum b=%g (first coordinate).\n",
	      x, bx);
      *err = addError(math_interpol2small, str, *err, __LINE__);
      return 0.0;
   }

   t = (x - ax)/dx;
   i = (int)(floor(t));
   testErrorRet(i+1>nx || i<0, math_interpoloutofrange,
		"2d Interpolation index (first coordinate) out of range",
		*err, __LINE__, 0.0);

   if (i+1==nx) i = i-1;
   dt = t - i;
   if (y < ay) {
      /* Linear extrapolation for second coordinate */
      // MKDEBUG: Check for lower==0??
      return ((1.-dt)*f[i][0] + dt*f[i+1][0]) + (y-ay)*lower;
   } else if (y > by) {
      return ((1.-dt)*f[i][ny-1] + dt*f[i+1][ny-1]) + (y-by)*upper;
   }
   s = (y - ay)/dy;
   j = (int)(floor(s));
   ds = s - j;
   return (1.-dt)*(1.-ds)*f[i][j] + (1.-dt)*ds*f[i][j+1] +
     dt*(1.-ds)*f[i+1][j] + dt*ds*f[i+1][j+1];
}

/* Calculates eigenvectors and -values. From Numerical Recipes */
#define NRANSI
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);

void jacobi_transform(const double *ainput, int n, double d[], double **v, int *nrot, error **err)
{
   int j,iq,ip,i;
   double tresh,theta,tau,t,sm,s,h,g,c,*b,*z, **a;

   /* Transform ainput (matrix as 1d-vector) to 2d-matrix, starting at [1][1] */
   a = sm2_matrix(1, n, 1, n, err);
   forwardError(*err, __LINE__,);
   for (i=0; i<n; i++) {
      for (j=0; j<n; j++) {
	 a[i+1][j+1] = ainput[i*n+j];
      }
   }

   b = sm2_vector(1,n, err);                     forwardError(*err, __LINE__,);
   z = sm2_vector(1,n, err);			 forwardError(*err, __LINE__,);
   for (ip=1;ip<=n;ip++) {
      for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
      v[ip][ip]=1.0;
   }
   for (ip=1;ip<=n;ip++) {
      b[ip]=d[ip]=a[ip][ip];
      z[ip]=0.0;
   }
   *nrot=0;
   for (i=1;i<=50;i++) {
      sm=0.0;
      for (ip=1;ip<=n-1;ip++) {
	 for (iq=ip+1;iq<=n;iq++)
	   sm += fabs(a[ip][iq]);
      }
      if (sm == 0.0) {
	 sm2_free_vector(z,1,n);
	 sm2_free_vector(b,1,n);
	 return;
      }
      if (i < 4)
	tresh=0.2*sm/(n*n);
      else
	tresh=0.0;
      for (ip=1;ip<=n-1;ip++) {
	 for (iq=ip+1;iq<=n;iq++) {
	    g=100.0*fabs(a[ip][iq]);
	    if (i > 4 && (double)(fabs(d[ip])+g) == (double)fabs(d[ip])
		&& (double)(fabs(d[iq])+g) == (double)fabs(d[iq]))
	      a[ip][iq]=0.0;
	    else if (fabs(a[ip][iq]) > tresh) {
	       h=d[iq]-d[ip];
	       if ((double)(fabs(h)+g) == (double)fabs(h))
		 t=(a[ip][iq])/h;
	       else {
		  theta=0.5*h/(a[ip][iq]);
		  t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
		  if (theta < 0.0) t = -t;
	       }
	       c=1.0/sqrt(1+t*t);
	       s=t*c;
	       tau=s/(1.0+c);
	       h=t*a[ip][iq];
	       z[ip] -= h;
	       z[iq] += h;
	       d[ip] -= h;
	       d[iq] += h;
	       a[ip][iq]=0.0;
	       for (j=1;j<=ip-1;j++) {
		  ROTATE(a,j,ip,j,iq)
		    }
	       for (j=ip+1;j<=iq-1;j++) {
		  ROTATE(a,ip,j,j,iq)
		    }
	       for (j=iq+1;j<=n;j++) {
		  ROTATE(a,ip,j,iq,j)
		    }
	       for (j=1;j<=n;j++) {
		  ROTATE(v,j,ip,j,iq)
		    }
	       ++(*nrot);
	    }
	 }
      }
      for (ip=1;ip<=n;ip++) {
	 b[ip] += z[ip];
	 d[ip]=b[ip];
	 z[ip]=0.0;
      }
   }

   *err = addError(math_tooManySteps, "Too many iterations in routine jacobi", *err, __LINE__);
   return;
}
#undef ROTATE
#undef NRANSI

#define NRANSI
#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 50

/* Indexing, from Numerical Recipes */
void indexx(unsigned long n, double arr[], unsigned long indx[], error **err)
{
   unsigned long i,indxt,ir=n,itemp,j,k,l=1;
   int jstack=0,*istack;
   double a;

   istack = malloc_err(sizeof(int)*(NSTACK+1), err);     forwardError(*err, __LINE__,);
   for (j=1;j<=n;j++) indx[j]=j;
   for (;;) {
      if (ir-l < M) {
	 for (j=l+1;j<=ir;j++) {
	    indxt=indx[j];
	    a=arr[indxt];
	    for (i=j-1;i>=l;i--) {
	       if (arr[indx[i]] <= a) break;
	       indx[i+1]=indx[i];
	    }
	    indx[i+1]=indxt;
	 }
	 if (jstack == 0) break;
	 ir=istack[jstack--];
	 l=istack[jstack--];
      } else {
	 k=(l+ir) >> 1;
	 SWAP(indx[k],indx[l+1]);
	 if (arr[indx[l]] > arr[indx[ir]]) {
	    SWAP(indx[l],indx[ir])
	      }
	 if (arr[indx[l+1]] > arr[indx[ir]]) {
	    SWAP(indx[l+1],indx[ir])
	      }
	 if (arr[indx[l]] > arr[indx[l+1]]) {
	    SWAP(indx[l],indx[l+1])
	      }
	 i=l+1;
	 j=ir;
	 indxt=indx[l+1];
	 a=arr[indxt];
	 for (;;) {
	    do i++; while (arr[indx[i]] < a);
	    do j--; while (arr[indx[j]] > a);
	    if (j < i) break;
	    SWAP(indx[i],indx[j])
	      }
	 indx[l+1]=indx[j];
	 indx[j]=indxt;
	 jstack += 2;
	 testErrorRet(jstack>NSTACK, math_stackTooSmall, "NSTACK too small in indexx.", *err, __LINE__,);
	 if (ir-i+1 >= j-l) {
	    istack[jstack]=ir;
	    istack[jstack-1]=i;
	    ir=j-1;
	 } else {
	    istack[jstack]=j-1;
	    istack[jstack-1]=l;
	    l=i;
	 }
      }
   }
   free(istack);
}
#undef M
#undef NSTACK
#undef SWAP
#undef NRANSI

#ifndef _NOFFTW_
/* ============================================================ *
 * Calculates the shear correlation function from the power     *
 * spectrum via a Hankel transform.			        *
 * ============================================================ */
void tpstat_via_hankel(void* self, double **xi, double *logthetamin, double *logthetamax, 
		       tpstat_t tpstat, double (*P_projected)(void *, double, int, int, error**),
		       int i_bin, int j_bin, error **err)
{

   double loglmax, loglmin, dlnl, lnrc, mu, q;
   int nc;
   double        l, kk, *lP, t, norm;
   fftw_plan     plan1, plan;
   fftw_complex *f_lP, *conv;
   fftw_complex  kernel;
   int           i, count;

   lP   = fftw_malloc(N_thetaH*sizeof(double));
   testErrorRet(lP==NULL, io_alloc, "Memory allocation of failed", *err, __LINE__,);

   f_lP = fftw_malloc((N_thetaH/2+1)*sizeof(fftw_complex));
   testErrorRet(f_lP==NULL, io_alloc, "Memory allocation of failed", *err, __LINE__,);

   conv = fftw_malloc((N_thetaH/2+1)*sizeof(fftw_complex));
   testErrorRet(conv==NULL, io_alloc, "Memory allocation of failed", *err, __LINE__,);

   plan  = fftw_plan_dft_r2c_1d(N_thetaH, lP, f_lP, FFTW_ESTIMATE);
   plan1 = fftw_plan_dft_c2r_1d(N_thetaH, conv, lP, FFTW_ESTIMATE);

   loglmax  = log(l_max);
   loglmin  = log(l_min);
   dlnl     = (loglmax-loglmin)/((double)N_thetaH-1.0);
   lnrc     = 0.5*(loglmax+loglmin);
   nc       = N_thetaH/2+1;
  
   /* Power spectrum on logarithmic bins */
   for(i=0; i<N_thetaH; i++) {
      l     = exp(lnrc+(i-nc)*dlnl);
      lP[i] = l*P_projected(self, l, i_bin, j_bin, err);
      forwardError(*err,__LINE__,);
   }

   /* Go to log-Fourier-space */
   fftw_execute(plan);

   /* q: bias, mu: order of Bessel function */
   switch (tpstat) {
      case tp_xipm :
      case tp_w :
	 q =  0.0;
	 mu = 0.0;
	 norm = 2.0*pi;
	 break;
      case tp_gt :
         q =  0.0; 
         mu = 2.0;
         norm = 2.0*pi;
         break;
      case tp_gsqr :
	 q = -2.0;
	 mu = 1.0;
	 norm = 2.0*pi/4.0;
	 break;
      case tp_map2_poly :
	 q = -4.0; 
	 mu = 4.0;
	 norm = 2.0*pi/24.0/24.0;
	 break;
      case tp_map2_gauss :
	 q = mu = 0.0;          /* unused */
	 norm = 4.0*2.0*pi;
	 break;
      case tp_xir :
	 q = mu = 0.0;          /* unused */
	 norm = 2.0*pi*pi;
	 break;
      default :
	 *err = addErrorVA(math_unknown, "Unknown tpstat_t %d", *err, __LINE__, tpstat);
	 return;
   }

   for (count=0; count<=(tpstat==tp_xipm?1:0); count++) {

      if (count==1) mu = 4.0;      /* For xim */

      /* Perform the convolution, negative sign for kernel (complex conj.!) */
      for(i=0; i<N_thetaH/2+1; i++) {
	 kk = 2*pi*i/(dlnl*N_thetaH);

	 switch (tpstat) {
	    case tp_xipm :
	    case tp_w :
	    case tp_gt  :
	       hankel_kernel_mu(kk, &kernel, q, mu, err);
	       break;
	    case tp_gsqr :
	    case tp_map2_poly :
	       hankel_kernel_mumu(kk, &kernel, q, mu, err);
	       break;
	    case tp_map2_gauss :
	       hankel_kernel_exp(kk, &kernel, err);
	       break;
	    case tp_xir :
	       hankel_kernel_tophat(kk, &kernel, err);
	       break;
	    default :
	       *err = addErrorVA(math_unknown, "unknown tpstat %d", *err, __LINE__, tpstat);
	       return;
	 }
	 forwardError(*err, __LINE__,);

	 /* Re and Im */
	 conv[i][0] = f_lP[i][0]*kernel[0]-f_lP[i][1]*kernel[1];
	 conv[i][1] = f_lP[i][1]*kernel[0]+f_lP[i][0]*kernel[1];
      }

      /* Force Nyquist- and 0-frequency-components to be double */
      conv[0][1] = 0;
      conv[N_thetaH/2][1] = 0;

      /* Go back to real space, i labels log-bins in theta */
      fftw_execute(plan1);
      for(i=0; i<N_thetaH; i++) {
	 t = exp((nc-i)*dlnl-lnrc);             /* t = 1/l */
	 xi[count][N_thetaH-i-1] = lP[i]/(t*N_thetaH*norm);
      }
   }

   *logthetamin = (nc-N_thetaH+1)*dlnl-lnrc;
   *logthetamax = nc*dlnl-lnrc;

   /* Clean up */
   fftw_free(conv);
   fftw_free(lP);
   fftw_free(f_lP);
   fftw_destroy_plan(plan);
   fftw_destroy_plan(plan1);
}

/* ============================================================ *
 * Convolution kernel for Hankel transform with bias q,		*
 * Bessel_mu (for xi_pm).					*
 * ============================================================ */
void hankel_kernel_mu(double k, fftw_complex *res, double q, double mu, error **err)
{
   fftw_complex a1, a2, g1, g2;
   double        mod, kln2, si, co, d1, d2, pref;

   /* arguments for complex Gamma */
   a1[0] = 0.5*(1.0+mu+q);
   a2[0] = 0.5*(1.0+mu-q);
   a1[1] = 0.5*k; a2[1]=-a1[1];

   cdgamma(a1,&g1);
   cdgamma(a2,&g2);

   kln2 = k*ln2;
   si   = sin(kln2);
   co   = cos(kln2);
   d1   = g1[0]*g2[0]+g1[1]*g2[1]; /* Re */
   d2   = g1[1]*g2[0]-g1[0]*g2[1]; /* Im */
   mod  = g2[0]*g2[0]+g2[1]*g2[1];
   pref = exp(ln2*q)/mod;

   (*res)[0] = pref*(co*d1-si*d2);
   (*res)[1] = pref*(si*d1+co*d2);

   testErrorRet(!finite((*res)[0])||!finite((*res)[1]), math_infnan, "Inf or nan", *err, __LINE__,);
}

/* ============================================================ *
 * Convolution kernel for Hankel transform with bias q,		*
 * Bessel_mu!2 (for <|g|!2>).					*
 * ============================================================ */
void hankel_kernel_mumu(double k, fftw_complex *res, double q, double mu, error **err)
{
   fftw_complex a1, a2, a3, a4, g1, g2, g3, g4;
   fftw_complex nom, den;
   double pref;

   testErrorRet(k<=-EPSILON2, math_underflow, "k must be larger than EPSILON2", *err, __LINE__,);

   /* arguments for complex Gamma */
   a1[0] = -0.5*q;              a1[1] = -0.5*k;
   a2[0] = 0.5*(1.0+2.0*mu+q);  a2[1] =  0.5*k;
   a3[0] = 0.5*(1.0-q);         a3[1] = -0.5*k;
   a4[0] = 0.5*(1.0+2.0*mu-q);  a4[1] = -0.5*k;

   cdgamma(a1, &g1);
   cdgamma(a2, &g2);
   cdgamma(a3, &g3);
   cdgamma(a4, &g4);

   sm2_Cmul(g1, g2, &nom);
   sm2_Cmul(g3, g4, &den);
   sm2_Cdiv(nom, den, res);

   pref = 1.0/(2.0*sqrt(pi));
   (*res)[0] *= pref;
   (*res)[1] *= pref;

   testErrorRet(!finite((*res)[0]) || !finite((*res)[1]), math_infnan, "Inf or nan", *err, __LINE__,);
}

void hankel_kernel_exp(double k, fftw_complex *res, error **err)
{
   fftw_complex a, g;

   a[0] = 2.5;   a[1] = 0.5*k;
   cdgamma(a, &g);
   (*res)[0] = g[0]*0.5;
   (*res)[1] = g[1]*0.5;

   testErrorRet(!finite((*res)[0]) || !finite((*res)[1]), math_infnan, "Inf or nan", *err, __LINE__,);
}

void hankel_kernel_tophat(double k, fftw_complex *res, error **err)
{
   (*res)[1] = 0.0;
   /* This looks good: */
   //if (fabs(k)<0.5*pi) (*res)[0] = 1.0/2.0/pi;
   /* This should be correct: */
   if (fabs(k)<pi) {
      (*res)[0] = pi;
   } else {
      (*res)[0] = 0.0;
   }

   forwardError(*err, __LINE__,); /* Avoid compiler warning */
}
#endif
