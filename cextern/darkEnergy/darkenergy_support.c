/*
 *  darkenergy_support.c
 *  Dark Energy Extension for V. Springel's Gadget-2
 *
 *  Created by Jan Michael Kratochvil - on 7/16/07.
 *  Copyright 2007 Jan Michael Kratochvil. All rights reserved.
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "darkenergy.h"
#include "darkenergy_support.h"


///////////////////////////////////////////
// INTEGRATION:
///////////////////////////////////////////


//#define NRANSI
//#include "nrutil.h"
#define MAXSTP 1000000
// orignial setting from Numerical Recipes in C was MAXSTP 10000
#define TINY 1.0e-30

// Macros for definitions of simple functions that are called by the integration subroutines from Numerical Recipes in C.
// They have been here reverse engineered based on Numerical Recipes in FORTRAN.
#define FMAX(A, B) (((A)>(B)) ? (A) : (B))
#define FMIN(A, B) (((A)>(B)) ? (B) : (A))
#define SIGN(A, B) ((((B)>=0.0) ? 1.0 : -1.0) * fabs(A))



void odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1,
	double hmin, int *nok, int *nbad,
	void (*derivs)(double, double [], double []),
	void (*rkqs)(double [], double [], int, double *, double, double, double [],
	double *, double *, void (*)(double, double [], double [])))
{
	int nstp,i;
	double xsav,x,hnext,hdid,h;
	double *yscal,*y,*dydx;

	yscal=Vector(nvar);
	y=Vector(nvar);
	dydx=Vector(nvar);

	//double yscal[nvar+1], y[nvar+1], dydx[nvar+1];
	//yscal=vector(1,nvar);
	//y=vector(1,nvar);
	//dydx=vector(1,nvar);
	x=x1;
	h=SIGN(h1,x2-x1);  // returns magnitude of h1 with the sign of x2-x1
	*nok = (*nbad) = kount = 0;
	for (i=1;i<=nvar;i++) y[i]=ystart[i];
	if (kmax > 0) xsav=x-dxsav*2.0;
	for (nstp=1;nstp<=MAXSTP;nstp++) {
		(*derivs)(x,y,dydx);
		for (i=1;i<=nvar;i++)
			yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
		if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) {
			kount++;
			xp[kount]=x;
			for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
			xsav=x;
		}
		if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
		(*rkqs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);
		if (hdid == h) (*nok)++; else (*nbad)++;
		if ((x-x2)*(x2-x1) >= 0.0) {
			for (i=1;i<=nvar;i++) ystart[i]=y[i];
			if (kmax) {
				kount++;
				xp[kount]=x;
				for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
			}
		//	free_vector(dydx,1,nvar);
		//	free_vector(y,1,nvar);
		//	free_vector(yscal,1,nvar);
		free_Vector(dydx);
		free_Vector(y);
		free_Vector(yscal);

			return;
		}
		if (fabs(hnext) <= hmin)
		{
			printf("Step size too small in odeint");
			exit(1);
		}
		h=hnext;
	}
	printf("Too many steps in routine odeint");
	exit(1);
}
#undef MAXSTP
#undef TINY
//#undef NRANSI



//#include <math.h>
//#define NRANSI
//#include "nrutil.h"
#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4

void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs)(double, double [], double []))
{
	void rkck(double y[], double dydx[], int n, double x, double h,
		double yout[], double yerr[], void (*derivs)(double, double [], double []));
	int i;
	double errmax,h,htemp,xnew,*yerr,*ytemp;
		
	//yerr=vector(1,n);
	//ytemp=vector(1,n);
	yerr=Vector(n);
	ytemp=Vector(n);
	h=htry;
	for (;;) {
		rkck(y,dydx,n,*x,h,ytemp,yerr,derivs);
		errmax=0.0;
		for (i=1;i<=n;i++) errmax=FMAX(errmax,fabs(yerr[i]/yscal[i]));
		errmax /= eps;
		if (errmax <= 1.0) break;
		htemp=SAFETY*h*pow(errmax,PSHRNK);
		h=(h >= 0.0 ? FMAX(htemp,0.1*h) : FMIN(htemp,0.1*h));
		xnew=(*x)+h;
		if (xnew == *x)
		{
			printf("stepsize underflow in rkqs");
			exit(1);
		}
	}
	if (errmax > ERRCON) *hnext=SAFETY*h*pow(errmax,PGROW);
	else *hnext=5.0*h;
	*x += (*hdid=h);
	for (i=1;i<=n;i++) y[i]=ytemp[i];
	//free_vector(ytemp,1,n);
	//free_vector(yerr,1,n);
	free_Vector(ytemp);
	free_Vector(yerr);
}
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON
//#undef NRANSI


//#define NRANSI
//#include "nrutil.h"

void rkck(double y[], double dydx[], int n, double x, double h, double yout[],
	double yerr[], void (*derivs)(double, double [], double []))
{
	int i;
	static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
		b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
		b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
		b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
		b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
		c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
		dc5 = -277.00/14336.0;
	double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
		dc4=c4-13525.0/55296.0,dc6=c6-0.25;
	double *ak2,*ak3,*ak4,*ak5,*ak6,*ytemp;

	//ak2=vector(1,n); // from original Numerical Recipes in C, replaced vector allocation subroutine by my own. 
	//ak3=vector(1,n);
	//ak4=vector(1,n);
	//ak5=vector(1,n);
	//ak6=vector(1,n);
	//ytemp=vector(1,n);
	ak2=Vector(n);
	ak3=Vector(n);
	ak4=Vector(n);
	ak5=Vector(n);
	ak6=Vector(n);
	ytemp=Vector(n);
	
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+b21*h*dydx[i];
	(*derivs)(x+a2*h,ytemp,ak2);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
	(*derivs)(x+a3*h,ytemp,ak3);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
	(*derivs)(x+a4*h,ytemp,ak4);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
	(*derivs)(x+a5*h,ytemp,ak5);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
	(*derivs)(x+a6*h,ytemp,ak6);
	for (i=1;i<=n;i++)
		yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
	for (i=1;i<=n;i++)
		yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
	//free_vector(ytemp,1,n);
	//free_vector(ak6,1,n);
	//free_vector(ak5,1,n);
	//free_vector(ak4,1,n);
	//free_vector(ak3,1,n);
	//free_vector(ak2,1,n);
	free_Vector(ytemp);
	free_Vector(ak6);
	free_Vector(ak5);
	free_Vector(ak4);
	free_Vector(ak3);
	free_Vector(ak2);

}
//#undef NRANSI


/////////////////////////////////////
// INTERPOLATION:
/////////////////////////////////////

void spline(double x[], double y[], int n, double yp1, double ypn, double y2[])
{
	int i,k;
	double p,qn,sig,un,*u;

	u=Vector(n-1);
	//u=vector(1,n-1);
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
	for (k=n-1;k>=1;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
	free_Vector(u);
	//free_vector(u,1,n-1);
}
//#undef NRANSI


void splint(double xa[], double ya[], double y2a[], int n, double x, double *y)
{
	//void nrerror(char error_text[]);
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
	if (h == 0.0) 
	{
		printf("Bad xa input to routine splint");
		//nrerror("Bad xa input to routine splint");
		exit(1);
	}
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}


///////////////////////////////////////
// ALLOCATION:
///////////////////////////////////////

// Subroutines for convenient initialization of arrays (1D vectors and 2D matrices)
// complying with the enumeation convention of numerical recipes, staring to count at 1.

double* Vector(int m) // Initializes 2D array (declare array by "double **arrayname;")
// with m rows and n columns. Counts start at 1 and run through m and n, respectively.
// This function was written to initialize 2D arrays for use with Numerical Recipes in C,
// which use this enumeration of components.
{
         int row;
		 double *vector = malloc(m * sizeof(double));
         assert(vector != NULL);
		 vector--;
         for (row = 1; row <= m; row++) {
			// Initialize new vector with zeroes:
               vector[row] = 0.0;
         } /* END for */
		return vector;
}

double** Matrix(int m, int n) // Initializes 2D array (declare array by "double **arrayname;")
// with m rows and n columns. Counts start at 1 and run through m and n, respectively.
// This function was written to initialize 2D arrays for use with Numerical Recipes in C,
// which use this enumeration of components.
{
         int row, col;
		 double **matrix = malloc(m * sizeof(double *));
         assert(matrix != NULL);
		 matrix--;
         for (row = 1; row <= m; row++) {
            matrix[row] = malloc(n * sizeof(double));
            assert(matrix[row] != NULL);
			matrix[row]--;
			// Initialize new matrix with zeroes:
            for (col = 1; col <= n; col++)
               matrix[row][col] = 0.0;
         } /* END for */
		return matrix;
}


// Memory freeing routines for above allocations of Vector and Matrix:

void free_Vector(double* vector)
{
	vector++; // compensate for offset "vector--;" given at allocation of vector, to conform with Numerical Recipes components enumeration, starting at 1 rather than 0.
	free(vector);
}

void free_Matrix(double** matrix, int m) // need to specify how many rows original matrix had in this routine! (could be improved) 
{
	 // int num_rows;
	 // num_rows = sizeof(matrix)/2; // DOESN'T WORK, ALWAYS RETURNS 2 (sizeof returns number of bytes; even so, size of ints, floats and doubles can be machine dependent!) 
     // printf("Extract number of rows in Matrix: %d\n", num_rows); // DOESN'T WORK - MANUAL NUMBER OF ROWS READ IN HAS BEEN IMPLEMENTED
	int row;
	for (row = 1; row <= m; row++)
	{
	matrix[row]++; // compensate for offset "matrix[row]--;" at matrix allocation.
	free(matrix[row]);
	}
	matrix++;
	free(matrix);
}

#undef FMAX
#undef FMIN
#undef SIGN
