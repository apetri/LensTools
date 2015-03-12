#ifndef __DARKENERGY_SUPPORT_H
#define __DARKENERGY_SUPPORT_H

/*
 *  darkenergy_support.h
 *  Dark Energy Extension for V. Springel's Gadget-2
 *
 *  Created by Jan Michael Kratochvil - on 7/16/07.
 *  Copyright 2007 Jan Michael Kratochvil. All rights reserved.
 *
 */


// will be included at the beginning of main.c
// list (empty) declaration here of all subroutines that are defined in integration.h, so they appear at the beginning of the main.c file, due to the include.

//////////////////////////////////
// INTEGRATION:
//////////////////////////////////

int odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1,
	double hmin, int *nok, int *nbad,
	void (*derivs)(double, double [], double []),
	void (*rkqs)(double [], double [], int, double *, double, double, double [],
	double *, double *, void (*)(double, double [], double [])));
	
void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs)(double, double [], double []));
	
void rkck(double y[], double dydx[], int n, double x, double h, double yout[],
	double yerr[], void (*derivs)(double, double [], double []));


////////////////////////
// INTERPOLATION:
////////////////////////

void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y);


////////////////////////////////
// ALLOCATION:
////////////////////////////////

double* Vector(int m);
double** Matrix(int m, int n);
void free_Vector(double* vector);
void free_Matrix(double** matrix, int m);

#endif
