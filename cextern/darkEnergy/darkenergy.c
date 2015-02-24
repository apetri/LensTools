/*
 *  darkenergy.c
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


//#include <nrutil.h>
//#define NRANSI

// Global variables for dark energy:
double w0, wa;


// Global variables for integration:
double dxsav; // distance at which steps are to be saved.
double *xp; // array of "times" at which output is saved ("time" = the evolution parameter of the differential equation).
double **yp; // array containing the values of the funtions at the output times given in array xp. 
int kmax; // maximal number of intermediate steps saved (last step is always saved)
int kount; // kounts through saved steps up to kmax
int nrhs;   // counts function evaluations (increased by one each time derivs is called, which contains the ODE's) 
int neqs; // number of differential equations (needs to be global for free_darkenergy(), otherwise just use as local variable in initialize_darkenergy().

double *y2;
double *ap;
double **DEp;

// This function specifies the ordinary differential equations.
// They must by first order ODE's, higher order equations must be rewritten as a system of ordinary first order equations.
// First derivative is on the left-hand side, the corresponding equation on the right-hand side.
// Functions are enumerated from 1 to neqs, where neqs specifies the number of equations and must be set properly within the (sub)routine that calls odeint. 
void derivs (double x, double y[], double dydx[])
{
	nrhs++; //counts function evaluations
	// dydx=f(y(x),x) first order ordinary differential equation:
	dydx[1]=(1+w(x))/(1+x); // arbitrary examples of first order differential equations being solved (this line and the next two).
}


void initialize_darkenergy (void) {
	
	int i;
	// int neqs; // number of differential equations
	//double ystart[neqs+1];
	double *ystart; // initial conditions array
	double x1, x2; // starting and end point of integration
	double eps, h1, hmin; // Performance control parameters for numerical integrator odeint.
	int nok, nbad; // counts number of good and bad steps (is passed on as a pointer to subroutines called by odeint, so can be modified by those correctly).

	// Number of ordinary differential equations to be solved:
	neqs=1;
	// The Differential Equations are specified in the function derivs.

	// Performance and Output Control Parameters for numerical integrator:
	// Performance:
	eps=pow(10,-18); // Precision, maximal allowed error
	h1=0.01; // guess for first stepsize
	hmin=0; // minimal stepsize (can be zero)
	// Output (output stored in (xp, yp[])):
	kmax=10000000; // maximum number of intermediate steps stored (first one and last one are always stored, and count towards the total number of steps specified by kmax).
	dxsav=0.0001; // steps saved only in intervals larger than this value.

	// Allocate arrays for differential equation ("time" parameter if the equation is xp, the functions are enumerated by yp[1-neqs]): 
	xp=Vector(kmax); // Initializes vector, ready for NR-C (Numerical Recipes in C) component enumeration 1 through kmax.
	yp=Matrix(neqs,kmax); // Initializes neqs x kmax matrix with NR-C enumeration. 
	ystart=Vector(neqs); // Initial conditions (position) for functions solving the differential equations (only one in example here).
	// WARNING: NEVER call xp[0], yp[0][...], or ystart[0] !!! Count starts at 1. (Otherwise you will overwrite some other variables!)

	// Allocate dark energy array:
	ap=Vector(kmax);
	DEp=Matrix(neqs,kmax);
		
	//Initial conditions (for first oder equation example here, only one starting value, no derivative, needed):
	ystart[1]=0.0; // function value of first ODE at starting point is 0, because it's an integral.
	x1=0.0; // starting point of integration is at redshift 0.
	x2=170.0; // end point of integration at this redshift (redshift needs to be larger than begin of N-body simulation, make higher if necessary).
	
	// printf("Kount before odeint: %d.\n", kount);
	
	// Call driver for numerical integrator with above parameters (the driver calls then further subroutines):
	odeint(ystart, neqs, x1, x2, eps, h1, hmin, &nok, &nbad, derivs, rkqs);

	// printf("Kount: %d.\n", kount);

	// Sample output to check that everything is o.k. and demonstrate how the integrator works:
	// Output should be only correct for writeouts from i=1 to i=kmax. The rest is included just as a reference.
	
	// printf("ystart, nok, nbad, nrhs: %e %d %d %d\n", ystart[1], nok, nbad, nrhs);

	// Before splining, replace redshift z by scale factor a (the name of the variable is xp), and integral by whole dark energy factor expression, reorder by ascending scale factor:
	for (i=1;i<=kount;i++)
	{
		ap[i]=1.0/(1.0+xp[kount+1-i]);
		DEp[1][i]=exp(3*yp[1][kount+1-i]);
		// printf("Eq 1: i, a (scale factor), xp (redshift), DEp (dark energy e-factor): %d --  %e %e %e\n", i, ap[i], xp[kount+1-i], DEp[1][i]);
	}
	//Now can spline this final expression as a function of scale factor:

	// Now do interpolation of above tabulated solutions:
	y2=Vector(kmax);
	// Initialize spline (need only do once):
	// Arguments for spline(): table of arguments of function (x), table of function values at those arguments (y(x)), number of points tabulated, first derivatives at first and last point, output: second derivatives of function at tabulated points).
	//spline(xp, yp[1], kmax, yp[1][1], yp[1][kmax], y2);
	spline(ap, DEp[1], kount, DEp[1][1], DEp[1][kount], y2);
	// Finding the first derivatives at start and end point is easy, because the ODE in derivs is always given as a first order equation, so it's always just the right-hand side of the equation, evaluated at [1] and [kmax] respectively.

/*
	FILE *output_file1;
	output_file1 = fopen ("DE_factor.txt", "w");
	
	for (i=1;i<=99;i++)
	{
		xx=i*0.01+0.01;
		// Evaluate spline:  
		// Arguments for splint(): table of arguments of function (x), table of function values at those arguments (y(x)), output table from function spline above (second derivatives probably), number of tabulated points, point at which splined function is to be evaluated (x), output needs to be supplied as address and is the value of the splined function at that point (y(x)).
		splint(ap, DEp[1], y2, kount, xx, &yy);
	
		fprintf(output_file1, "%e %e \n", xx, yy);
	}
	
	fclose(output_file1);
*/	

	

	free_Vector(xp);
	free_Matrix(yp, neqs);
	free_Vector(ystart);
	
	// Do not free those until the very end of the whole Gadget run:
	//free_Vector(y2);
	//free_Vector(ap);
	//free_Matrix(DEp, neqs);
		
    printf("Finished initializing dark energy.\n");
    return;
}

double w(double z)
{
	double ww; //, w0, wa;
	
	//w0=-1.0; wa=0.0;
	ww=w0+(z/(1+z))*wa; // example of a redshift-dependent dark energy equation of state parameter w(z).
	return ww;

}

double DarkEnergy (double a)
{
	double yy;
	splint(ap, DEp[1], y2, kount, a, &yy);
	return yy;
}


double free_darkenergy(void) // frees any dark energy calculation arrays not freed so far.
{
	free_Vector(y2);
	free_Vector(ap);
	free_Matrix(DEp, neqs);
	return 0.0;
}

//#undef NRANSI
