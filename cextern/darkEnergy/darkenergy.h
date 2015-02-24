#ifndef __DARKENERGY_H
#define __DARKENERGY_H

/*
 *  darkenergy.h
 *  Dark Energy Extension for V. Springel's Gadget-2
 *
 *  Created by Jan Michael Kratochvil - on 7/16/07.
 *  Copyright 2007 Jan Michael Kratochvil. All rights reserved.
 *
 */


// Must be included in all .c files that are not main.c
// Contains all references to global variables and functions which are declared in main.c

extern double w0, wa; // parameters for dark energy (must be global only if don't want to change them hardcoded in function w(a) upon compilation time.

extern int kmax, kount, nrhs, neqs;
extern double dxsav, *xp, **yp;
extern double *y2;
extern double *ap;
extern double **DEp;

void derivs (double x, double y[], double dydx[]);
void initialize_darkenergy (void);

double w(double z);
double DarkEnergy (double a);

double free_darkenergy(void);

#endif
