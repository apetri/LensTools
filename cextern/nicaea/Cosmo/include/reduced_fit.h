#ifndef __REDUCED_FIT_H
#define __REDUCED_FIT_H


#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/* Error codes */
#define reduced_base           -2100
#define reduced_par            -1 + reduced_base
#define reduced_limit          -2 + reduced_base
#define reduced_fourier_limit  -3 + reduced_base
#define reduced_realsp_limit   -4 + reduced_base

/* Matrix dimensions */
#define M_PAR  8
#define N_POLY 8
#define N_PL   4
#define N_B    4
#define N_C    4


extern const double B_fit[M_PAR][N_PL][N_B];
extern const double C_fit[M_PAR][N_POLY][N_C];
extern const double limits_lower[M_PAR];
extern const double limits_upper[M_PAR];

/* Limits for asymptotic power-law fits, y=ln(l) */
#define Y_LOW  2.0
#define Y_UP   11.5

/* Numerical derivatives */
#define FH 0.01


double h_piece(double logell, const double b[]);

double sum_B_a(int alpha, int i, double a);
double sum_C_a(int alpha, int i, double a);
double Q(int alpha, double logell, double a);
double h_piece(double logell, const double b[]);
double Q(int alpha, double logell, double a);
double sum_a_for_Pg1(double logell, double a_min, int Na, double da, const double *fmn, const double **dfmn_dp,
		     const double dpar[M_PAR]);
int check_limits(const double dpar[M_PAR]);


#endif
