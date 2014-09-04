#ifndef __DESIGN_H
#define __DESIGN_H

#include <gsl/gsl_matrix.h>

//function prototypes
double distance(double *x,double *y,int D,double p);
double cost(gsl_matrix *data,int Npoints,int D,double p,double lambda);
double diagonalCost(int Npoints,double lambda);
double swap(double **data,int Npoints,int D,double p,double lambda,int i1,int i2, int d);
void swapBack(double **data,int i1,int i2,int d);

//main design sampler prototype
double sampler(int Npoints,int D,double p,double lambda,int seed,int maxIterations,double **data,double *costValues);

#endif