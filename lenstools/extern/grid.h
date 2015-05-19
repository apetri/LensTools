#ifndef __GRID_H
#define __GRID_H

#include <math.h>

int grid2d(double *x,double *y,double *s,double *map,int Nobjects,int Npixel,double map_size);
int grid3d(float *positions,float *weights,double *radius,double *concentration,int Npart,double leftX,double leftY,double leftZ,double sizeX,double sizeY,double sizeZ,int nx,int ny,int nz,float *grid,double(*kernel)(double,double,double,double));
int adaptiveSmoothing(int NumPart,float *positions,float *weights,double *rp,double *concentration,double *binning0, double *binning1,double center,int direction0,int direction1,int normal,int size0,int size1,int projectAll,double *lensingPlane,double(*kernel)(double,double,double,double));

static inline double quadraticKernel(double dsquared,double w,double rv,double c){
	return (1.0/pow(rv,2)) * pow(1.0 - dsquared/pow(rv,2),2);
}

double nfwKernel(double dsquared,double w,double rv,double c);


#endif