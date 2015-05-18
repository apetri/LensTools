#ifndef __GRID_H
#define __GRID_H

#include <math.h>

int grid2d(double *x,double *y,double *s,double *map,int Nobjects,int Npixel,double map_size);
int grid3d(float *positions,float *weights,int Npart,double leftX,double leftY,double leftZ,double sizeX,double sizeY,double sizeZ,int nx,int ny,int nz,float *grid);
int adaptiveSmoothing(int NumPart,float *positions,double *rp,double *binning0, double *binning1,double center,int direction0,int direction1,int normal,int size0,int size1,int projectAll,double *lensingPlane,double(*kernel)(double,double,double));

static inline double quadraticKernel(double dsquared,double rs,double w){
	return (1.0/pow(rs,2)) * pow(1.0 - dsquared/pow(rs,2),2);
}


#endif