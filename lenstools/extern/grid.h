#ifndef __GRID_H
#define __GRID_H

#include <math.h>

int grid3d(float *positions,int Npart,double leftX,double leftY,double leftZ,double sizeX,double sizeY,double sizeZ,int nx,int ny,int nz,float *grid);
int adaptiveSmoothing(int NumPart,float *positions,double *rp,double *binning0, double *binning1,double center,int direction0,int direction1,int normal,int size0,int size1,int projectAll,double *lensingPlane);
static inline double quadraticKernel(double d,double rp){
	return (1.0/pow(rp,2)) * pow(1.0 - d/pow(rp,2),2);
}


#endif