#ifndef __GRID_H
#define __GRID_H

int grid3d(float *positions,int Npart,double leftX,double leftY,double leftZ,double sizeX,double sizeY,double sizeZ,int nx,int ny,int nz,float *grid);
int adaptiveSmoothing(int NumPart,float *positions,double *rp,double *binning0, double *binning1,int direction0,int direction1,int normal,int size0,int size1,double *lensingPlane);


#endif