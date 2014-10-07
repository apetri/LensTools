#ifndef __DIFFERENTIALS_H
#define __DIFFERENTIALS_H

void gradient_xy(double *map,double *grad_map_x,double *grad_map_y,long map_size,int Npoints,int *x_points,int *y_points);
void hessian(double *,double *,double *,double *,long);

#endif