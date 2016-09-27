#ifndef __DIFFERENTIALS_H
#define __DIFFERENTIALS_H

void gradient_xy(double *map,double *grad_map_x,double *grad_map_y,long map_size,int Npoints,int *x_points,int *y_points);
void hessian(double *map,double *hess_xx_map,double *hess_yy_map,double *hess_xy_map,long map_size,int Npoints, int *x_points,int *y_points);
void gradLaplacian(double *map,double *grad_map_x,double *grad_map_y,long map_size,int Npoints,int *x_points,int *y_points);

#endif