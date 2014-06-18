#ifndef __DIFFERENTIALS_H
#define __DIFFERENTIALS_H

float average_map(float *,long);
void average_subtract(float *,long);
float variance_map(float *,long);
void gradient(float *,float *,long);
void gradient_xy(float *,float *,float *,long);
void hessian(float *,float *,float *,float *,long);

#endif