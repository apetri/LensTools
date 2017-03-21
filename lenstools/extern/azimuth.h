#ifndef __AZIMUTH_H
#define __AZIMUTH_H

#include <complex.h>

int azimuthal_rfft2(double _Complex *ft_map1,double _Complex *ft_map2,long size_x,long size_y,double map_angle_degrees,int Nvalues,double *lvalues,double *power_l,double *scale);
int azimuthal_rfft3(double _Complex *ft_map1,double _Complex *ft_map2,long size_x,long size_y,long size_z,double kpixX,double kpixY,double kpixZ,int Nvalues,double *kvalues,double *power_k,long *hits);

int bispectrum(double _Complex *ft_map1,double _Complex *ft_map2,double _Complex *ft_map3,long size_x,long size_y,double map_angle_degrees,int Nvalues,double *lvalues,double *bispectrum_l,int (*k1tok2)(int,int,int*,int*,void*),void *args);

int k1tok2_equilateral(int kx1,int ky1,int *kx2,int *ky2,void *args);
int bispectrum_equilateral(double _Complex *ft_map1,double _Complex *ft_map2,double _Complex *ft_map3,long size_x,long size_y,double map_angle_degrees,int Nvalues,double *lvalues,double *bispectrum_l);


int k1tok2_folded(int kx1,int ky1,int *kx2,int *ky2,void *args);
int bispectrum_folded(double _Complex *ft_map1,double _Complex *ft_map2,double _Complex *ft_map3,long size_x,long size_y,double map_angle_degrees,double folding_ratio,int Nvalues,double *lvalues,double *bispectrum_l);


#endif