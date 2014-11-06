#ifndef __AZIMUTH_H
#define __AZIMUTH_H

#include <complex.h>

int azimuthal_rfft2(double _Complex *ft_map1,double _Complex *ft_map2,long size_x,long size_y,double map_angle_degrees,int Nvalues,double *lvalues,double *power_l);
int azimuthal_rfft3(double _Complex *ft_map1,double _Complex *ft_map2,long size_x,long size_y,long size_z,double kpixX,double kpixY,double kpixZ,int Nvalues,double *kvalues,double *power_k,long *hits);

#endif