#ifndef __PEAKS_H
#define __PEAKS_H

void peak_count(double *map,unsigned char *mask,long map_size, double sigma, int Nthresh, double *thresholds, double *peaks);
int peak_locations(double *map,unsigned char *mask,long map_size, double sigma, int Nthresh, double *thresholds,double *values,int *locations_x,int *locations_y);

#endif