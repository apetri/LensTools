#ifndef __INTERFACE_H
#define __INTERFACE_H

/* note the underscore! */
int f77main_(double * hI, double * omegamI, double * omegavI, double * omegakI, double * omegaQI, double * wQI, double * wQpI, double * z1, double * z2, double * D1, double *D2, double *zmaxact, double *zminact, int *iwmodeI);
double Dplus_Interface(double z2, double z1, double omegamI, double omegaQI, double omegavI, double wQI, double wQpI, double hI, double *vel_prefac_lam);

#endif
