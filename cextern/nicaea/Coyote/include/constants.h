/* ============================================================ *
 * constants.h							*
 * Definition of constants for the Coyote emulator of the non-  *
 * linear power spectrum.					*
 * ============================================================ */

#ifndef __CONSTANTS_H
#define __CONSTANTS_H

#define m    37
#define neta 6000
#define p    5
#define peta 5
#define rs   6

extern const double kemu[1000];
extern const double mean[neta];
extern const double K[neta][peta];
extern const double x[m][p], xmin[p], xrange[p], aemu[rs];
extern const double lamws[peta], lamz[peta];
extern const double beta[peta][p];
extern const double w_coyote[peta][m];
extern const double KrigBasis[peta][m];

#endif