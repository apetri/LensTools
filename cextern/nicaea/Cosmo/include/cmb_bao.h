/* ============================================================ *
 * cmb_bao.h							*
 * Martin Kilbinger						*
 * ============================================================ */

#ifndef __CMB_BAO_H
#define __CMB_BAO_H


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "cosmo.h"
#include "errorlist.h"
#include "io.h"
#include "maths.h"
#include "mvdens.h"


double z_star(cosmo *self);
double acoustic_scale(cosmo *self, error **err);
double shift_parameter(cosmo *self, error **err);
double D_V(cosmo *self, double a, error **err);

double chi2_cmbDP(cosmo *model, mvdens *g, error **err);
double chi2_bao(cosmo *model, mvdens *g, error **err);
double chi2_bao_only(cosmo *model, mvdens *g, error **err);
double chi2_bao_A(cosmo *model, mvdens *g, const double *z_BAO, error **err);
double chi2_bao_d_z(cosmo *model, mvdens *g, const double *z_BAO, error **err);
double chi2_bao_D_V_ratio(cosmo *model, mvdens *g, const double *z_BAO, error **err);


#endif
