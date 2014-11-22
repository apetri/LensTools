#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "errorlist.h"
#include "io.h"
#include "maths.h"
#include "cosmo.h"
#include "coyote.h"
#include "cmb_bao.h"
#include "lensing.h"
#include "nofz.h"

int main(){


	cosmo_lens *model;
	int Nnz[] = {2};
	const nofz_t nofz[] = {single};
	double par_nz[] = {2.0,2.0};
	nonlinear_t computation_type=smith03;
	transfer_t transfer_function=eisenhu;
	growth_t growth=growth_de;
	de_param_t dark_energy=linder;
	norm_t norm_mode=norm_s8;
	tomo_t tomography=tomo_all;
	reduced_t sreduced=reduced_none;
	double Q_MAG_SIZE=1.0;
	error *myerr=NULL,**err;
	err=&myerr;

	model=init_parameters_lens(0.26,0.74,-1.0,0.0,NULL,0,0.7,0.046,0.0,0.0,1.0,0.96,1,Nnz,nofz,par_nz,computation_type,transfer_function,growth,dark_energy,norm_mode,tomo_all,sreduced,Q_MAG_SIZE,err);
	printf("%le\n",Pshear(model,10000.0,0,0,err));
	quitOnError(*err,__LINE__,stderr);
	
	free_parameters_lens(&model);

	return 0;

}