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
	
	double Om=0.26;
	double Ode=0.74;
	double w0=-1.0;
	double w1=0.0;
	double *poly_W=NULL;
	int NPoly=0;
	double H100=0.7;
	double Omegab=0.046;
	double Omeganu=0.0;
	double Neff=0.0;
	double si8=0.8;
	double ns=0.960;
	double nzbins=1;

	int Nnz[] = {3};
	const nofz_t nofz[] = {hist};
	double par_nz[] = {2.0,2.01,10};
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

	model=init_parameters_lens(Om,Ode,w0,w1,poly_W,NPoly,H100,Omegab,Omeganu,Neff,si8,ns,nzbins,Nnz,nofz,par_nz,computation_type,transfer_function,growth,dark_energy,norm_mode,tomo_all,sreduced,Q_MAG_SIZE,err);
	printf("%le\n",Pshear(model,10000.0,0,0,err));
	if(isError(*err)){
		printError(stderr,*err);
		exit(getErrorValue(*err));
	}
	
	free_parameters_lens(&model);

	return 0;

}