/*Python wrapper module for the NICAEA weak lensing code by M.Kilbinger;

The module is called _nicaea and it defines the methods below (see docstrings)
*/

#include <stdlib.h>
#include <getopt.h>

#include <Python.h>
#include <numpy/arrayobject.h>

#include "errorlist.h"
#include "cosmo.h"
#include "lensing.h"
#include "nofz.h"

//Python module docstrings 
static char module_docstring[] = "This module provides a python interface to the NICAEA computations";
static char shearPowerSpectrum_docstring[] = "Compute the shear power spectrum";

//Method declarations
static PyObject *_nicaea_shearPowerSpectrum(PyObject *self,PyObject *args);

//_nicaea method definitions
static PyMethodDef module_methods[] = {

	{"shearPowerSpectrum",_nicaea_shearPowerSpectrum,METH_VARARGS,shearPowerSpectrum_docstring},
	{NULL,NULL,0,NULL}

} ;

//_nicaea constructor
PyMODINIT_FUNC init_nicaea(void){

	PyObject *m = Py_InitModule3("_nicaea",module_methods,module_docstring);
	if(m==NULL) return;

	/*Load numpy functionality*/
	import_array();

}

//////////////////////////////////////////////////
/*Function implementations using backend C code*/
/////////////////////////////////////////////////

//shearPowerSpectrum() implementation
static PyObject *_nicaea_shearPowerSpectrum(PyObject *self,PyObject *args){

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

	char stringerr[4192];

	model=init_parameters_lens(Om,Ode,w0,w1,poly_W,NPoly,H100,Omegab,Omeganu,Neff,si8,ns,nzbins,Nnz,nofz,par_nz,computation_type,transfer_function,growth,dark_energy,norm_mode,tomo_all,sreduced,Q_MAG_SIZE,err);
	printf("%le\n",Pshear(model,10000.0,0,0,err));
	if(isError(*err)){
		stringError(stringerr,*err);
		PyErr_SetString(PyExc_ValueError,stringerr);
		free_parameters_lens(&model);
		return NULL;
	}
	
	free_parameters_lens(&model);

	Py_RETURN_NONE;

}