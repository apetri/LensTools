/*Python wrapper module for the NICAEA weak lensing code by M.Kilbinger;

The module is called _nicaea and it defines the methods below (see docstrings)
*/

#include <stdlib.h>
#include <getopt.h>
#include <assert.h>
#include <string.h>

#include <Python.h>
#include <numpy/arrayobject.h>

#include "errorlist.h"
#include "cosmo.h"
#include "lensing.h"
#include "nofz.h"

//Python module docstrings 
static char module_docstring[] = "This module provides a python interface to the NICAEA computations";
static char shearPowerSpectrum_docstring[] = "Compute the shear power spectrum";

//Useful method for parsing Nicaea class attributes into cosmo_lens structs
static cosmo_lens *parse_model(PyObject *args, error **err);

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

/////////////////////////////////
/*Implementation of parse_model*/
/////////////////////////////////

static cosmo_lens *parse_model(PyObject *args, error **err){

	//Cosmological parameters
	double Om,Ode,w0,w1,H100,Omegab,Omeganu,Neff,si8,ns;
	int nzbins;
	int i;
	char *distr_type;

	//Multipoles/angles, redshift distribution and other settings
	PyObject *spec_obj,*Nnz_obj,*nofz_obj,*par_nz_obj,*settings_dict;

	//Parse the input tuple
	if(!PyArg_ParseTuple(args,"ddddddddddiOOOOO",&Om,&Ode,&w0,&w1,&H100,&Omegab,&Omeganu,&Neff,&si8,&ns,&nzbins,&spec_obj,&Nnz_obj,&nofz_obj,&par_nz_obj,&settings_dict)){
		fprintf(stderr,"Input tuple format doesn't match signature!");
		return NULL;
	}

	///////////////////////////////////////////////////////////////////////
	/////////////////////Redshift info/////////////////////////////////////
	///////////////////////////////////////////////////////////////////////

	//Parse Nnz_obj and par_nz_obj into numpy arrays
	PyObject *Nnz_array = PyArray_FROM_OTF(Nnz_obj,NPY_INT32,NPY_IN_ARRAY);
	PyObject *par_nz_array = PyArray_FROM_OTF(par_nz_obj,NPY_DOUBLE,NPY_IN_ARRAY);
	
	if(Nnz_array==NULL || par_nz_array==NULL){
		Py_XDECREF(Nnz_array);
		Py_XDECREF(par_nz_array);
		return NULL;
	} 

	int *Nnz = (int *)PyArray_DATA(Nnz_array);
	double *par_nz = (double *)PyArray_DATA(par_nz_array);

	//Safety
	assert(nzbins==(int)PyArray_DIM(Nnz_array,0));

	//Parse redshift distribution information
	nofz_t nofz[nzbins];
	for(i=0;i<nzbins;i++){

		PyArg_Parse(PyList_GetItem(nofz_obj,i),"s",&distr_type);
		fprintf(stderr,"%d=%s\n",i,distr_type);
		
		if(strcmp(distr_type,"ludo")==0){
			nofz[i]=ludo;
		} else if(strcmp(distr_type,"jonben")==0){
			nofz[i]=jonben;
		} else if(strcmp(distr_type,"ymmk")==0){
			nofz[i]=ymmk;
		} else if(strcmp(distr_type,"hist")==0){
			nofz[i]=hist;
		} else if(strcmp(distr_type,"single")==0){
			nofz[i]=single;
		} else{
			fprintf(stderr,"Distribution %s not implemented",distr_type);
			Py_DECREF(Nnz_array);
			Py_DECREF(par_nz_array);
			return NULL;
		}


	}

	//Set these to default for the moment
	nonlinear_t nonlinear_type=smith03;
	transfer_t transfer_function=eisenhu;
	growth_t growth=growth_de;
	de_param_t dark_energy=linder;
	norm_t norm_mode=norm_s8;
	tomo_t tomography=tomo_all;
	reduced_t sreduced=reduced_none;
	double Q_MAG_SIZE=1.0;

	//cosmo model object
	cosmo_lens *model=init_parameters_lens(Om,Ode,w0,w1,NULL,0,H100,Omegab,Omeganu,Neff,si8,ns,nzbins,Nnz,nofz,par_nz,nonlinear_type,transfer_function,growth,dark_energy,norm_mode,tomography,sreduced,Q_MAG_SIZE,err); 


	//cleanup
	Py_DECREF(Nnz_array);
	Py_DECREF(par_nz_array);

	return model;

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

	model=parse_model(args,err);
	
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