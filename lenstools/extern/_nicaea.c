/*Python wrapper module for the NICAEA weak lensing code by M.Kilbinger;

The module is called _nicaea and it defines the methods below (see docstrings)
*/

#include <stdio.h>
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

//Useful methods for parsing Nicaea class attributes into cosmo_lens structs
static cosmo_lens *parse_model(PyObject *args, error **err);
static int translate(int Nobjects, char *string_dictionary[],char *string);

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
/*Implementation of translate*///
/////////////////////////////////

static int translate(int Nobjects,char *string_dictionary[],char *string){

	int i;
	char error_buf[256];

	for(i=0;i<Nobjects;i++){
		if(strcmp(string_dictionary[i],string)==0) return i;
	}

	sprintf(error_buf,"Setting %s not implemented",string);
	PyErr_SetString(PyExc_ValueError,error_buf);
	return -1;

}

/////////////////////////////////
/*Implementation of parse_model*/
/////////////////////////////////

static cosmo_lens *parse_model(PyObject *args, error **err){

	//Type translators
	char *distribution_strings[5] = {"ludo","jonben","ymmk","hist","single"};
	const nofz_t distribution_types[5] = {ludo,jonben,ymmk,hist,single};

	//Cosmological parameters
	double Om,Ode,w0,w1,H100,Omegab,Omeganu,Neff,si8,ns;
	int nzbins;
	int i,j;
	char *distr_type,error_buf[256];

	//Multipoles/angles, redshift distribution and other settings
	PyObject *spec_obj,*Nnz_obj,*nofz_obj,*par_nz_obj,*settings_dict,*extra_obj;

	//Parse the input tuple
	if(!PyArg_ParseTuple(args,"ddddddddddiOOOOOO",&Om,&Ode,&w0,&w1,&H100,&Omegab,&Omeganu,&Neff,&si8,&ns,&nzbins,&spec_obj,&Nnz_obj,&nofz_obj,&par_nz_obj,&settings_dict,&extra_obj)){
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
		
		if((j=translate(5,distribution_strings,distr_type))==-1){
			
			Py_DECREF(Nnz_array);
			Py_DECREF(par_nz_array);
			return NULL;
		
		} else{

			fprintf(stderr,"%d=distr_type[%d]\n",i,j);
			nofz[i]=distribution_types[j];
		}

	}

	//Set these to default for the moment
	/*TODO parse from settings dictionary*/

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

	//counters
	int l;

	//cosmological model handler
	cosmo_lens *model;
	
	//multipoles and power spectrum
	PyObject *ell_array,*power_spectrum_array;

	//NICAEA error handlers
	error *myerr=NULL,**err;
	err=&myerr;

	//Convert NICAEA errors to string
	char stringerr[4096];

	//Build a cosmo_lens instance parsing the input tuple
	model=parse_model(args,err);
	if(model==NULL){
		return NULL;
	}

	//Read in the multipoles
	ell_array=PyArray_FROM_OTF(PyTuple_GetItem(args,11),NPY_DOUBLE,NPY_IN_ARRAY);
	if(ell_array==NULL){
		free_parameters_lens(&model);
		return NULL;
	}

	int Nl=(int)PyArray_DIM(ell_array,0);
	double *ell=(double *)PyArray_DATA(ell_array);

	//Allocate numpy array for the power spectrum
	npy_intp ell_dims[] = {(npy_intp)Nl};
	power_spectrum_array = PyArray_ZEROS(1,ell_dims,NPY_DOUBLE,0);
	if(power_spectrum_array==NULL){
		Py_DECREF(ell_array);
		free_parameters_lens(&model);
		return NULL;
	}

	//Call NICAEA to measure the power spectrum
	double *power_spectrum=(double *)PyArray_DATA(power_spectrum_array);
	for(l=0;l<Nl;l++){
		
		power_spectrum[l]=Pshear(model,ell[l],0,0,err);
		
		if(isError(*err)){
			stringError(stringerr,*err);
			PyErr_SetString(PyExc_RuntimeError,stringerr);
			free_parameters_lens(&model);
			Py_DECREF(ell_array);
			Py_DECREF(power_spectrum_array);
			return NULL;
		}

	}
	
	//Computation succeeded, cleanup and return
	Py_DECREF(ell_array);
	return power_spectrum_array;

}