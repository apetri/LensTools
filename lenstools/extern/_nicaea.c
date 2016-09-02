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

#include "lenstoolsPy3.h"

#include "errorlist.h"
#include "cosmo.h"
#include "lensing.h"
#include "nofz.h"

#include "nicaea_wrappers.h"

#define SPEC_IN_TUPLE 11

#ifndef IS_PY3K
static struct module_state _state;
#endif

//Python module docstrings 
static char module_docstring[] = "This module provides a python interface to the NICAEA computations";
static char shearPowerSpectrum_docstring[] = "Compute the shear power spectrum";
static char shear2pt_docstring[] = "Compute the shear correlation function";

//Useful methods for parsing Nicaea class attributes into cosmo_lens structs
static PyObject *extra_args(PyObject *args);
static int translate(int Nobjects, char *string_dictionary[],char *string);
static cosmo_lens *parse_model(PyObject *args, error **err);

//output memory allocator
static PyObject *alloc_output(PyObject *spec,cosmo_lens *model);

//NICAEA method wrapper
static PyObject *_nicaea_Wrapper(PyObject *args,double (*nicaea_method)(cosmo_lens*,double,int,int,error**));

//Method declarations
static PyObject *_nicaea_shearPowerSpectrum(PyObject *self,PyObject *args);
static PyObject *_nicaea_shear2pt(PyObject *self,PyObject *args);

//_nicaea method definitions
static PyMethodDef module_methods[] = {

	{"shearPowerSpectrum",_nicaea_shearPowerSpectrum,METH_VARARGS,shearPowerSpectrum_docstring},
	{"shear2pt",_nicaea_shear2pt,METH_VARARGS,shear2pt_docstring},
	{NULL,NULL,0,NULL}

} ;

//_nicaea constructor

#ifdef IS_PY3K

static struct PyModuleDef moduledef = {

	PyModuleDef_HEAD_INIT,
	"_nicaea",
	module_docstring,
	sizeof(struct module_state),
	module_methods,
	NULL,
	myextension_trasverse,
	myextension_clear,
	NULL

};

#define INITERROR return NULL

PyMODINIT_FUNC
PyInit__nicaea(void)

#else

#define INITERROR return

void
init_nicaea(void)
#endif

{
#ifdef IS_PY3K
	PyObject *m = PyModule_Create(&moduledef);
#else
	PyObject *m = Py_InitModule3("_nicaea",module_methods,module_docstring);
#endif

	if(m==NULL)
		INITERROR;
	struct module_state *st = GETSTATE(m);

	st->error = PyErr_NewException("_nicaea.Error", NULL, NULL);
    if (st->error == NULL) {
        Py_DECREF(m);
        INITERROR;
    }

	/*Load numpy functionality*/
	import_array();

	/*Return*/
#ifdef IS_PY3K
	return m;
#endif

}

//////////////////////////////////////////////
//////////////////////////////////////////////
//////////////////////////////////////////////

//////////////////////////////////
/*Implementation of extra_args*///
//////////////////////////////////
static PyObject *extra_args(PyObject *args){

	int args_size=(int)PyTuple_Size(args);
	return PyTuple_GetItem(args,args_size-1);
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
	char *distribution_strings[Nnofz_t] = {"ludo", "jonben", "ymmk", "ymmk0const", "hist", "single"};
	const nofz_t distribution_types[Nnofz_t] = {ludo, jonben, ymmk, ymmk0const, hist, single};

	char *nonlinear_strings[Nnonlinear_t] = {"linear", "pd96", "smith03", "smith03_de", "coyote10", "coyote13", "halodm", "smith03_revised"}; 
	const nonlinear_t nonlinear_types[Nnonlinear_t] = {linear, pd96, smith03, smith03_de, coyote10, coyote13, halodm, smith03_revised};

	char *transfer_strings[Ntransfer_t] = {"bbks", "eisenhu", "eisenhu_osc", "be84"};
	const transfer_t transfer_types[Ntransfer_t] = {bbks, eisenhu, eisenhu_osc, be84}; 

	char *growth_strings[Ngrowth_t] = {"heath", "growth_de"};
	const growth_t growth_types[Ngrowth_t] = {heath, growth_de};

	char *de_param_strings[Nde_param_t] = {"jassal", "linder", "earlyDE", "poly_DE"};
	const de_param_t de_param_types[Nde_param_t] = {jassal, linder, earlyDE, poly_DE};

	char *norm_strings[2] = {"norm_s8", "norm_as"};
	const norm_t norm_types[2] = {norm_s8, norm_as};

	char *tomo_strings[Ntomo_t] = {"tomo_all", "tomo_auto_only", "tomo_cross_only"};
	const tomo_t tomo_types[Ntomo_t] = {tomo_all, tomo_auto_only, tomo_cross_only};

	char *reduced_strings[Nreduced_t] = {"none", "reduced_K10"};
	const reduced_t reduced_types[Nreduced_t] = {reduced_none, reduced_K10};

	//Intrinsic alignment interface not implemented yet (defaults to none)
	ia_t IA=ia_none;
	ia_terms_t IA_TERMS=ia_undef;
	double A_IA=0.0;

	//Cosmological parameters
	double Om,Ode,w0,w1,H100,Omegab,Omeganu,Neff,si8,ns;
	int nzbins;
	int i,j;
	char *distr_type,*setting_string;

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
		
		if((j=translate(Nnofz_t,distribution_strings,distr_type))==-1){
			
			Py_DECREF(Nnz_array);
			Py_DECREF(par_nz_array);
			return NULL;
		
		} else{

			nofz[i]=distribution_types[j];
		}

	}

	/*Parse these computation settings from the settings dictionary*/

	//nonlinear_t
	nonlinear_t nonlinear_type;
	PyArg_Parse(PyDict_GetItemString(settings_dict,"snonlinear"),"s",&setting_string);
	if((j=translate(Nnonlinear_t,nonlinear_strings,setting_string))==-1){
		Py_DECREF(Nnz_array);
		Py_DECREF(par_nz_array);
		return NULL;
	} else{
		nonlinear_type=nonlinear_types[j];
	}

	//transfer_t
	transfer_t transfer_function;
	PyArg_Parse(PyDict_GetItemString(settings_dict,"stransfer"),"s",&setting_string);
	if((j=translate(Ntransfer_t,transfer_strings,setting_string))==-1){
		Py_DECREF(Nnz_array);
		Py_DECREF(par_nz_array);
		return NULL;
	} else{
		transfer_function=transfer_types[j];
	}

	//growth_t
	growth_t growth;
	PyArg_Parse(PyDict_GetItemString(settings_dict,"sgrowth"),"s",&setting_string);
	if((j=translate(Ngrowth_t,growth_strings,setting_string))==-1){
		Py_DECREF(Nnz_array);
		Py_DECREF(par_nz_array);
		return NULL;
	} else{
		growth=growth_types[j];
	}


	//de_param_t
	de_param_t dark_energy;
	PyArg_Parse(PyDict_GetItemString(settings_dict,"sde_param"),"s",&setting_string);
	if((j=translate(Nde_param_t,de_param_strings,setting_string))==-1){
		Py_DECREF(Nnz_array);
		Py_DECREF(par_nz_array);
		return NULL;
	} else{
		dark_energy=de_param_types[j];
	}

	//norm_t
	norm_t norm_mode;
	PyArg_Parse(PyDict_GetItemString(settings_dict,"normmode"),"s",&setting_string);
	if((j=translate(2,norm_strings,setting_string))==-1){
		Py_DECREF(Nnz_array);
		Py_DECREF(par_nz_array);
		return NULL;
	} else{
		norm_mode=norm_types[j];
	}

	//tomo_t
	tomo_t tomography;
	PyArg_Parse(PyDict_GetItemString(settings_dict,"stomo"),"s",&setting_string);
	if((j=translate(Ntomo_t,tomo_strings,setting_string))==-1){
		Py_DECREF(Nnz_array);
		Py_DECREF(par_nz_array);
		return NULL;
	} else{
		tomography=tomo_types[j];
	}

	//reduced_t
	reduced_t sreduced;
	PyArg_Parse(PyDict_GetItemString(settings_dict,"sreduced"),"s",&setting_string);
	if((j=translate(Nreduced_t,reduced_strings,setting_string))==-1){
		Py_DECREF(Nnz_array);
		Py_DECREF(par_nz_array);
		return NULL;
	} else{
		sreduced=reduced_types[j];
	}
	
	double Q_MAG_SIZE = PyFloat_AsDouble(PyDict_GetItemString(settings_dict,"q_mag_size"));

	//cosmo model object
	cosmo_lens *model=init_parameters_lens(Om,Ode,w0,w1,NULL,0,H100,Omegab,Omeganu,Neff,si8,ns,nzbins,Nnz,nofz,par_nz,nonlinear_type,transfer_function,growth,dark_energy,norm_mode,tomography,sreduced,Q_MAG_SIZE,IA,IA_TERMS,A_IA,err);

	//cleanup
	Py_DECREF(Nnz_array);
	Py_DECREF(par_nz_array);

	return model;

}

///////////////////////////////////
/*Implementation of alloc_output*/
//////////////////////////////////

//Allocate the appropriate python object for the function output (typically a Nl x Nbins array)
//spec is a numpy array with the multipoles/angles

static PyObject *alloc_output(PyObject *spec,cosmo_lens *model){

	int Ns = (int)PyArray_DIM(spec,0);
	int Nbins = model->redshift->Nzbin;
	tomo_t tomo = model->tomo;
	int Nz;

	//Number of redshift entries depends on tomography type (auto only, cross only or auto anc cross)
	if(tomo==tomo_auto_only){
		
		Nz=Nbins;
	
	} else if(tomo==tomo_cross_only){
		
		Nz=Nbins*(Nbins-1)/2;
		if(Nz==0){
			PyErr_SetString(PyExc_ValueError,"There is nothing to compute, you selected tomo_cross_only with only one redshift bin!!");
			return NULL;
		}

	} else if(tomo==tomo_all){

		Nz=Nbins*(Nbins+1)/2;

	} else return NULL;

	//Once we have the number of redshift components, we can allocate the array
	npy_intp output_dims[] = {(npy_intp)Ns,(npy_intp)Nz};
	PyObject *output_array = PyArray_ZEROS(2,output_dims,NPY_DOUBLE,0);

	//check if this succeeded
	if(output_array==NULL){
		return NULL;
	} else{
		return output_array;
	}


}

/////////////////////////////////////
/*Implementation of _nicaea_wrapper*/
/////////////////////////////////////

//shearPowerSpectrum() implementation
static PyObject *_nicaea_Wrapper(PyObject *args,double (*nicaea_method)(cosmo_lens*,double,int,int,error**)){

	//counters
	int l,i,j,b;

	//cosmological model handler
	cosmo_lens *model;
	
	//specifications (multipoles, angles) and output
	PyObject *spec_array,*output_array;

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
	spec_array=PyArray_FROM_OTF(PyTuple_GetItem(args,SPEC_IN_TUPLE),NPY_DOUBLE,NPY_IN_ARRAY);
	if(spec_array==NULL){
		free_parameters_lens(&model);
		return NULL;
	}

	//allocate the numpy array which will hold the output
	output_array = alloc_output(spec_array,model);
	if(output_array==NULL){
		Py_DECREF(spec_array);
		free_parameters_lens(&model);
		return NULL;
	}

	int Nl=(int)PyArray_DIM(output_array,0);
	int Nztot=(int)PyArray_DIM(output_array,1);
	int Nzbin=model->redshift->Nzbin;
	tomo_t tomo=model->tomo;
	
	//Data pointer to specification and output array
	double *spec=(double *)PyArray_DATA(spec_array);
	double *output=(double *)PyArray_DATA(output_array);

	//Call the NICAEA method

	//cycle over specification (angles or multipoles)
	for(l=0;l<Nl;l++){
		
		//cycle over redshift bins
		if(tomo==tomo_auto_only){
			
			//diagonal only
			for(i=0;i<Nzbin;i++){
				output[l*Nzbin+i]=nicaea_method(model,spec[l],i,i,err);
				if(isError(*err)) break;
			} 
		} else if(tomo==tomo_cross_only){

			b=0;
			//cross (off diagonal)
			for(i=0;i<Nzbin;i++){
				for(j=i+1;j<Nzbin;j++){

					output[l*Nztot+b]=nicaea_method(model,spec[l],i,j,err);
					b+=1;
					if(isError(*err)) break;
				}
				if(isError(*err)) break;
			}


		} else{

			b=0;
			//all (off and on diagonal)
			for(i=0;i<Nzbin;i++){
				for(j=i;j<Nzbin;j++){

					output[l*Nztot+b]=nicaea_method(model,spec[l],i,j,err);
					b+=1;
					if(isError(*err)) break;
				}
				if(isError(*err)) break;
			}

		}
		
		
		if(isError(*err)){
			stringError(stringerr,*err);
			PyErr_SetString(PyExc_RuntimeError,stringerr);
			free_parameters_lens(&model);
			Py_DECREF(spec_array);
			Py_DECREF(output_array);
			return NULL;
		}

	}
	
	//Computation succeeded, cleanup and return
	Py_DECREF(spec_array);
	return output_array;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////
/*Function implementations using backend C code*/
/////////////////////////////////////////////////
static PyObject *_nicaea_shearPowerSpectrum(PyObject *self,PyObject *args){ 

	return _nicaea_Wrapper(args,Pshear);
}

static PyObject *_nicaea_shear2pt(PyObject *self,PyObject *args){

#ifdef IS_PY3K
	int pm=(int)PyLong_AsLong(extra_args(args));
#else
	int pm=(int)PyInt_AsLong(extra_args(args));
#endif

	if(pm==1){
		
		return _nicaea_Wrapper(args,xi_plus);
	
	} else if(pm==-1){
		
		return _nicaea_Wrapper(args,xi_minus);
	
	} else{
		
		PyErr_SetString(PyExc_ValueError,"Only +1 and -1 allowed for pm!");
		return NULL;
	}


}

