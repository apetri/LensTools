/*Python wrapper module for gadget snapshot reading;
the method used is the same as in 
http://dan.iel.fm/posts/python-c-extensions/ 

The module is called _nbody and it defines the methods below (see docstrings)
*/

#include <stdio.h>

#include <Python.h>
#include <numpy/arrayobject.h>

#include "lenstoolsPy3.h"
#include "grid.h"

#ifndef IS_PY3K
static struct module_state _state;
#endif

//Python module docstrings
static char module_docstring[] = "This module provides a python interface for operations on Nbody simulation snapshots";
static char grid3d_docstring[] = "Put the snapshot particles on a regularly spaced grid";
static char grid3d_nfw_docstring[] = "Put the snapshot particles on a regularly spaced grid, but give each particle a NFW profile";
static char adaptive_docstring[] = "Put the snapshot particles on a regularly spaced grid using adaptive smoothing";

//Useful
static PyObject *_apply_kernel2d(PyObject *args,double(*kernel)(double,double,double,double));
static PyObject *_apply_kernel3d(PyObject *args,double(*kernel)(double,double,double,double));

//Method declarations
static PyObject * _nbody_grid3d(PyObject *self,PyObject *args);
static PyObject *_nbody_grid3d_nfw(PyObject *self,PyObject *args);
static PyObject * _nbody_adaptive(PyObject *self,PyObject *args);

//_nbody method definitions
static PyMethodDef module_methods[] = {

	{"grid3d",_nbody_grid3d,METH_VARARGS,grid3d_docstring},
	{"grid3d_nfw",_nbody_grid3d_nfw,METH_VARARGS,grid3d_nfw_docstring},
	{"adaptive",_nbody_adaptive,METH_VARARGS,adaptive_docstring},
	{NULL,NULL,0,NULL}

} ;

//_nbody constructor

#ifdef IS_PY3K

static struct PyModuleDef moduledef = {

	PyModuleDef_HEAD_INIT,
	"_nbody",
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
PyInit__nbody(void)

#else

#define INITERROR return

void
init_nbody(void)
#endif

{
#ifdef IS_PY3K
	PyObject *m = PyModule_Create(&moduledef);
#else
	PyObject *m = Py_InitModule3("_nbody",module_methods,module_docstring);
#endif

	if(m==NULL)
		INITERROR;
	struct module_state *st = GETSTATE(m);

	st->error = PyErr_NewException("_nbody.Error", NULL, NULL);
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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

//apply_kernel2d() implementation
static PyObject *_apply_kernel2d(PyObject *args,double(*kernel)(double,double,double,double)){

	PyObject *positions_obj,*weights_obj,*rp_obj,*concentration_obj,*binning_obj,*projectAll;
	double center,*concentration;
	int direction0,direction1,normal;
	float *weights;

	//Parse argument tuple
	if(!PyArg_ParseTuple(args,"OOOOOdiiiO",&positions_obj,&weights_obj,&rp_obj,&concentration_obj,&binning_obj,&center,&direction0,&direction1,&normal,&projectAll)){
		return NULL;
	}

	//Parse arrays
	PyObject *positions_array = PyArray_FROM_OTF(positions_obj,NPY_FLOAT32,NPY_IN_ARRAY);
	PyObject *rp_array = PyArray_FROM_OTF(rp_obj,NPY_DOUBLE,NPY_IN_ARRAY);
	PyObject *binning0_array = PyArray_FROM_OTF(PyList_GetItem(binning_obj,0),NPY_DOUBLE,NPY_IN_ARRAY);
	PyObject *binning1_array = PyArray_FROM_OTF(PyList_GetItem(binning_obj,1),NPY_DOUBLE,NPY_IN_ARRAY);

	//Check if anything went wrong
	if(positions_array==NULL || rp_array==NULL || binning0_array==NULL || binning1_array==NULL){
		
		Py_XDECREF(positions_array);
		Py_XDECREF(rp_array);
		Py_XDECREF(binning0_array);
		Py_XDECREF(binning1_array);

		return NULL;
	}

	//Parse the weights too if provided
	PyObject *weights_array,*concentration_array;

	if(weights_obj!=Py_None){

		weights_array = PyArray_FROM_OTF(weights_obj,NPY_FLOAT32,NPY_IN_ARRAY);
		if(weights_array==NULL){

			Py_DECREF(positions_array);
			Py_DECREF(rp_array);
			Py_DECREF(binning0_array);
			Py_DECREF(binning1_array);

			return NULL;

		}

		weights = (float *)PyArray_DATA(weights_array);

	} else{

		weights = NULL;

	}

	if(concentration_obj!=Py_None){

		concentration_array = PyArray_FROM_OTF(concentration_obj,NPY_DOUBLE,NPY_IN_ARRAY);
		if(concentration_array==NULL){

			Py_DECREF(positions_array);
			Py_DECREF(rp_array);
			Py_DECREF(binning0_array);
			Py_DECREF(binning1_array);
			if(weights) Py_DECREF(weights_array);

			return NULL;
		}

		concentration = (double *)PyArray_DATA(concentration_array);

	} else{

		concentration = NULL;

	}


	//Compute the number of particles
	int NumPart = (int)PyArray_DIM(positions_array,0);

	//Allocate space for lensing plane
	npy_intp dims[] =  {PyArray_DIM(binning0_array,0)-1,PyArray_DIM(binning1_array,0)-1};
	int size0 = (int)dims[0];
	int size1 = (int)dims[1];
	
	PyObject *lensingPlane_array = PyArray_ZEROS(2,dims,NPY_DOUBLE,0);
	
	if(lensingPlane_array==NULL){

		Py_DECREF(positions_array);
		Py_DECREF(rp_array);
		Py_DECREF(binning0_array);
		Py_DECREF(binning1_array);

		if(weights) Py_DECREF(weights_array);
		if(concentration) Py_DECREF(concentration_array);

		return NULL;

	}

	//Get data pointers
	float *positions = (float *)PyArray_DATA(positions_array);
	double *rp = (double *)PyArray_DATA(rp_array);
	double *binning0 = (double *)PyArray_DATA(binning0_array);
	double *binning1 = (double *)PyArray_DATA(binning1_array);
	double *lensingPlane = (double *)PyArray_DATA(lensingPlane_array);

	//Compute the adaptive smoothing using C backend
	adaptiveSmoothing(NumPart,positions,weights,rp,concentration,binning0,binning1,center,direction0,direction1,normal,size0,size1,PyObject_IsTrue(projectAll),lensingPlane,kernel);

	//Cleanup
	Py_DECREF(positions_array);
	Py_DECREF(rp_array);
	Py_DECREF(binning0_array);
	Py_DECREF(binning1_array);

	if(weights) Py_DECREF(weights_array);
	if(concentration) Py_DECREF(concentration_array);

	//Return
	return lensingPlane_array;

}

//apply_kernel3d() implementation
static PyObject *_apply_kernel3d(PyObject *args,double(*kernel)(double,double,double,double)){


	PyObject *positions_obj,*bins_obj,*weights_obj,*radius_obj,*concentration_obj;
	float *weights;
	double *radius,*concentration;

	//parse input tuple
	if(!PyArg_ParseTuple(args,"OOOOO",&positions_obj,&bins_obj,&weights_obj,&radius_obj,&concentration_obj)){
		return NULL;
	}

	//interpret parsed objects as arrays
	PyObject *positions_array = PyArray_FROM_OTF(positions_obj,NPY_FLOAT32,NPY_IN_ARRAY);
	PyObject *binsX_array = PyArray_FROM_OTF(PyTuple_GET_ITEM(bins_obj,0),NPY_DOUBLE,NPY_IN_ARRAY);
	PyObject *binsY_array = PyArray_FROM_OTF(PyTuple_GET_ITEM(bins_obj,1),NPY_DOUBLE,NPY_IN_ARRAY);
	PyObject *binsZ_array = PyArray_FROM_OTF(PyTuple_GET_ITEM(bins_obj,2),NPY_DOUBLE,NPY_IN_ARRAY);

	//check if anything failed
	if(positions_array==NULL || binsX_array==NULL || binsY_array==NULL || binsZ_array==NULL){
		
		Py_XDECREF(positions_array);
		Py_XDECREF(binsX_array);
		Py_XDECREF(binsY_array);
		Py_XDECREF(binsZ_array);

		return NULL;
	}

	//check if weights,radius,concentration are provided
	PyObject *weights_array,*radius_array,*concentration_array;
	
	if(weights_obj!=Py_None){
		
		weights_array = PyArray_FROM_OTF(weights_obj,NPY_FLOAT32,NPY_IN_ARRAY);
		
		if(weights_array==NULL){

			Py_DECREF(positions_array);
			Py_DECREF(binsX_array);
			Py_DECREF(binsY_array);
			Py_DECREF(binsZ_array);

			return NULL;

		}

		//Data pointer
		weights = (float *)PyArray_DATA(weights_array);

	} else{

		weights = NULL;
	}


	if(radius_obj!=Py_None){
		
		radius_array = PyArray_FROM_OTF(radius_obj,NPY_DOUBLE,NPY_IN_ARRAY);
		
		if(radius_array==NULL){

			Py_DECREF(positions_array);
			Py_DECREF(binsX_array);
			Py_DECREF(binsY_array);
			Py_DECREF(binsZ_array);
			if(weights) Py_DECREF(weights_array);

			return NULL;

		}

		//Data pointer
		radius = (double *)PyArray_DATA(radius_array);

	} else{

		radius = NULL;
	}


	if(concentration_obj!=Py_None){
		
		concentration_array = PyArray_FROM_OTF(concentration_obj,NPY_DOUBLE,NPY_IN_ARRAY);
		
		if(concentration_array==NULL){

			Py_DECREF(positions_array);
			Py_DECREF(binsX_array);
			Py_DECREF(binsY_array);
			Py_DECREF(binsZ_array);
			if(weights) Py_DECREF(weights_array);
			if(radius) Py_DECREF(radius_array);

			return NULL;

		}

		//Data pointer
		concentration = (double *)PyArray_DATA(concentration_array);

	} else{

		concentration = NULL;
	}

	//Get data pointers
	float *positions_data = (float *)PyArray_DATA(positions_array);
	double *binsX_data = (double *)PyArray_DATA(binsX_array);
	double *binsY_data = (double *)PyArray_DATA(binsY_array);
	double *binsZ_data = (double *)PyArray_DATA(binsZ_array);

	//Get info about the number of bins
	int NumPart = (int)PyArray_DIM(positions_array,0);
	int nx = (int)PyArray_DIM(binsX_array,0) - 1;
	int ny = (int)PyArray_DIM(binsY_array,0) - 1;
	int nz = (int)PyArray_DIM(binsZ_array,0) - 1;

	//Allocate the new array for the grid
	PyObject *grid_array;
		
	npy_intp gridDims[] = {(npy_intp) nx,(npy_intp) ny,(npy_intp) nz};
	grid_array = PyArray_ZEROS(3,gridDims,NPY_FLOAT32,0);

	if(grid_array==NULL){

		Py_DECREF(positions_array);
		Py_DECREF(binsX_array);
		Py_DECREF(binsY_array);
		Py_DECREF(binsZ_array);

		if(weights) Py_DECREF(weights_array);
		if(radius) Py_DECREF(radius_array);
		if(concentration) Py_DECREF(concentration_array);

		return NULL;

	}

	//Get a data pointer
	float *grid_data = (float *)PyArray_DATA(grid_array);

	//Snap the particles on the grid
	grid3d(positions_data,weights,radius,concentration,NumPart,binsX_data[0],binsY_data[0],binsZ_data[0],binsX_data[1] - binsX_data[0],binsY_data[1] - binsY_data[0],binsZ_data[1] - binsZ_data[0],nx,ny,nz,grid_data,kernel);

	//return the grid
	Py_DECREF(positions_array);
	Py_DECREF(binsX_array);
	Py_DECREF(binsY_array);
	Py_DECREF(binsZ_array);

	if(weights) Py_DECREF(weights_array);
	if(radius) Py_DECREF(radius_array);
	if(concentration) Py_DECREF(concentration_array);
	
	return grid_array;


}

//grid3d() implementation
static PyObject *_nbody_grid3d(PyObject *self,PyObject *args){

	return _apply_kernel3d(args,NULL);

}

//grid3d_nfw() implementation
static PyObject *_nbody_grid3d_nfw(PyObject *self,PyObject *args){

	return _apply_kernel3d(args,nfwKernel);

}


//adaptive() implementation
static PyObject * _nbody_adaptive(PyObject *self,PyObject *args){

	return _apply_kernel2d(args,quadraticKernel);

}