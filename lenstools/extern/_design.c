/*Python wrapper module for parameter space sampling;
the method used is the same as in 
http://dan.iel.fm/posts/python-c-extensions/ 

The module is called _design and it defines the methods below (see docstrings)
*/

#include <gsl/gsl_matrix.h>

#include <Python.h>
#include <numpy/arrayobject.h>

#include "lenstoolsPy3.h"
#include "design.h"

#ifndef IS_PY3K
static struct module_state _state;
#endif

//Python module docstrings 
static char module_docstring[] = "This module provides a python interface for sampling a N-dimensional parameter space";
static char diagonalCost_docstring[] = "Compute the cost function for a diagonal design";
static char cost_docstring[] = "Compute the cost function for a given";
static char sample_docstring[] = "Sample the N-dimensional parameter space and record cost function changes";

//Method declarations
static PyObject *_design_diagonalCost(PyObject *self,PyObject *args);
static PyObject *_design_cost(PyObject *self,PyObject *args);
static PyObject *_design_sample(PyObject *self,PyObject *args);

//_design method definitions
static PyMethodDef module_methods[] = {

	{"diagonalCost",_design_diagonalCost,METH_VARARGS,diagonalCost_docstring},
	{"cost",_design_cost,METH_VARARGS,cost_docstring},
	{"sample",_design_sample,METH_VARARGS,sample_docstring},
	{NULL,NULL,0,NULL}

} ;

//_design constructor

#ifdef IS_PY3K

static struct PyModuleDef moduledef = {

	PyModuleDef_HEAD_INIT,
	"_design",
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
PyInit__design(void)

#else

#define INITERROR return

void
init_design(void)
#endif

{
#ifdef IS_PY3K
	PyObject *m = PyModule_Create(&moduledef);
#else
	PyObject *m = Py_InitModule3("_design",module_methods,module_docstring);
#endif

	if(m==NULL)
		INITERROR;
	struct module_state *st = GETSTATE(m);

	st->error = PyErr_NewException("_design.Error", NULL, NULL);
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

//////////////////////////////////////////////////
/*Function implementations using backend C code*/
/////////////////////////////////////////////////

//diagonalCost() implementation
static PyObject *_design_diagonalCost(PyObject *self,PyObject *args){

	int Npoints;
	double lambda,costValue;

	//Interpret arguments
	if(!PyArg_ParseTuple(args,"id",&Npoints,&lambda)){
		return NULL;
	}

	//Compute cost function for a diagonal design
	costValue = diagonalCost(Npoints,lambda);

	//Throw error if cost is not positive
	if(costValue<=0.0){
		PyErr_SetString(PyExc_RuntimeError,"Cost function is not positive!");
		return NULL;
	}

	//Build return value and return
	PyObject *ret = Py_BuildValue("d",costValue);
	return ret;

}

//cost() implementation
static PyObject *_design_cost(PyObject *self,PyObject *args){

	int i,j;
	double p,lambda,costValue;
	PyObject *data_obj;

	/*Parse the input tuple*/
	if(!PyArg_ParseTuple(args,"Odd",&data_obj,&p,&lambda)){
		return NULL;
	}

	/*Interpret the points in the design as a numpy array*/
	PyObject *data_array = PyArray_FROM_OTF(data_obj,NPY_DOUBLE,NPY_IN_ARRAY);
	if(data_array==NULL){
		return NULL;
	}

	/*Get a pointer to the array data*/
	double *data = (double *)PyArray_DATA(data_array);

	/*squeeze the dimensions of the parameter space and the number of points*/
	int Npoints = (int)PyArray_DIM(data_array,0);
	int Ndim = (int)PyArray_DIM(data_array,1);

	/*Wrap a gsl matrix object around the data_array*/
	gsl_matrix *m = gsl_matrix_alloc(Npoints,Ndim);
	
	/*Copy the elements into it*/
	for(i=0;i<Npoints;i++){
		for(j=0;j<Ndim;j++){
			gsl_matrix_set(m,i,j,data[i*Ndim+j]);
		}
	}

	/*Call the C backend to compute the cost function*/
	costValue = cost(m,Npoints,Ndim,p,lambda);

	/*Release the gsl resources*/
	gsl_matrix_free(m);

	/*Throw error if cost is not positive*/
	if(costValue<=0.0){
		PyErr_SetString(PyExc_RuntimeError,"Cost function is not positive!");
		Py_DECREF(data_array);
		return NULL;
	}

	/*Build the return value and return it*/
	PyObject *ret = Py_BuildValue("d",costValue);
	Py_DECREF(data_array);

	return ret;

}

//sampler() implementation
static PyObject *_design_sample(PyObject *self,PyObject *args){

	int i,j,maxIterations,seed;
	double p,lambda;
	PyObject *data_obj,*cost_obj;

	/*Parse the input tuple*/
	if(!PyArg_ParseTuple(args,"OddiiO",&data_obj,&p,&lambda,&maxIterations,&seed,&cost_obj)){
		return NULL;
	}

	/*Interpret the data_obj as an array*/
	PyObject *data_array = PyArray_FROM_OTF(data_obj,NPY_DOUBLE,NPY_IN_ARRAY);
	PyObject *cost_array = PyArray_FROM_OTF(cost_obj,NPY_DOUBLE,NPY_IN_ARRAY);
	
	if(data_array==NULL || cost_array==NULL){
		
		Py_XDECREF(data_array);
		Py_XDECREF(cost_array);

		return NULL;
	}

	/*Get a pointer to the array data*/
	double *data = (double *)PyArray_DATA(data_array);

	/*Get a pointer to the cost values*/
	double *cost_values = (double *)PyArray_DATA(cost_array);

	/*squeeze the dimensions of the parameter space and the number of points*/
	int Npoints = (int)PyArray_DIM(data_array,0);
	int Ndim = (int)PyArray_DIM(data_array,1);

	/*Wrap a gsl matrix object around the data_array*/
	gsl_matrix *m = gsl_matrix_alloc(Npoints,Ndim);
	
	/*Copy the elements into it*/
	for(i=0;i<Npoints;i++){
		for(j=0;j<Ndim;j++){
			gsl_matrix_set(m,i,j,data[i*Ndim+j]);
		}
	}

	/*Spread the points in the parameter space looking for the cost function minimum*/
	double deltaPerc = sample(Npoints,Ndim,p,lambda,seed,maxIterations,m,cost_values);

	/*Copy the points positions from the gsl matrix to the data array*/
	for(i=0;i<Npoints;i++){
		for(j=0;j<Ndim;j++){
			data[i*Ndim+j] = gsl_matrix_get(m,i,j);
		}
	}

	/*Release the gsl resources*/
	gsl_matrix_free(m);

	/*Release the other resources*/
	Py_DECREF(data_array);
	Py_DECREF(cost_array);

	/*Build the return value*/
	PyObject *ret = Py_BuildValue("d",deltaPerc);
	return ret;

}