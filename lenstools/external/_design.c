/*Python wrapper module for parameter space sampling;
the method used is the same as in 
http://dan.iel.fm/posts/python-c-extensions/ 

The module is called _design and it defines the methods below (see docstrings)
*/

#include <gsl/gsl_matrix.h>

#include <Python.h>
#include <numpy/arrayobject.h>

#include "design.h"

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
PyMODINIT_FUNC init_design(void){

	PyObject *m = Py_InitModule3("_design",module_methods,module_docstring);
	if(m==NULL) return;

	/*Load numpy functionality*/
	import_array();

}

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

	return NULL;
}