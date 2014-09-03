/*Python wrapper module for parameter space sampling;
the method used is the same as in 
http://dan.iel.fm/posts/python-c-extensions/ 

The module is called _design and it defines the methods below (see docstrings)
*/

#include <Python.h>
#include <numpy/arrayobject.h>

#include "design.h"

//Python module docstrings 
static char module_docstring[] = "This module provides a python interface for sampling a N-dimensional parameter space";
static char diagonalCost_docstring[] = "Compute the cost function for a diagonal design";
static char sampler_docstring[] = "Sample the N-dimensional parameter space and record cost function changes";

//Method declarations
static PyObject *_design_diagonalCost(PyObject *self,PyObject *args);
static PyObject *_design_sampler(PyObject *self,PyObject *args);

//_design method definitions
static PyMethodDef module_methods[] = {

	{"diagonalCost",_design_diagonalCost,METH_VARARGS,diagonalCost_docstring},
	{"sampler",_design_sampler,METH_VARARGS,sampler_docstring},
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

//sampler() implementation
static PyObject *_design_sampler(PyObject *self,PyObject *args){

	return NULL;
}