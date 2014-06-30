/*Python wrapper module for peak count code interpolation code;
the method used is the same as in 
http://dan.iel.fm/posts/python-c-extensions/ 

The module is called _peaks and it defines a method called
count()
*/

#include <Python.h>
#include <numpy/arrayobject.h>

#include "peaks.h"

//Python module docstrings 
static char module_docstring[] = "This module provides a python interface for counting peaks in a 2D map";
static char count_docstring[] = "Calculate the peak counts in a 2D map";

//count() method declaration
static PyObject *_peaks_count(PyObject *self,PyObject *args);

//_peaks method definitions
static PyMethodDef module_methods[] = {

	{"count",_peaks_count,METH_VARARGS,count_docstring},
	{NULL,NULL,0,NULL}

} ;

//_peaks constructor
PyMODINIT_FUNC init_peaks(void){

	PyObject *m = Py_InitModule3("_peaks",module_methods,module_docstring);
	if(m==NULL) return;

	/*Load numpy functionality*/
	import_array();

}

//count() implementation
static PyObject *_peaks_count(PyObject *self,PyObject *args){

	PyObject *map_obj,*thresholds_obj; 
	double sigma;

	/*Parse the input tuple*/
	if(!PyArg_ParseTuple(args,"OOd",&map_obj,&thresholds_obj,&sigma)){ 
		return NULL;
	}

	/*Interpret the inputs as a numpy arrays*/
	PyObject *map_array = PyArray_FROM_OTF(map_obj,NPY_DOUBLE,NPY_IN_ARRAY);
	PyObject *thresholds_array = PyArray_FROM_OTF(thresholds_obj,NPY_DOUBLE,NPY_IN_ARRAY);

	/*Get the size of the map (in pixels)*/
	int Nside = (int)PyArray_DIM(map_array,0);
	/*Get the number of thresholds to output*/
	int Nthreshold = (int)PyArray_DIM(thresholds_array,0);

	/*Prepare a new array object that will hold the peak counts*/
	npy_intp dims[] = {(npy_intp) Nthreshold - 1};
	PyObject *peaks_array = PyArray_SimpleNew(1,dims,NPY_DOUBLE);

	/*Throw exception if this failed*/
	if(peaks_array==NULL){
		
		Py_DECREF(map_array);
		Py_DECREF(thresholds_array);

		return NULL;
	}

	/*Call the underlying C function that counts the peaks*/
	peak_count((double *)PyArray_DATA(map_array),Nside,sigma,Nthreshold,(double *)PyArray_DATA(thresholds_array),(double *)PyArray_DATA(peaks_array));

	/*Clean up and return*/
	Py_DECREF(map_array);
	Py_DECREF(thresholds_array);

	return peaks_array;

}