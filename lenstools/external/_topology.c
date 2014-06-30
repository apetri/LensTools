/*Python wrapper module for peak count code interpolation code;
the method used is the same as in 
http://dan.iel.fm/posts/python-c-extensions/ 

The module is called _peaks and it defines a method called
count()
*/

#include <Python.h>
#include <numpy/arrayobject.h>

#include "peaks.h"
#include "differentials.h"

//Python module docstrings 
static char module_docstring[] = "This module provides a python interface for counting peaks in a 2D map";
static char peakCount_docstring[] = "Calculate the peak counts in a 2D map";
static char gradient_docstring[] = "Compute the gradient of a 2D image"; 

//method declarations
static PyObject *_topology_peakCount(PyObject *self,PyObject *args);
static PyObject *_topology_gradient(PyObject *self,PyObject *args);

//_topology method definitions
static PyMethodDef module_methods[] = {

	{"peakCount",_topology_peakCount,METH_VARARGS,peakCount_docstring},
	{"gradient",_topology_gradient,METH_VARARGS,gradient_docstring},
	{NULL,NULL,0,NULL}

} ;

//_topology constructor
PyMODINIT_FUNC init_topology(void){

	PyObject *m = Py_InitModule3("_topology",module_methods,module_docstring);
	if(m==NULL) return;

	/*Load numpy functionality*/
	import_array();

}

//peakCount() implementation
static PyObject *_topology_peakCount(PyObject *self,PyObject *args){

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

//gradient() implementation
static PyObject *_topology_gradient(PyObject *self,PyObject *args){

	PyObject *map_obj;

	/*Parse the input*/
	if(!PyArg_ParseTuple(args,"O",&map_obj)){ 
		return NULL;
	}

	/*Interpret the input as a numpy array*/
	PyObject *map_array = PyArray_FROM_OTF(map_obj,NPY_DOUBLE,NPY_IN_ARRAY);

	/*Get the size of the map (in pixels)*/
	int Nside = (int)PyArray_DIM(map_array,0);

	/*Prepare the new array objects that will hold the gradients along x and y*/
	npy_intp dims[] = {(npy_intp) Nside, (npy_intp) Nside};
	PyObject *gradient_x_array = PyArray_SimpleNew(2,dims,NPY_DOUBLE);
	PyObject *gradient_y_array = PyArray_SimpleNew(2,dims,NPY_DOUBLE);

	/*Throw exception if this failed*/
	if(gradient_x_array==NULL || gradient_y_array==NULL){
		
		Py_DECREF(map_array);
		return NULL;
	}

	/*Call the underlying C function that computes the gradient*/
	gradient_xy((double *)PyArray_DATA(map_array),(double *)PyArray_DATA(gradient_x_array),(double *)PyArray_DATA(gradient_y_array),Nside);

	/*Prepare a tuple container for the output*/
	PyObject *gradient_output = PyTuple_New(2);
	
	if(PyTuple_SetItem(gradient_output,0,gradient_x_array)){

		Py_DECREF(map_array);
		Py_DECREF(gradient_x_array);
		Py_DECREF(gradient_y_array);
		Py_DECREF(gradient_output);
		return NULL;

	}

	if(PyTuple_SetItem(gradient_output,1,gradient_y_array)){

		Py_DECREF(map_array);
		Py_DECREF(gradient_x_array);
		Py_DECREF(gradient_y_array);
		Py_DECREF(gradient_output);
		return NULL;

	}

	/*Done, now return*/
	return gradient_output;

}