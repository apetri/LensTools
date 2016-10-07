/*Python wrapper module for gadget snapshot reading;
the method used is the same as in 
http://dan.iel.fm/posts/python-c-extensions/ 

The module is called _pixelize and it defines the methods below (see docstrings)
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
static char module_docstring[] = "This module provides a python interface for two dimensional pixelizations of catalogs";
static char grid2d_docstring[] = "Construct a 2D pixelization of a scalar quantity in a catalog";

//Method declarations
static PyObject *_pixelize_grid2d(PyObject *self,PyObject *args);

//_pixelize method definitions
static PyMethodDef module_methods[] = {

	{"grid2d",_pixelize_grid2d,METH_VARARGS,grid2d_docstring},
	{NULL,NULL,0,NULL}

} ;

//_pixelize constructor

#ifdef IS_PY3K

static struct PyModuleDef moduledef = {

	PyModuleDef_HEAD_INIT,
	"_pixelize",
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
PyInit__pixelize(void)

#else

#define INITERROR return

void
init_pixelize(void)
#endif

{
#ifdef IS_PY3K
	PyObject *m = PyModule_Create(&moduledef);
#else
	PyObject *m = Py_InitModule3("_pixelize",module_methods,module_docstring);
#endif

	if(m==NULL)
		INITERROR;
	struct module_state *st = GETSTATE(m);

	st->error = PyErr_NewException("_pixelize.Error", NULL, NULL);
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


//grid2d() implementation
static PyObject *_pixelize_grid2d(PyObject *self,PyObject *args){

	PyObject *x_obj,*y_obj,*s_obj,*map_obj;
	double map_size;
	int err;

	//parse input tuple
	if(!PyArg_ParseTuple(args,"OOOdO",&x_obj,&y_obj,&s_obj,&map_size,&map_obj)) return NULL;

	//interpret arrays
	PyObject *x_array = PyArray_FROM_OTF(x_obj,NPY_DOUBLE,NPY_IN_ARRAY);
	PyObject *y_array = PyArray_FROM_OTF(y_obj,NPY_DOUBLE,NPY_IN_ARRAY);
	PyObject *s_array = PyArray_FROM_OTF(s_obj,NPY_DOUBLE,NPY_IN_ARRAY);
	PyObject *map_array = PyArray_FROM_OTF(map_obj,NPY_DOUBLE,NPY_IN_ARRAY);

	//check if anything failed
	if(x_array==NULL || y_array==NULL || s_array==NULL || map_array==NULL){
		
		Py_XDECREF(x_array);
		Py_XDECREF(y_array);
		Py_XDECREF(s_array);
		Py_XDECREF(map_array);

		return NULL;
	}

	//get the number of objects in the catalog and the number of pixels
	int Nobjects = (int)PyArray_DIM(x_array,0);
	int Npixel = (int)PyArray_DIM(map_array,0);

	//get the data pointers
	double *x = (double *)PyArray_DATA(x_array);
	double *y = (double *)PyArray_DATA(y_array);
	double *s = (double *)PyArray_DATA(s_array);
	double *map = (double *)PyArray_DATA(map_array);

	//call the C backend for the gridding procedure
	err = grid2d(x,y,s,map,Nobjects,Npixel,map_size);

	if(err){

		if(err==1) PyErr_SetString(PyExc_MemoryError,"A call to malloc failed");
		
		Py_DECREF(x_array);
		Py_DECREF(y_array);
		Py_DECREF(s_array);
		Py_DECREF(map_array);

		return NULL;

	}

	//cleanup
	Py_DECREF(x_array);
	Py_DECREF(y_array);
	Py_DECREF(s_array);
	Py_DECREF(map_array);

	//return None
	Py_RETURN_NONE;


}