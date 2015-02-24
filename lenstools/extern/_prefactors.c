#include <Python.h>
#include <numpy/arrayobject.h>

#include "interface.h"

static char module_docstring[] = "This module computes the growth and velocity prefactors for generating Gadget ICs";
static char prefactors_docstring[] = "Numerical values of the prefactors";

//Method declarations
static PyObject *_prefactors_prefactors(PyObject *self,PyObject *args);

static PyMethodDef module_methods[] = {

	{"prefactors",_prefactors_prefactors,METH_VARARGS,prefactors_docstring},
	{NULL,NULL,0,NULL}

} ;


//_prefactors constructor
PyMODINIT_FUNC init_prefactors(void){

	PyObject *m = Py_InitModule3("_prefactors",module_methods,module_docstring);
	if(m==NULL) return;

	/*Load numpy functionality*/
	import_array();

}

//Method implementations
static PyObject *_prefactors_prefactors(PyObject *self,PyObject *args){



	Py_RETURN_NONE;

}