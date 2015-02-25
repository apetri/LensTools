#include <Python.h>
#include <numpy/arrayobject.h>

#include "interface.h"

static char module_docstring[] = "This module computes the growth and velocity prefactors for generating Gadget ICs";
static char velocity_docstring[] = "Numerical values of the prefactors";

//Method declarations
static PyObject *_prefactors_velocity(PyObject *self,PyObject *args);

static PyMethodDef module_methods[] = {

	{"velocity",_prefactors_velocity,METH_VARARGS,velocity_docstring},
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
static PyObject *_prefactors_velocity(PyObject *self,PyObject *args){

	double z,Om,Ol,Ov,w0,wa,h,d2,d2minus,delz,zmaxact,zminact,prefactor;
	int iwmode,err;

	//Interpret input tuple
	if(!PyArg_ParseTuple(args,"ddddddddddddi",&z,&Om,&Ol,&Ov,&w0,&wa,&h,&d2,&d2minus,&delz,&zmaxact,&zminact,&iwmode)){
		PyErr_SetString(PyExc_ValueError,"Badly formatted input tuple in _prefactors_velocity");
		return NULL;
	}


	//Compute prefactor with C subroutine
	prefactor = fvel_Interface(0.0,z,Om,Ol,Ov,w0,wa,h,d2,d2minus,delz,zmaxact,zminact,iwmode,&err);

	if(err==1) PyErr_SetString(PyExc_RuntimeError,"Step size too small in odeint");
	if(err==2) PyErr_SetString(PyExc_RuntimeError,"Too many steps in routine odeint");

	//Return the prefactor to the python user
	return Py_BuildValue("d",prefactor);

}