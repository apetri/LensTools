#ifndef __LENSTOOLS_PY3
#define __LENSTOOLS_PY3

#include "stdlib.h"
#include "Python.h"
#include "bytesobject.h"
#include <numpy/arrayobject.h>

/* NumPy 2.x compatibility macros */
#ifndef NPY_ARRAY_WRITEABLE  
#define NPY_ARRAY_WRITEABLE NPY_WRITEABLE
#endif

/* For NumPy 2.x compatibility - handle PyArray_DATA casting */
#if NPY_API_VERSION >= 0x00000018  /* NumPy 2.0+ */
#define PYARRAY_DATA_CAST(arr, type) ((type *)PyArray_DATA((PyArrayObject*)(arr)))
#else
#define PYARRAY_DATA_CAST(arr, type) ((type *)PyArray_DATA(arr))
#endif

struct module_state {
	PyObject *error;
};

#if PY_MAJOR_VERSION >= 3

#define IS_PY3K
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))

#else
#define GETSTATE(m) (&_state)

#endif

#ifdef IS_PY3K

static int myextension_trasverse(PyObject *m, visitproc visit, void *arg){
	Py_VISIT(GETSTATE(m)->error);
	return 0;
}

static int myextension_clear(PyObject *m){
	Py_CLEAR(GETSTATE(m)->error);
	return 0;
}

#endif

#endif