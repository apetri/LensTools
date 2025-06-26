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

/* NumPy 2.x deprecated constants compatibility */
#ifndef NPY_IN_ARRAY
#define NPY_IN_ARRAY NPY_ARRAY_IN_ARRAY
#endif

#ifndef NPY_OUT_ARRAY
#define NPY_OUT_ARRAY NPY_ARRAY_OUT_ARRAY
#endif

/* For NumPy 2.x compatibility - handle PyArray function casting */
#define PYARRAY_DATA_CAST(arr, type) ((type *)PyArray_DATA((const PyArrayObject*)(arr)))
#define PYARRAY_DIM_CAST(arr, idim) PyArray_DIM((const PyArrayObject*)(arr), idim)
#define PYARRAY_SIZE_CAST(arr) PyArray_SIZE((const PyArrayObject*)(arr))

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