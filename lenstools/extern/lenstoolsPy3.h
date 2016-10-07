#ifndef __LENSTOOLS_PY3
#define __LENSTOOLS_PY3

#include "stdlib.h"
#include "Python.h"
#include "bytesobject.h"

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