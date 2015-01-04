#include <Python.h>
#include <arrayobject.h>

#include "countpairs.h"
#include "countpairs_rp_pi.h"
#include "countpairs_wp.h"

#ifndef MAXLEN
#define MAXLEN 1000
#endif


//Docstrings for the methods
static char module_docstring[]             =	"This module provides an interface for calculating correlation functions using C.";
static char countpairs_docstring[]         =	"Calculate the 3-D \\xi auto/cross-correlation function given two sets of X1/Y1/Z1 and X2/Y2/Z2 arrays.";
static char countpairs_rp_pi_docstring[]   =	"Calculate the 2-D DD(rp,pi) auto/cross-correlation function given two sets of X1/Y1/Z1 and X2/Y2/Z2 arrays.";
static char countpairs_wp_docstring[]      =	"Calculate the projected auto-correlation function wp (assumes PERIODIC) given one set of X1/Y1/Z1 arrays.";

/* function proto-type*/
static PyObject *countpairs_countpairs(PyObject *self, PyObject *args);
static PyObject *countpairs_countpairs_rp_pi(PyObject *self, PyObject *args);
static PyObject *countpairs_countpairs_wp(PyObject *self, PyObject *args);
PyMODINIT_FUNC init_countpairs(void);

static PyMethodDef module_methods[] = {
	{"countpairs"       ,countpairs_countpairs       ,METH_VARARGS,countpairs_docstring},
	{"countpairs_rp_pi" ,countpairs_countpairs_rp_pi ,METH_VARARGS,countpairs_rp_pi_docstring},
	{"countpairs_wp"    ,countpairs_countpairs_wp    ,METH_VARARGS,countpairs_wp_docstring},
	{NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC init_countpairs(void)
{
	PyObject *m = Py_InitModule3("_countpairs", module_methods, module_docstring);
	if (m == NULL)
		return;

	/* Load `numpy` functionality. */
	import_array();
}



static PyObject *countpairs_countpairs(PyObject *self, PyObject *args)
{
	(void) self;//to suppress the unused variable warning. Terrible hack
	PyObject *x1_obj, *y1_obj, *z1_obj, *x2_obj,*y2_obj,*z2_obj;
	struct timeval t0,t1;
	gettimeofday(&t0,NULL);
	int autocorr=0;
	int nthreads=4;
	char binfile[MAXLEN];
	
	if (!PyArg_ParseTuple(args, "iisOOOOOO",&autocorr,&nthreads,&binfile,&x1_obj,&y1_obj,&z1_obj,&x2_obj,&y2_obj,&z2_obj))
		return NULL;

	/* Interpret the input objects as numpy arrays. */
#ifdef DOUBLE_PREC	
	PyObject *x1_array = PyArray_FROM_OTF(x1_obj, NPY_DOUBLE, NPY_IN_ARRAY);
	PyObject *y1_array = PyArray_FROM_OTF(y1_obj, NPY_DOUBLE, NPY_IN_ARRAY);
	PyObject *z1_array = PyArray_FROM_OTF(z1_obj, NPY_DOUBLE,	NPY_IN_ARRAY);
#else
	PyObject *x1_array = PyArray_FROM_OTF(x1_obj, NPY_FLOAT, NPY_IN_ARRAY);
	PyObject *y1_array = PyArray_FROM_OTF(y1_obj, NPY_FLOAT, NPY_IN_ARRAY);
	PyObject *z1_array = PyArray_FROM_OTF(z1_obj, NPY_FLOAT, NPY_IN_ARRAY);
#endif
	
	if (x1_array == NULL || y1_array == NULL || z1_array == NULL) {
		Py_XDECREF(x1_array);
		Py_XDECREF(y1_array);
		Py_XDECREF(z1_array);
		return NULL;
	}


	/* Interpret the input objects as numpy arrays. */
#ifdef DOUBLE_PREC	
	PyObject *x2_array = PyArray_FROM_OTF(x2_obj, NPY_DOUBLE, NPY_IN_ARRAY);
	PyObject *y2_array = PyArray_FROM_OTF(y2_obj, NPY_DOUBLE, NPY_IN_ARRAY);
	PyObject *z2_array = PyArray_FROM_OTF(z2_obj, NPY_DOUBLE,	NPY_IN_ARRAY);
#else
	PyObject *x2_array = PyArray_FROM_OTF(x2_obj, NPY_FLOAT, NPY_IN_ARRAY);
	PyObject *y2_array = PyArray_FROM_OTF(y2_obj, NPY_FLOAT, NPY_IN_ARRAY);
	PyObject *z2_array = PyArray_FROM_OTF(z2_obj, NPY_FLOAT, NPY_IN_ARRAY);
#endif
	
	if (x2_array == NULL || y2_array == NULL || z2_array == NULL) {
		Py_XDECREF(x2_array);
		Py_XDECREF(y2_array);
		Py_XDECREF(z2_array);
		return NULL;
	}

	
  /* How many data points are there? */
	const int ND1 = (int)PyArray_DIM(x1_array, 0);
	const int ND2 = (int)PyArray_DIM(x2_array, 0);

	/* Get pointers to the data as C-types. */
	DOUBLE *x1 = (DOUBLE *)PyArray_DATA(x1_array);
	DOUBLE *y1 = (DOUBLE *)PyArray_DATA(y1_array);
	DOUBLE *z1 = (DOUBLE *)PyArray_DATA(z1_array);

	DOUBLE *x2 = (DOUBLE *)PyArray_DATA(x2_array);
	DOUBLE *y2 = (DOUBLE *)PyArray_DATA(y2_array);
	DOUBLE *z2 = (DOUBLE *)PyArray_DATA(z2_array);
		

	results_countpairs *results = countpairs(ND1,x1,y1,z1,
																					 ND2,x2,y2,z2,
#ifdef USE_OMP
																					 nthreads,
#endif
																					 autocorr,
																					 binfile);
		
	/* Clean up. */
	Py_DECREF(x1_array);Py_DECREF(y1_array);Py_DECREF(z1_array);
	Py_DECREF(x2_array);Py_DECREF(y2_array);Py_DECREF(z2_array);

	double value=0.0;
	gettimeofday(&t1,NULL);
  /* Build the output tuple */
	PyObject *ret = Py_BuildValue("d", value);

	free_results(&results);
	return ret;
}


static PyObject *countpairs_countpairs_rp_pi(PyObject *self, PyObject *args)
{
	(void) self;//to suppress the unused variable warning. Terrible hack
	PyObject *x1_obj, *y1_obj, *z1_obj, *x2_obj,*y2_obj,*z2_obj;
	struct timeval t0,t1;
	char binfile[MAXLEN];
	gettimeofday(&t0,NULL);
	int autocorr=0;
	int nthreads=4;
	double pimax;
	
	if (!PyArg_ParseTuple(args, "iidsOOOOOO",&autocorr,&nthreads,&pimax,&binfile,&x1_obj,&y1_obj,&z1_obj,&x2_obj,&y2_obj,&z2_obj))
		return NULL;

	/* Interpret the input objects as numpy arrays. */
#ifdef DOUBLE_PREC	
	PyObject *x1_array = PyArray_FROM_OTF(x1_obj, NPY_DOUBLE, NPY_IN_ARRAY);
	PyObject *y1_array = PyArray_FROM_OTF(y1_obj, NPY_DOUBLE, NPY_IN_ARRAY);
	PyObject *z1_array = PyArray_FROM_OTF(z1_obj, NPY_DOUBLE,	NPY_IN_ARRAY);
#else
	PyObject *x1_array = PyArray_FROM_OTF(x1_obj, NPY_FLOAT, NPY_IN_ARRAY);
	PyObject *y1_array = PyArray_FROM_OTF(y1_obj, NPY_FLOAT, NPY_IN_ARRAY);
	PyObject *z1_array = PyArray_FROM_OTF(z1_obj, NPY_FLOAT, NPY_IN_ARRAY);
#endif
	
	if (x1_array == NULL || y1_array == NULL || z1_array == NULL) {
		Py_XDECREF(x1_array);
		Py_XDECREF(y1_array);
		Py_XDECREF(z1_array);
		return NULL;
	}


	/* Interpret the input objects as numpy arrays. */
#ifdef DOUBLE_PREC	
	PyObject *x2_array = PyArray_FROM_OTF(x2_obj, NPY_DOUBLE, NPY_IN_ARRAY);
	PyObject *y2_array = PyArray_FROM_OTF(y2_obj, NPY_DOUBLE, NPY_IN_ARRAY);
	PyObject *z2_array = PyArray_FROM_OTF(z2_obj, NPY_DOUBLE,	NPY_IN_ARRAY);
#else
	PyObject *x2_array = PyArray_FROM_OTF(x2_obj, NPY_FLOAT, NPY_IN_ARRAY);
	PyObject *y2_array = PyArray_FROM_OTF(y2_obj, NPY_FLOAT, NPY_IN_ARRAY);
	PyObject *z2_array = PyArray_FROM_OTF(z2_obj, NPY_FLOAT, NPY_IN_ARRAY);
#endif
	
	if (x2_array == NULL || y2_array == NULL || z2_array == NULL) {
		Py_XDECREF(x2_array);
		Py_XDECREF(y2_array);
		Py_XDECREF(z2_array);
		return NULL;
	}

	
  /* How many data points are there? */
	const int ND1 = (int)PyArray_DIM(x1_array, 0);
	const int ND2 = (int)PyArray_DIM(x2_array, 0);

	/* Get pointers to the data as C-types. */
	DOUBLE *x1 = (DOUBLE *)PyArray_DATA(x1_array);
	DOUBLE *y1 = (DOUBLE *)PyArray_DATA(y1_array);
	DOUBLE *z1 = (DOUBLE *)PyArray_DATA(z1_array);

	DOUBLE *x2 = (DOUBLE *)PyArray_DATA(x2_array);
	DOUBLE *y2 = (DOUBLE *)PyArray_DATA(y2_array);
	DOUBLE *z2 = (DOUBLE *)PyArray_DATA(z2_array);
		

	results_countpairs_rp_pi *results = countpairs_rp_pi(ND1,x1,y1,z1,
																											 ND2,x2,y2,z2,
#ifdef USE_OMP
																											 nthreads,
#endif
																											 autocorr,
																											 binfile,
																											 pimax);
	
	/* Clean up. */
	Py_DECREF(x1_array);Py_DECREF(y1_array);Py_DECREF(z1_array);
	Py_DECREF(x2_array);Py_DECREF(y2_array);Py_DECREF(z2_array);

	double value=0.0;
	gettimeofday(&t1,NULL);
  /* Build the output tuple */
	PyObject *ret = Py_BuildValue("d", value);

	free_results_rp_pi(&results);
	return ret;
}

static PyObject *countpairs_countpairs_wp(PyObject *self, PyObject *args)
{
	(void) self;//to suppress the unused variable warning. Terrible hack
	PyObject *x1_obj, *y1_obj, *z1_obj;
	double boxsize,pimax;
	char binfile[MAXLEN];
	struct timeval t0,t1;
	gettimeofday(&t0,NULL);
	int nthreads=4;
	
	if (!PyArg_ParseTuple(args, "ddisOOO",&boxsize,&pimax,&nthreads,&binfile,&x1_obj,&y1_obj,&z1_obj))
		return NULL;

	/* Interpret the input objects as numpy arrays. */
#ifdef DOUBLE_PREC	
	PyObject *x1_array = PyArray_FROM_OTF(x1_obj, NPY_DOUBLE, NPY_IN_ARRAY);
	PyObject *y1_array = PyArray_FROM_OTF(y1_obj, NPY_DOUBLE, NPY_IN_ARRAY);
	PyObject *z1_array = PyArray_FROM_OTF(z1_obj, NPY_DOUBLE,	NPY_IN_ARRAY);
#else
	PyObject *x1_array = PyArray_FROM_OTF(x1_obj, NPY_FLOAT, NPY_IN_ARRAY);
	PyObject *y1_array = PyArray_FROM_OTF(y1_obj, NPY_FLOAT, NPY_IN_ARRAY);
	PyObject *z1_array = PyArray_FROM_OTF(z1_obj, NPY_FLOAT, NPY_IN_ARRAY);
#endif
	
	if (x1_array == NULL || y1_array == NULL || z1_array == NULL) {
		Py_XDECREF(x1_array);
		Py_XDECREF(y1_array);
		Py_XDECREF(z1_array);
		return NULL;
	}


  /* How many data points are there? */
	const int ND1 = (int)PyArray_DIM(x1_array, 0);

	/* Get pointers to the data as C-types. */
	DOUBLE *x1 = (DOUBLE *)PyArray_DATA(x1_array);
	DOUBLE *y1 = (DOUBLE *)PyArray_DATA(y1_array);
	DOUBLE *z1 = (DOUBLE *)PyArray_DATA(z1_array);
	results_countpairs_wp *results = countpairs_wp(ND1,x1,y1,z1,
																								 boxsize,
#ifdef USE_OMP
																								 nthreads,
#endif
																								 binfile,
																								 pimax);
	
	/* Clean up. */
	Py_DECREF(x1_array);Py_DECREF(y1_array);Py_DECREF(z1_array);

	double value=0.0;
	gettimeofday(&t1,NULL);
  /* Build the output tuple */
	PyObject *ret = Py_BuildValue("d", value);

	free_results_wp(&results);
	return ret;
}
