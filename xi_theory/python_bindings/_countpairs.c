/* File: _countpairs.c */
/*
		This file is a part of the Corrfunc package
		Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
		License: MIT LICENSE. See LICENSE file under the top-level
		directory at https://bitbucket.org/manodeep/corrfunc/
*/
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <arrayobject.h>

#include "countpairs.h"
#include "countpairs_rp_pi.h"
#include "countpairs_wp.h"
#include "countpairs_xi.h"

//for the vpf
#include "countspheres.h"


struct module_state {
	PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#define INITERROR return NULL

#else //python2 follows
#define GETSTATE(m) (&_state)
static struct module_state _state;
#define INITERROR return 
#endif 

static PyObject *error_out(PyObject *module) {
	struct module_state *st = GETSTATE(module);
	PyErr_SetString(st->error, "something bad happened");
	return NULL;
}


//Docstrings for the methods
static char module_docstring[]             =	"This module provides an interface for calculating correlation functions using C.";
static char countpairs_docstring[]         =	"Calculate the 3-D \\xi auto/cross-correlation function given two sets of X1/Y1/Z1 and X2/Y2/Z2 arrays.";
static char countpairs_rp_pi_docstring[]   =	"Calculate the 2-D DD(rp,pi) auto/cross-correlation function given two sets of X1/Y1/Z1 and X2/Y2/Z2 arrays.";
static char countpairs_wp_docstring[]      =	"Calculate the projected auto-correlation function wp (assumes PERIODIC) given one set of X1/Y1/Z1 arrays.";
static char countpairs_xi_docstring[]      =	"Calculate the 3-d auto-correlation function xi (assumes PERIODIC) given one set of X1/Y1/Z1 arrays.";
static char countspheres_vpf_docstring[]   =	"Calculate the counts-in-spheres given one set of X1/Y1/Z1 arrays.";

/* function proto-type*/
static PyObject *countpairs_countpairs(PyObject *self, PyObject *args);
static PyObject *countpairs_countpairs_rp_pi(PyObject *self, PyObject *args);
static PyObject *countpairs_countpairs_wp(PyObject *self, PyObject *args);
static PyObject *countpairs_countpairs_xi(PyObject *self, PyObject *args);
static PyObject *countpairs_countspheres_vpf(PyObject *self, PyObject *args);

PyMODINIT_FUNC init_countpairs(void);

static PyMethodDef module_methods[] = {
	{"error_out", (PyCFunction)error_out, METH_NOARGS, NULL},
	{"countpairs"       ,countpairs_countpairs       ,METH_VARARGS,countpairs_docstring},
	{"countpairs_rp_pi" ,countpairs_countpairs_rp_pi ,METH_VARARGS,countpairs_rp_pi_docstring},
	{"countpairs_wp"    ,countpairs_countpairs_wp    ,METH_VARARGS,countpairs_wp_docstring},
	{"countpairs_xi"    ,countpairs_countpairs_xi    ,METH_VARARGS,countpairs_xi_docstring},
	{"countspheres_vpf" ,countpairs_countspheres_vpf ,METH_VARARGS,countspheres_vpf_docstring},
	{NULL, NULL, 0, NULL}
};


#if PY_MAJOR_VERSION >= 3

static int _countpairs_traverse(PyObject *m, visitproc visit, void *arg) {
	Py_VISIT(GETSTATE(m)->error);
	return 0;
}

static int _countpairs_clear(PyObject *m) {
	Py_CLEAR(GETSTATE(m)->error);
	return 0;
}


static struct PyModuleDef moduledef = {
	PyModuleDef_HEAD_INIT,
	"_countpairs",
	NULL,
	sizeof(struct module_state),
	module_methods,
	NULL,
	_countpairs_traverse,
	_countpairs_clear,
  NULL
};


PyObject *PyInit_countpairs(void)

#else
PyMODINIT_FUNC init_countpairs(void)
#endif
{


#if PY_MAJOR_VERSION >= 3
	PyObject *module = PyModule_Create(&moduledef);
#else
	PyObject *module = Py_InitModule3("_countpairs", module_methods, module_docstring);
#endif
	
	if (module == NULL) {
		INITERROR;
	}

	struct module_state *st = GETSTATE(module);

	st->error = PyErr_NewException("_countpairs.Error", NULL, NULL);
	if (st->error == NULL) {
		Py_DECREF(module);
		INITERROR;
	}

	
	/* Load `numpy` functionality. */
	import_array();


#if PY_MAJOR_VERSION >= 3
	return module;
#endif

}


static PyObject *countpairs_countpairs(PyObject *self, PyObject *args)
{
	(void) self;//to suppress the unused variable warning. Terrible hack
	PyObject *x1_obj, *y1_obj, *z1_obj, *x2_obj,*y2_obj,*z2_obj;
	int autocorr=0;
	int nthreads=4;
	char *binfile;
	
	if (!PyArg_ParseTuple(args, "iisOOOOOO",&autocorr,&nthreads,&binfile,&x1_obj,&y1_obj,&z1_obj,&x2_obj,&y2_obj,&z2_obj))
		return NULL;

	/* Interpret the input objects as numpy arrays. */
#ifdef DOUBLE_PREC	
	PyObject *x1_array = PyArray_FROM_OTF(x1_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
	PyObject *y1_array = PyArray_FROM_OTF(y1_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
	PyObject *z1_array = PyArray_FROM_OTF(z1_obj, NPY_DOUBLE,	NPY_ARRAY_IN_ARRAY);
#else
	PyObject *x1_array = PyArray_FROM_OTF(x1_obj, NPY_FLOAT, NPY_ARRAY_IN_ARRAY);
	PyObject *y1_array = PyArray_FROM_OTF(y1_obj, NPY_FLOAT, NPY_ARRAY_IN_ARRAY);
	PyObject *z1_array = PyArray_FROM_OTF(z1_obj, NPY_FLOAT, NPY_ARRAY_IN_ARRAY);
#endif
	
	if (x1_array == NULL || y1_array == NULL || z1_array == NULL) {
		Py_XDECREF(x1_array);
		Py_XDECREF(y1_array);
		Py_XDECREF(z1_array);
		Py_RETURN_NONE;
	}


	/* Interpret the input objects as numpy arrays. */
#ifdef DOUBLE_PREC	
	PyObject *x2_array = PyArray_FROM_OTF(x2_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
	PyObject *y2_array = PyArray_FROM_OTF(y2_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
	PyObject *z2_array = PyArray_FROM_OTF(z2_obj, NPY_DOUBLE,	NPY_ARRAY_IN_ARRAY);
#else
	PyObject *x2_array = PyArray_FROM_OTF(x2_obj, NPY_FLOAT, NPY_ARRAY_IN_ARRAY);
	PyObject *y2_array = PyArray_FROM_OTF(y2_obj, NPY_FLOAT, NPY_ARRAY_IN_ARRAY);
	PyObject *z2_array = PyArray_FROM_OTF(z2_obj, NPY_FLOAT, NPY_ARRAY_IN_ARRAY);
#endif
	
	if (x2_array == NULL || y2_array == NULL || z2_array == NULL) {
		Py_XDECREF(x2_array);
		Py_XDECREF(y2_array);
		Py_XDECREF(z2_array);
		Py_RETURN_NONE;
	}

	
  /* How many data points are there? */
	const int64_t ND1 = (int64_t)PyArray_DIM(x1_array, 0);
	const int64_t ND2 = (int64_t)PyArray_DIM(x2_array, 0);

	/* Get pointers to the data as C-types. */
	DOUBLE *X1 = (DOUBLE *)PyArray_DATA(x1_array);
	DOUBLE *Y1 = (DOUBLE *)PyArray_DATA(y1_array);
	DOUBLE *Z1 = (DOUBLE *)PyArray_DATA(z1_array);

	DOUBLE *X2 = (DOUBLE *)PyArray_DATA(x2_array);
	DOUBLE *Y2 = (DOUBLE *)PyArray_DATA(y2_array);
	DOUBLE *Z2 = (DOUBLE *)PyArray_DATA(z2_array);
		

	results_countpairs *results = countpairs(ND1,X1,Y1,Z1,
																					 ND2,X2,Y2,Z2,
#ifdef USE_OMP
																					 nthreads,
#endif
																					 autocorr,
																					 binfile);
		
	/* Clean up. */
	Py_DECREF(x1_array);Py_DECREF(y1_array);Py_DECREF(z1_array);
	Py_DECREF(x2_array);Py_DECREF(y2_array);Py_DECREF(z2_array);

	/* Build the output list */
	PyObject *ret = PyList_New(0);
	DOUBLE rlow=results->rupp[0];
	for(int i=1;i<results->nbin;i++) {
		PyObject *item = NULL;
		const DOUBLE rpavg = results->rpavg[i];
		
#ifdef DOUBLE_PREC
	  item = Py_BuildValue("(dddk)", rlow,results->rupp[i],rpavg,results->npairs[i]);
#else
		item = Py_BuildValue("(fffk)", rlow,results->rupp[i],rpavg,results->npairs[i]);
#endif
		PyList_Append(ret, item);
		Py_XDECREF(item);
		rlow=results->rupp[i];
	}

	free_results(&results);
	return ret;
}


static PyObject *countpairs_countpairs_rp_pi(PyObject *self, PyObject *args)
{
	(void) self;//to suppress the unused variable warning. Terrible hack
	PyObject *x1_obj, *y1_obj, *z1_obj, *x2_obj,*y2_obj,*z2_obj;
	int autocorr=0;
	int nthreads=4;
	double pimax;
	char *binfile;
	
	if (!PyArg_ParseTuple(args, "iidsOOOOOO",&autocorr,&nthreads,&pimax,&binfile,&x1_obj,&y1_obj,&z1_obj,&x2_obj,&y2_obj,&z2_obj))
		return NULL;

	/* Interpret the input objects as numpy arrays. */
#ifdef DOUBLE_PREC	
	PyObject *x1_array = PyArray_FROM_OTF(x1_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
	PyObject *y1_array = PyArray_FROM_OTF(y1_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
	PyObject *z1_array = PyArray_FROM_OTF(z1_obj, NPY_DOUBLE,	NPY_ARRAY_IN_ARRAY);
#else
	PyObject *x1_array = PyArray_FROM_OTF(x1_obj, NPY_FLOAT, NPY_ARRAY_IN_ARRAY);
	PyObject *y1_array = PyArray_FROM_OTF(y1_obj, NPY_FLOAT, NPY_ARRAY_IN_ARRAY);
	PyObject *z1_array = PyArray_FROM_OTF(z1_obj, NPY_FLOAT, NPY_ARRAY_IN_ARRAY);
#endif
	
	if (x1_array == NULL || y1_array == NULL || z1_array == NULL) {
		Py_XDECREF(x1_array);
		Py_XDECREF(y1_array);
		Py_XDECREF(z1_array);
		return NULL;
	}


	/* Interpret the input objects as numpy arrays. */
#ifdef DOUBLE_PREC	
	PyObject *x2_array = PyArray_FROM_OTF(x2_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
	PyObject *y2_array = PyArray_FROM_OTF(y2_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
	PyObject *z2_array = PyArray_FROM_OTF(z2_obj, NPY_DOUBLE,	NPY_ARRAY_IN_ARRAY);
#else
	PyObject *x2_array = PyArray_FROM_OTF(x2_obj, NPY_FLOAT, NPY_ARRAY_IN_ARRAY);
	PyObject *y2_array = PyArray_FROM_OTF(y2_obj, NPY_FLOAT, NPY_ARRAY_IN_ARRAY);
	PyObject *z2_array = PyArray_FROM_OTF(z2_obj, NPY_FLOAT, NPY_ARRAY_IN_ARRAY);
#endif
	
	if (x2_array == NULL || y2_array == NULL || z2_array == NULL) {
		Py_XDECREF(x2_array);
		Py_XDECREF(y2_array);
		Py_XDECREF(z2_array);
		return NULL;
	}

	
  /* How many data points are there? */
	const int64_t ND1 = (int64_t)PyArray_DIM(x1_array, 0);
	const int64_t ND2 = (int64_t)PyArray_DIM(x2_array, 0);

	/* Get pointers to the data as C-types. */
	DOUBLE *X1 = (DOUBLE *)PyArray_DATA(x1_array);
	DOUBLE *Y1 = (DOUBLE *)PyArray_DATA(y1_array);
	DOUBLE *Z1 = (DOUBLE *)PyArray_DATA(z1_array);

	DOUBLE *X2 = (DOUBLE *)PyArray_DATA(x2_array);
	DOUBLE *Y2 = (DOUBLE *)PyArray_DATA(y2_array);
	DOUBLE *Z2 = (DOUBLE *)PyArray_DATA(z2_array);
		

	results_countpairs_rp_pi *results = countpairs_rp_pi(ND1,X1,Y1,Z1,
																											 ND2,X2,Y2,Z2,
#ifdef USE_OMP
																											 nthreads,
#endif
																											 autocorr,
																											 binfile,
																											 pimax);
	
	/* Clean up. */
	Py_DECREF(x1_array);Py_DECREF(y1_array);Py_DECREF(z1_array);
	Py_DECREF(x2_array);Py_DECREF(y2_array);Py_DECREF(z2_array);

	/* Build the output list */
	PyObject *ret = PyList_New(0);//create an empty list
	DOUBLE rlow=results->rupp[0];
	const DOUBLE dpi = pimax/(DOUBLE)results->npibin ;
	
	for(int i=1;i<results->nbin;i++) {
		for(int j=0;j<results->npibin;j++) {
			const int bin_index = i*(results->npibin + 1) + j;
			PyObject *item = NULL;
			const DOUBLE rpavg = results->rpavg[bin_index];
#ifdef DOUBLE_PREC
			item = Py_BuildValue("(ddddk)", rlow,results->rupp[i],rpavg,(j+1)*dpi,results->npairs[bin_index]);
#else
			item = Py_BuildValue("(ffffk)", rlow,results->rupp[i],rpavg,(j+1)*dpi,results->npairs[bin_index]);
#endif
			PyList_Append(ret, item);
			Py_XDECREF(item);
		}
		rlow=results->rupp[i];
	}
	free_results_rp_pi(&results);
	return ret;
}

static PyObject *countpairs_countpairs_wp(PyObject *self, PyObject *args)
{
	(void) self;//to suppress the unused variable warning. Terrible hack
	PyObject *x1_obj, *y1_obj, *z1_obj;
	double boxsize,pimax;
	int nthreads=4;
	char *binfile;
	
	if (!PyArg_ParseTuple(args, "ddisOOO",&boxsize,&pimax,&nthreads,&binfile,&x1_obj,&y1_obj,&z1_obj))
		return NULL;

	/* Interpret the input objects as numpy arrays. */
#ifdef DOUBLE_PREC	
	PyObject *x1_array = PyArray_FROM_OTF(x1_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
	PyObject *y1_array = PyArray_FROM_OTF(y1_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
	PyObject *z1_array = PyArray_FROM_OTF(z1_obj, NPY_DOUBLE,	NPY_ARRAY_IN_ARRAY);
#else
	PyObject *x1_array = PyArray_FROM_OTF(x1_obj, NPY_FLOAT, NPY_ARRAY_IN_ARRAY);
	PyObject *y1_array = PyArray_FROM_OTF(y1_obj, NPY_FLOAT, NPY_ARRAY_IN_ARRAY);
	PyObject *z1_array = PyArray_FROM_OTF(z1_obj, NPY_FLOAT, NPY_ARRAY_IN_ARRAY);
#endif
	
	if (x1_array == NULL || y1_array == NULL || z1_array == NULL) {
		Py_XDECREF(x1_array);
		Py_XDECREF(y1_array);
		Py_XDECREF(z1_array);
		return NULL;
	}


  /* How many data points are there? */
	const int64_t ND1 = (int64_t)PyArray_DIM(x1_array, 0);

	/* Get pointers to the data as C-types. */
	DOUBLE *X1 = (DOUBLE *)PyArray_DATA(x1_array);
	DOUBLE *Y1 = (DOUBLE *)PyArray_DATA(y1_array);
	DOUBLE *Z1 = (DOUBLE *)PyArray_DATA(z1_array);
	results_countpairs_wp *results = countpairs_wp(ND1,X1,Y1,Z1,
																								 boxsize,
#ifdef USE_OMP
																								 nthreads,
#endif
																								 binfile,
																								 pimax);
	
	/* Clean up. */
	Py_DECREF(x1_array);Py_DECREF(y1_array);Py_DECREF(z1_array);
/* 	for(int i=1;i<results->nbin;i++) { */
/* 		const DOUBLE rpavg = results->rpavg[i]; */
/* 		fprintf(stderr,"%lf %lf %lf %lf %"PRIu64"\n",results->rupp[i-1],results->rupp[i],rpavg,results->wp[i],results->npairs[i]); */
/* 	} */

	
	/* Build the output list */
	PyObject *ret = PyList_New(0);
	DOUBLE rlow=results->rupp[0];
	for(int i=1;i<results->nbin;i++) {
		PyObject *item = NULL;
		const DOUBLE rpavg = results->rpavg[i];

#ifdef DOUBLE_PREC		
	  item = Py_BuildValue("(ddddk)", rlow,results->rupp[i],rpavg,results->wp[i],results->npairs[i]);
#else
		item = Py_BuildValue("(ffffk)", rlow,results->rupp[i],rpavg,results->wp[i],results->npairs[i]);
#endif//DOUBLE_PREC

		PyList_Append(ret, item);
		Py_XDECREF(item);
		rlow=results->rupp[i];
	}
	free_results_wp(&results);
	return ret;
}


static PyObject *countpairs_countpairs_xi(PyObject *self, PyObject *args)
{
	(void) self;//to suppress the unused variable warning. Terrible hack
	PyObject *x1_obj, *y1_obj, *z1_obj;
	double boxsize;
	int nthreads=4;
	char *binfile;
	
	if (!PyArg_ParseTuple(args, "disOOO",&boxsize,&nthreads,&binfile,&x1_obj,&y1_obj,&z1_obj))
		return NULL;

	/* Interpret the input objects as numpy arrays. */
#ifdef DOUBLE_PREC	
	PyObject *x1_array = PyArray_FROM_OTF(x1_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
	PyObject *y1_array = PyArray_FROM_OTF(y1_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
	PyObject *z1_array = PyArray_FROM_OTF(z1_obj, NPY_DOUBLE,	NPY_ARRAY_IN_ARRAY);
#else
	PyObject *x1_array = PyArray_FROM_OTF(x1_obj, NPY_FLOAT, NPY_ARRAY_IN_ARRAY);
	PyObject *y1_array = PyArray_FROM_OTF(y1_obj, NPY_FLOAT, NPY_ARRAY_IN_ARRAY);
	PyObject *z1_array = PyArray_FROM_OTF(z1_obj, NPY_FLOAT, NPY_ARRAY_IN_ARRAY);
#endif
	
	if (x1_array == NULL || y1_array == NULL || z1_array == NULL) {
		Py_XDECREF(x1_array);
		Py_XDECREF(y1_array);
		Py_XDECREF(z1_array);
		return NULL;
	}


  /* How many data points are there? */
	const int64_t ND1 = (int64_t)PyArray_DIM(x1_array, 0);

	/* Get pointers to the data as C-types. */
	DOUBLE *X1 = (DOUBLE *)PyArray_DATA(x1_array);
	DOUBLE *Y1 = (DOUBLE *)PyArray_DATA(y1_array);
	DOUBLE *Z1 = (DOUBLE *)PyArray_DATA(z1_array);
	results_countpairs_xi *results = countpairs_xi(ND1,X1,Y1,Z1,
																								 boxsize,
#ifdef USE_OMP
																								 nthreads,
#endif
																								 binfile);
	
	/* Clean up. */
	Py_DECREF(x1_array);Py_DECREF(y1_array);Py_DECREF(z1_array);
/* 	for(int i=1;i<results->nbin;i++) { */
/* 		const DOUBLE rpavg = results->rpavg[i]; */
/* 		fprintf(stderr,"%lf %lf %lf %lf %"PRIu64"\n",results->rupp[i-1],results->rupp[i],rpavg,results->xi[i],results->npairs[i]); */
/* 	} */

	
	/* Build the output list */
	PyObject *ret = PyList_New(0);
	DOUBLE rlow=results->rupp[0];
	for(int i=1;i<results->nbin;i++) {
		PyObject *item = NULL;
		const DOUBLE rpavg = results->rpavg[i];

#ifdef DOUBLE_PREC		
	  item = Py_BuildValue("(ddddk)", rlow,results->rupp[i],rpavg,results->xi[i],results->npairs[i]);
#else
		item = Py_BuildValue("(ffffk)", rlow,results->rupp[i],rpavg,results->xi[i],results->npairs[i]);
#endif//DOUBLE_PREC

		PyList_Append(ret, item);
		Py_XDECREF(item);
		rlow=results->rupp[i];
	}
	free_results_xi(&results);
	return ret;
}

static PyObject *countpairs_countspheres_vpf(PyObject *self, PyObject *args)
{
	(void) self;//to suppress the unused variable warning. Terrible hack
	PyObject *x1_obj, *y1_obj, *z1_obj;
	double rmax;
	int nbin,nc,num_pN;
	unsigned long seed=-1;

	
	if (!PyArg_ParseTuple(args, "diiikOOO",&rmax,&nbin,&nc,&num_pN,&seed,&x1_obj,&y1_obj,&z1_obj))
		return NULL;

	/* Interpret the input objects as numpy arrays. */
#ifdef DOUBLE_PREC	
	PyObject *x1_array = PyArray_FROM_OTF(x1_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
	PyObject *y1_array = PyArray_FROM_OTF(y1_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
	PyObject *z1_array = PyArray_FROM_OTF(z1_obj, NPY_DOUBLE,	NPY_ARRAY_IN_ARRAY);
#else
	PyObject *x1_array = PyArray_FROM_OTF(x1_obj, NPY_FLOAT, NPY_ARRAY_IN_ARRAY);
	PyObject *y1_array = PyArray_FROM_OTF(y1_obj, NPY_FLOAT, NPY_ARRAY_IN_ARRAY);
	PyObject *z1_array = PyArray_FROM_OTF(z1_obj, NPY_FLOAT, NPY_ARRAY_IN_ARRAY);
#endif
	
	if (x1_array == NULL || y1_array == NULL || z1_array == NULL) {
		Py_XDECREF(x1_array);
		Py_XDECREF(y1_array);
		Py_XDECREF(z1_array);
		return NULL;
	}


  /* How many data points are there? */
	const int64_t ND1 = (int64_t)PyArray_DIM(x1_array, 0);

	/* Get pointers to the data as C-types. */
	DOUBLE *X1 = (DOUBLE *)PyArray_DATA(x1_array);
	DOUBLE *Y1 = (DOUBLE *)PyArray_DATA(y1_array);
	DOUBLE *Z1 = (DOUBLE *)PyArray_DATA(z1_array);

	/* Do the VPF calculation */
	results_countspheres *results = countspheres(ND1, X1, Y1, Z1,
																							 rmax, nbin, nc,
																							 num_pN,
																							 seed);

  /* Clean up. */
	Py_DECREF(x1_array);Py_DECREF(y1_array);Py_DECREF(z1_array);

	
	/* Build the output list (of lists, since num_pN is determined at runtime) */
	PyObject *ret = PyList_New(0);
	const DOUBLE rstep = rmax/(DOUBLE)nbin ;
	for(int ibin=0;ibin<results->nbin;ibin++) {
		const double r=(ibin+1)*rstep;
		PyObject *item = PyList_New(0);
		PyObject *this_val = Py_BuildValue("d",r);
		PyList_Append(item, this_val);
		Py_XDECREF(this_val);
		for(int i=0;i<num_pN;i++) {
#ifdef DOUBLE_PREC					
			this_val = Py_BuildValue("d",(results->pN)[ibin][i]);
#else
			this_val = Py_BuildValue("d",(results->pN)[ibin][i]);
#endif
			PyList_Append(item, this_val);
			Py_XDECREF(this_val);
		}
		PyList_Append(ret, item);
		Py_XDECREF(item);
	}
	
	free_results_countspheres(&results);
	return ret;
}
