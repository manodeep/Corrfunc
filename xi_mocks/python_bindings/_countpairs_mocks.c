/* File: _countpairs_mocks.c */
/*
		This file is a part of the Corrfunc package
		Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
		License: MIT LICENSE. See LICENSE file under the top-level
		directory at https://bitbucket.org/manodeep/corrfunc/
*/

#include <Python.h>
#include <arrayobject.h>

#include "countpairs_rp_pi_mocks.h"
#include "countpairs_theta_mocks.h"

//for the vpf
#include "countspheres_mocks.h"

//Docstrings for the methods
static char module_docstring[]             =	"This module provides an interface for calculating correlation functions on MOCKS (spherical geometry) using C.";
static char countpairs_rp_pi_docstring[]   =	"Calculate the 2-D DD(rp,pi) auto/cross-correlation function given two sets of ra/dec/cz and ra/dec/cz arrays.";
static char countpairs_theta_docstring[]   =	"Calculate DD(theta) auto/cross-correlation function given two sets of ra/dec/cz and ra/dec/cz arrays.";
static char countspheres_vpf_docstring[]   =	"Calculate the counts-in-spheres given one set of X1/Y1/Z1 arrays.";

/* function proto-type*/
static PyObject *countpairs_countpairs_rp_pi_mocks(PyObject *self, PyObject *args);
static PyObject *countpairs_countpairs_theta_mocks(PyObject *self, PyObject *args);
static PyObject *countpairs_countspheres_vpf_mocks(PyObject *self, PyObject *args);

PyMODINIT_FUNC init_countpairs_mocks(void);

static PyMethodDef module_methods[] = {
	{"countpairs_rp_pi_mocks" ,countpairs_countpairs_rp_pi_mocks ,METH_VARARGS,countpairs_rp_pi_docstring},
	{"countpairs_theta_mocks" ,countpairs_countpairs_theta_mocks ,METH_VARARGS,countpairs_theta_docstring},
	{"countspheres_vpf_mocks" ,countpairs_countspheres_vpf_mocks ,METH_VARARGS,countspheres_vpf_docstring},
	{NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC init_countpairs_mocks(void)
{
	PyObject *m = Py_InitModule3("_countpairs_mocks", module_methods, module_docstring);
	if (m == NULL)
		return;

	/* Load `numpy` functionality. */
	import_array();
}



static PyObject *countpairs_countpairs_rp_pi_mocks(PyObject *self, PyObject *args)
{
	(void) self;//to suppress the unused variable warning. Terrible hack

	//x1->ra (phi), y1-> declination (theta1), z1->cz (cz1)
	//x2->ra (ph2), y2-> declination (theta2), z2->cz (cz2)
	PyObject *x1_obj, *y1_obj, *z1_obj, *x2_obj,*y2_obj,*z2_obj;
	int autocorr=0;
	int nthreads=4;
	int cosmology=1;
	double pimax;
	char *binfile;
	
	if (!PyArg_ParseTuple(args, "iiidsOOOOOO",&autocorr,&cosmology,&nthreads,&pimax,&binfile,&x1_obj,&y1_obj,&z1_obj,&x2_obj,&y2_obj,&z2_obj))
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
	const int64_t ND1 = (int64_t)PyArray_DIM(x1_array, 0);
	const int64_t ND2 = (int64_t)PyArray_DIM(x2_array, 0);

	/* Get pointers to the data as C-types. */
	DOUBLE *phiD1   = (DOUBLE *)PyArray_DATA(x1_array);
	DOUBLE *thetaD1 = (DOUBLE *)PyArray_DATA(y1_array);
	DOUBLE *czD1    = (DOUBLE *)PyArray_DATA(z1_array);

	DOUBLE *phiD2   = (DOUBLE *)PyArray_DATA(x2_array);
	DOUBLE *thetaD2 = (DOUBLE *)PyArray_DATA(y2_array);
	DOUBLE *czD2    = (DOUBLE *)PyArray_DATA(z2_array);
		

  results_countpairs_mocks *results  = countpairs_mocks(ND1,phiD1,thetaD1,czD1,
																												ND2,phiD2,thetaD2,czD2,
#ifdef USE_OMP
																												nthreads,
#endif
																												autocorr,
																												binfile,
																												pimax,
																												cosmology);
	
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
	free_results_mocks(&results);
	return ret;
}

static PyObject *countpairs_countpairs_theta_mocks(PyObject *self, PyObject *args)
{
	(void) self;//to suppress the unused variable warning. Terrible hack
	PyObject *x1_obj, *y1_obj, *x2_obj,*y2_obj;
	int nthreads=4;
	char *binfile;
	int autocorr=0;
	int cosmology=1;
	
	if (!PyArg_ParseTuple(args, "iiisOOOO",&autocorr,&cosmology,&nthreads,&binfile,&x1_obj,&y1_obj,&x2_obj,&y2_obj))
		return NULL;

	/* Interpret the input objects as numpy arrays. */
#ifdef DOUBLE_PREC	
	PyObject *x1_array = PyArray_FROM_OTF(x1_obj, NPY_DOUBLE, NPY_IN_ARRAY);
	PyObject *y1_array = PyArray_FROM_OTF(y1_obj, NPY_DOUBLE, NPY_IN_ARRAY);
#else
	PyObject *x1_array = PyArray_FROM_OTF(x1_obj, NPY_FLOAT, NPY_IN_ARRAY);
	PyObject *y1_array = PyArray_FROM_OTF(y1_obj, NPY_FLOAT, NPY_IN_ARRAY);
#endif
	
	if (x1_array == NULL || y1_array == NULL) {
		Py_XDECREF(x1_array);
		Py_XDECREF(y1_array);
		return NULL;
	}

	/* Interpret the input objects as numpy arrays. */
#ifdef DOUBLE_PREC	
	PyObject *x2_array = PyArray_FROM_OTF(x2_obj, NPY_DOUBLE, NPY_IN_ARRAY);
	PyObject *y2_array = PyArray_FROM_OTF(y2_obj, NPY_DOUBLE, NPY_IN_ARRAY);
#else
	PyObject *x2_array = PyArray_FROM_OTF(x2_obj, NPY_FLOAT, NPY_IN_ARRAY);
	PyObject *y2_array = PyArray_FROM_OTF(y2_obj, NPY_FLOAT, NPY_IN_ARRAY);
#endif
	
	if (x2_array == NULL || y2_array == NULL) {
		Py_XDECREF(x2_array);
		Py_XDECREF(y2_array);
		return NULL;
	}

	
  /* How many data points are there? */
	const int64_t ND1 = (int64_t)PyArray_DIM(x1_array, 0);
	const int64_t ND2 = (int64_t)PyArray_DIM(x2_array, 0);

	/* Get pointers to the data as C-types. */
	DOUBLE *phiD1   = (DOUBLE *)PyArray_DATA(x1_array);
	DOUBLE *thetaD1 = (DOUBLE *)PyArray_DATA(y1_array);

	DOUBLE *phiD2   = (DOUBLE *)PyArray_DATA(x2_array);
	DOUBLE *thetaD2 = (DOUBLE *)PyArray_DATA(y2_array);
	

  results_countpairs_theta *results = countpairs_theta_mocks(ND1,phiD1,thetaD1,
																														 ND2,phiD2,thetaD2,
#ifdef USE_OMP
																														 nthreads,
#endif
																														 autocorr,
																														 binfile) ;
	
	
	/* Clean up. */
	Py_DECREF(x1_array);Py_DECREF(y1_array);
	Py_DECREF(x2_array);Py_DECREF(y2_array);
	
  /*---Output-Pairs-------------------------------------*/
	/* DOUBLE theta_low = results->theta_upp[0]; */
	/* for(int i=1;i<results->nbin;i++) { */
	/* 	fprintf(stdout,"%10"PRIu64" %20.8lf %20.8lf %20.8lf \n",results->npairs[i],results->theta_avg[i],theta_low,results->theta_upp[i]); */
	/* 	theta_low=results->theta_upp[i]; */
	/* } */
	
	
	/* Build the output list */
	PyObject *ret = PyList_New(0);
	DOUBLE rlow=results->theta_upp[0];
	for(int i=1;i<results->nbin;i++) {
		PyObject *item = NULL;
		const DOUBLE theta_avg = results->theta_avg[i];
#ifdef DOUBLE_PREC		
	  item = Py_BuildValue("(dddk)", rlow,results->theta_upp[i],theta_avg,results->npairs[i]);
#else
		item = Py_BuildValue("(fffk)", rlow,results->theta_upp[i],theta_avg,results->npairs[i]);
#endif//DOUBLE_PREC

		PyList_Append(ret, item);
		Py_XDECREF(item);
		rlow=results->theta_upp[i];
	}
	free_results_countpairs_theta(&results);
	return ret;
}



static PyObject *countpairs_countspheres_vpf_mocks(PyObject *self, PyObject *args)
{
	(void) self;//to suppress the unused variable warning. Terrible hack

	//x1->ra (phi), y1-> declination (theta1), z1->cz (cz1)
	//x2->ra (ph2), y2-> declination (theta2), z2->cz (cz2)
	PyObject *x1_obj, *y1_obj, *z1_obj, *x2_obj,*y2_obj,*z2_obj;
	int cosmology=1;
	double rmax;
	int nbin,num_spheres, num_pN;
	char *centers_file;
	int threshold_neighbors;
	
	if (!PyArg_ParseTuple(args, "diiiisiOOOOOO",&rmax,&nbin,&num_spheres,&num_pN,&threshold_neighbors,&centers_file,&cosmology,&x1_obj,&y1_obj,&z1_obj,&x2_obj,&y2_obj,&z2_obj))
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
	const int64_t ND1 = (int64_t)PyArray_DIM(x1_array, 0);
	const int64_t ND2 = (int64_t)PyArray_DIM(x2_array, 0);

	/* Get pointers to the data as C-types. */
	DOUBLE *phiD1   = (DOUBLE *)PyArray_DATA(x1_array);
	DOUBLE *thetaD1 = (DOUBLE *)PyArray_DATA(y1_array);
	DOUBLE *czD1    = (DOUBLE *)PyArray_DATA(z1_array);

	DOUBLE *phiD2   = (DOUBLE *)PyArray_DATA(x2_array);
	DOUBLE *thetaD2 = (DOUBLE *)PyArray_DATA(y2_array);
	DOUBLE *czD2    = (DOUBLE *)PyArray_DATA(z2_array);


	results_countspheres_mocks *results = countspheres_mocks(ND1, phiD1,thetaD1, czD1,
																													 ND2, phiD2,thetaD2, czD2,
																													 threshold_neighbors,
																													 rmax, nbin, num_spheres,
																													 num_pN,
																													 centers_file,
																													 cosmology);



  //Output the results
	/* const DOUBLE rstep = rmax/(DOUBLE)nbin ; */
	/* for(int ibin=0;ibin<results->nbin;ibin++) { */
	/*      const double r=(ibin+1)*rstep; */
	/*      fprintf(stdout,"%10.2"DOUBLE_FORMAT" ", r); */
	/*      for(int i=0;i<num_pN;i++) { */
	/*              fprintf(stdout," %10.4e", (results->pN)[ibin][i]); */
	/*      } */
	/*      fprintf(stdout,"\n"); */
	/* } */

	
	/* Clean up. */
	Py_DECREF(x1_array);Py_DECREF(y1_array);Py_DECREF(z1_array);
	Py_DECREF(x2_array);Py_DECREF(y2_array);Py_DECREF(z2_array);


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
	free_results_countspheres_mocks(&results);
	

	return ret;
}
