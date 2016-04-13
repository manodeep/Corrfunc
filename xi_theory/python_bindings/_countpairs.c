/* File: _countpairs.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>

/* Now, include the numpy header*/
#include <arrayobject.h>

//for correlation functions
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
//python3 follows
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#define INITERROR return NULL
PyObject *PyInit__countpairs(void);

#else
//python2 follows
#define GETSTATE(m) (&_state)
static struct module_state _state;
#define INITERROR return
PyMODINIT_FUNC init_countpairs(void);

#endif


#ifdef DOUBLE_PREC
#define ELEMENT_SIZE NPY_DOUBLE
#else
#define ELEMENT_SIZE NPY_FLOAT
#endif

/* #ifndef PIMAX_UNICODE */
#define PI_UNICODE    "\u03C0"
/* #define XI_UNICODE    "\u039E" */
#define XI_UNICODE    "xi"
#define PIMAX_UNICODE "pimax"
/* #define RP_UNICODE    "r\u209a" */
/* #endif */



//Docstrings for the methods
static char module_docstring[]             =    "This module provides an interface for calculating clustering statistics using python extensions written in C.\n"
    "\n"
    "countpairs       : Calculate the 3-D xi auto/cross-correlation function given two sets of X1/Y1/Z1 and X2/Y2/Z2 arrays.\n"
    "countpairs_rp_pi : Calculate the 2-D DD(rp,pi) auto/cross-correlation function given two sets of X1/Y1/Z1 and X2/Y2/Z2 arrays.\n"
    "countpairs_wp    : Calculate the projected auto-correlation function wp (assumes PERIODIC) given one set of X1/Y1/Z1 arrays.\n"
    "countpairs_xi    : Calculate the 3-d auto-correlation function xi (assumes PERIODIC) given one set of X1/Y1/Z1 arrays.\n"
    "countpairs_vpf   : Calculate the counts-in-spheres given one set of X1/Y1/Z1 arrays.\n"
    "\n\n"
    "See `Corrfunc/call_correlation_functions.py` for example calls to each function.\n";

static char error_out_docstring[]          =  "Error-handler for the module.";

/* function proto-type*/
static PyObject *countpairs_countpairs(PyObject *self, PyObject *args);
static PyObject *countpairs_countpairs_rp_pi(PyObject *self, PyObject *args);
static PyObject *countpairs_countpairs_wp(PyObject *self, PyObject *args);
static PyObject *countpairs_countpairs_xi(PyObject *self, PyObject *args);
static PyObject *countpairs_countspheres_vpf(PyObject *self, PyObject *args);
static PyObject *countpairs_error_out(PyObject *module);


/* Inline documentation for the methods so that help(function) has something reasonably useful*/
static PyMethodDef module_methods[] = {
    {"countpairs_error_out"  ,(PyCFunction) countpairs_error_out        ,METH_NOARGS, error_out_docstring},
    {"countpairs"            ,(PyCFunction) countpairs_countpairs       ,METH_VARARGS,
     "countpairs(autocorr, nthreads, binfile, X1, Y1, Z1, X2, Y2, Z2)\n"
     "\n"
     "Calculate the 3-D pair-counts, "XI_UNICODE"(r), auto/cross-correlation function given two sets of points\n"
     "represented by X1/Y1/Z1 and X2/Y2/Z2 arrays.\n"
     "Note, that this module only returns pair counts and not the actual correlation function\n"
     ""XI_UNICODE"(r). See the xi_mocks/wtheta/wtheta.c for computing "XI_UNICODE"(r) from DD(r).\n"
     "\n"
     "parameters (all parameters are required)\n"
     "----------------------------------------\n"
     "autocorr: boolean\n"
     "    Flag for auto/cross-correlation. If autocorr is not 0, the X2/Y2/Z2 arrays\n"
     "    are not used (but must still be passed as valid arrays).\n"
     "nthreads: integer\n"
     "    The number of OpenMP threads to use. Has no effect if OpenMP was not used\n"
     "    during library compilation. \n"
     "binfile: filename\n"
     "    Filename containing the radial bins for the correlation function. The file\n"
     "    is expected to contain white-space separated ``rmin  rmax`` with the bin\n"
     "    edges. Units must be the same as the XYZ positions (see the ``bins`` file\n"
     "    in the tests directory for a sample). For usual logarithmic bins, ``logbins``\n"
     "   in the root directory of this package will create a compatible ``binfile``.\n"
     "X1/Y1/Z1: float/double (default float)\n"
     "    The array of X/Y/Z positions for the first set of points.\n"
     "    Float/double must be consistent with the library otherwise there will be a\n"
     "    RunTimeError from numpy. Specifically, if ``DOUBLE_PREC`` is not enabled,\n"
     "    then the arrays should be compatible with ``numpy.float32``; otherwise,\n"
     "    use ``numpy.float``. \n"
     "X2/Y2/Z2: float/double (default float)\n"
     "    Same as X1/Y1/Z1\n"
     "\n"
     "returns\n"
     "-------\n"
     "a Python list containing [rmin, rmax, ravg, npairs] for each radial bin\n"
     "specified in the ``binfile``. If ``OUTPUT_RPAVG`` is not defined in\n"
     "``theory.options`` then ``ravg`` will be set to 0.0 for all bins. ``npairs``\n"
     "contains the number of pairs in that bin and can be used to compute the\n"
     "actual "XI_UNICODE"(r) by combining with RR counts.\n"
    },
    {"countpairs_rp_pi"      ,(PyCFunction) countpairs_countpairs_rp_pi ,METH_VARARGS,
     "countpairs_rp_pi(autocorr, nthreads, "PIMAX_UNICODE", binfile, X1, Y1, Z1, X2, Y2, Z2)\n"
     "\n"
     "Calculate the 2-D pair-counts, "XI_UNICODE"(rp, "PI_UNICODE"), auto/cross-correlation function given two\n"
     "sets of X1/Y1/Z1 and X2/Y2/Z2 arrays. Assumes redshift-space distortions have already\n"
     "been applied. Uses z-direction as the line of sight ("PI_UNICODE") distance. Default bins in "PI_UNICODE"\n"
     "are set to 1.0 in position units.\n"
     "\n"
     "Note, that this module only returns pair counts and not the actual correlation function\n"
     ""XI_UNICODE"(rp, "PI_UNICODE"). See the xi_theory/xi_rp_pi/wprp.c for computing wp(rp) from DD(rp, "PI_UNICODE").\n"
     "\n"
     "parameters (all parameters are required)\n"
     "----------------------------------------\n"
     "autocorr: boolean\n"
     "    Flag for auto/cross-correlation. If autocorr is not 0, the X2/Y2/Z2 arrays\n"
     "    are not used (but must still be passed as valid arrays).\n"
     "nthreads: integer\n"
     "    The number of OpenMP threads to use. Has no effect if OpenMP was not used\n"
     "    during library compilation. \n"
     ""PIMAX_UNICODE": double\n"
     "    The max. integration distance along the z-direction ("PI_UNICODE"). Pairs with z\n"
     "    separation larger than "PIMAX_UNICODE" will not be counted. Units must have same as\n"
     "    those of positions specified in the X1/Y1/Z1 arrays.\n"
     "binfile: filename\n"
     "    Filename containing the radial bins for the correlation function. The file\n"
     "    is expected to contain white-space separated ``rmin  rmax`` with the bin\n"
     "    edges.  Units must be the same as the XYZ positions (see the ``bins`` file\n"
     "    in the tests directory for a sample). For usual logarithmic bins, ``logbins``\n"
     "    in the root directory of this package will create a compatible ``binfile``.\n"
     "X1/Y1/Z1: float/double (default float)\n"
     "    The array of X/Y/Z positions for the first set of points.\n"
     "    Float/double must be consistent with the library otherwise there will\n"
     "    be a RunTimeError from numpy. Specifically, if ``DOUBLE_PREC`` is not\n"
     "    enabled, then the arrays should be compatible with ``numpy.float32``;\n"
     "    otherwise, use ``numpy.float``.\n"
     "X2/Y2/Z2: float/double (default float)\n"
     "    Same as X1/Y1/Z1\n"
     "\n"
     "returns\n"
     "-------\n"
     "a Python list containing [rmin, rmax, ravg, "PI_UNICODE", npairs] \n"
     "for each "PI_UNICODE"-bin (up to "PIMAX_UNICODE") for each radial bin specified in\n"
     "the ``binfile``. For instance, for a ``"PIMAX_UNICODE"`` of 40.0 Mpc/h, each radial\n"
     "bin will be split into 40 "PI_UNICODE" bins (default "PI_UNICODE" bin is 1.0). Thus, the\n"
     "total number of items in the list is {(int) ``"PIMAX_UNICODE"`` * number of rp bins}.\n"
     "If ``OUTPUT_RPAVG`` is not defined in ``theory.options`` \n"
     "then ``ravg`` will be set to 0.0 for all bins. "PI_UNICODE" for each bin\n"
     "is the upper limit of the "PI_UNICODE" values that were considered in that (rp, "PI_UNICODE") bin.\n"
     "``npairs``contains the number of pairs in that bin and can be used to compute the\n"
     "actual "XI_UNICODE"(rp, "PI_UNICODE") by combining with RR counts.\n"
     "\n"
     "example\n"
     "-------\n"
     "import numpy as np\n"
     "from Corrfunc._countpairs import countpairs\n"
     "x,y,z = np.genfromtxt('/path/to/ascii/galaxy/file/(x y z)',dtype=np.float32,unpack=True)\n"
     "autocorr=1\n"
     "nthreads=4\n"
     "DD = countpairs(autocorr,nthreads,'../xi_theory/tests/bins',x,y,z,x,y,z)\n"
     "\n"
    },
    {"countpairs_wp"         ,(PyCFunction) countpairs_countpairs_wp    ,METH_VARARGS,
     "countpairs_wp(boxsize, "PIMAX_UNICODE", nthreads, binfile, X1, Y1, Z1)\n"
     "\n"
     "Calculates the projected 2-pt auto-correlation function, wp(rp), on periodic boxes from X1/Y1/Z1.\n"
     "Assumes redshift-space distortions have already been applied. Uses z-direction as the line\n"
     "of sight ("PI_UNICODE") distance. *Always* uses ``PERIODIC`` boundary conditions.\n"
     "\n"
     "Note, that this module returns the actual correlation function using the natural estimator.\n"
     "Analytic randoms are used to compute wp(rp) from the pair counts. If you need a different estimator,\n"
     "Landy-Szalay for instance, then you should use the module countpairs_rp_pi and use the pair counts\n"
     "to obtain the Landy-Szalay estimator for wp(rp).\n"
     "\n"
     "parameters (all parameters are required)\n"
     "----------------------------------------\n"
     "boxsize: double\n"
     "    The size of the cosmological box from which the points are generated. Should have\n"
     "    same units as the positions.\n"
     ""PIMAX_UNICODE": double\n"
     "    The max. integration distance along the z-direction ("PI_UNICODE"). Pairs with z\n"
     "    separation larger than "PIMAX_UNICODE" will not be counted. Units must have same as\n"
     "    those of positions specified in the X1/Y1/Z1 arrays.\n"
     "nthreads: integer\n"
     "    The number of OpenMP threads to use. Has no effect if OpenMP was not used\n"
     "    during library compilation. \n"
     "binfile: filename\n"
     "    Filename containing the radial bins for the correlation function. The file\n"
     "    is expected to contain white-space separated ``rmin  rmax`` with the bin\n"
     "    edges.  Units must be the same as the XYZ positions (see the ``bins`` file\n"
     "    in the tests directory for a sample). For usual logarithmic bins, ``logbins``\n"
     "    in the root directory of this package will create a compatible ``binfile``.\n"
     "X1/Y1/Z1: float/double (default float)\n"
     "    The array of X/Y/Z positions for the first set of points.\n"
     "    Float/double must be consistent with the library otherwise there will\n"
     "    be a RunTimeError from numpy. Specifically, if ``DOUBLE_PREC`` is not\n"
     "    enabled, then the arrays should be compatible with ``numpy.float32``;\n"
     "    otherwise, use ``numpy.float``.\n"
     "\n"
     "returns\n"
     "-------\n"
     "a Python list containing [rmin, rmax, ravg, wp, npairs] for each radial bin\n"
     "specified in the ``binfile``. If ``OUTPUT_RPAVG`` is not defined in\n"
     "``theory.options`` then ``ravg`` will be set to 0.0 for all bins. ``wp``\n"
     "contains the projected correlation function while ``npairs`` contains the\n"
     "number of pairs in that bin.\n"
    },
    {"countpairs_xi"         ,(PyCFunction) countpairs_countpairs_xi    ,METH_VARARGS,
     "countpairs_xi(boxsize, nthreads, binfile, X1, Y1, Z1)\n"
     "\n"
     "Calculates the 3-D 2-pt auto-correlation function, xi(r), on periodic boxes from\n"
     "X1/Y1/Z1. Assumes redshift-space distortions have already been applied. *Always* uses\n"
     "``PERIODIC`` boundary conditions.\n"
     "\n"
     "Note, that this module returns the actual correlation function using the natural estimator.\n"
     "Analytic randoms are used to compute xi(r) from the pair counts. If you need a different estimator,\n"
     "Landy-Szalay for instance, then you should use the module countpairs and use the pair counts\n"
     "to obtain the Landy-Szalay estimator for xi(r).\n"
     "\n"
     "parameters (all parameters are required)\n"
     "----------------------------------------\n"
     "boxsize: double\n"
     "    The size of the cosmological box from which the points are generated. Should have\n"
     "    same units as the positions.\n"
     "nthreads: integer\n"
     "    The number of OpenMP threads to use. Has no effect if OpenMP was not used\n"
     "    during library compilation. \n"
     "binfile: filename\n"
     "    Filename containing the radial bins for the correlation function. The file\n"
     "    is expected to contain white-space separated ``rmin  rmax`` with the bin\n"
     "    edges.  Units must be the same as the XYZ positions (see the ``bins`` file\n"
     "    in the tests directory for a sample). For usual logarithmic bins, ``logbins``\n"
     "    in the root directory of this package will create a compatible ``binfile``.\n"
     "X1/Y1/Z1: float/double (default float)\n"
     "    The array of X/Y/Z positions for the first set of points.\n"
     "    Float/double must be consistent with the library otherwise there will\n"
     "    be a RunTimeError from numpy. Specifically, if ``DOUBLE_PREC`` is not\n"
     "    enabled, then the arrays should be compatible with ``numpy.float32``;\n"
     "    otherwise, use ``numpy.float``.\n"
     "\n"
     "returns\n"
     "-------\n"
     "a Python list containing [rmin, rmax, ravg, xi, npairs] for each radial bin\n"
     "specified in the ``binfile``. If ``OUTPUT_RPAVG`` is not defined in\n"
     "``theory.options`` then ``ravg`` will be set to 0.0 for all bins. ``xi``\n"
     "contains the 3-D autocorrelation function while ``npairs`` contains the\n"
     "number of pairs in that bin.\n"
    },
    {"countspheres_vpf"      ,(PyCFunction) countpairs_countspheres_vpf ,METH_VARARGS,
     "countspheres_vpf(rmax, nbin, ncenters, num_pN, seed, X1, Y1, Z1)\n"
     "\n"
     "Calculates the fraction of random spheres that contain exactly *N* points, pN(r).\n"
     "\n"
     "parameters (all parameters are required)\n"
     "----------------------------------------\n"
     "rmax: double\n"
     "    Largest sphere radius to calculate. Should have same units as the positions.\n"
     "nbin: integer\n"
     "    Number of linear bins to divide rmax into. Typical value is numerically equal\n"
     "    to rmax, so that each radial bin increments by 1.0.\n"
     "ncenters: integer\n"
     "    Number of random spheres to place on the point distribution. Too small a value\n"
     "    will give noisy results, while too large a value only increases runtime without\n"
     "    producing any additional information.\n"
     "num_pN: integer (>=1)\n"
     "    The number of counts in a cell to consider. The void probability function, vpf (or p0),\n"
     "    is returned for num_pN=1. Note that num_pN does not refer to the actual counts in the\n"
     "    cell; num_pN refers to the total number of such counts to return. For num_pN=2, the\n"
     "    code will return the fraction of cells that contain exactly 0 and exactly 1 points. \n"
     "    More explicitly, the columns in the results look like the following:\n"
     "    num_pN = 1 -> p0 \n"
     "    num_pN = 2 -> p0 p1\n"
     "    num_pN = 3 -> p0 p1 p2\n"
     "    and so on...(note that p0 is the vpf).\n"
     "seed: unsigned long\n"
     "    Random number seed for the RNG. Used version is the GSL Mersenne-Twister (mt19937).\n"
     "X1/Y1/Z1: float/double (default float)\n"
     "    The array of X/Y/Z positions for the first set of points.\n"
     "    Float/double must be consistent with the library otherwise there will\n"
     "    be a RunTimeError from numpy. Specifically, if ``DOUBLE_PREC`` is not\n"
     "    enabled, then the arrays should be compatible with ``numpy.float32``;\n"
     "    otherwise, use ``numpy.float``.\n"
     "\n"
     "returns\n"
     "-------\n"
     "a Python list containing [rmax, p0, p1,..., p(num_pN-1)] for each radial bin.\n"
    },
    {NULL, NULL, 0, NULL}
};

static PyObject *countpairs_error_out(PyObject *module) {
    (void) module;//to avoid unused warning
    struct module_state *st = GETSTATE(module);
    PyErr_SetString(st->error, "something bad happened");
    return NULL;
}


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
    module_docstring,
    sizeof(struct module_state),
    module_methods,
    NULL,
    _countpairs_traverse,
    _countpairs_clear,
    NULL
};


PyObject *PyInit__countpairs(void)

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

    st->error = PyErr_NewException("_countpairs.error", NULL, NULL);
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

    /* Interpret the input objects as numpy arrays. ELEMENT_SIZE macro expands to NPY_DOUBLE or NPY_FLOAT*/
    PyObject *x1_array = PyArray_FROM_OTF(x1_obj, ELEMENT_SIZE, NPY_ARRAY_IN_ARRAY);
    PyObject *y1_array = PyArray_FROM_OTF(y1_obj, ELEMENT_SIZE, NPY_ARRAY_IN_ARRAY);
    PyObject *z1_array = PyArray_FROM_OTF(z1_obj, ELEMENT_SIZE, NPY_ARRAY_IN_ARRAY);

    if (x1_array == NULL || y1_array == NULL || z1_array == NULL) {
        Py_XDECREF(x1_array);
        Py_XDECREF(y1_array);
        Py_XDECREF(z1_array);
        Py_RETURN_NONE;
    }


    /* Interpret the input objects as numpy arrays. */
    PyObject *x2_array = PyArray_FROM_OTF(x2_obj, ELEMENT_SIZE, NPY_ARRAY_IN_ARRAY);
    PyObject *y2_array = PyArray_FROM_OTF(y2_obj, ELEMENT_SIZE, NPY_ARRAY_IN_ARRAY);
    PyObject *z2_array = PyArray_FROM_OTF(z2_obj, ELEMENT_SIZE, NPY_ARRAY_IN_ARRAY);

    if (x2_array == NULL || y2_array == NULL || z2_array == NULL) {
        Py_XDECREF(x2_array);
        Py_XDECREF(y2_array);
        Py_XDECREF(z2_array);
        Py_RETURN_NONE;
    }


    /* How many data points are there? */
    const int64_t ND1 = (int64_t)PyArray_DIM((const PyArrayObject *) x1_array, 0);
    const int64_t ND2 = (int64_t)PyArray_DIM((const PyArrayObject *) x2_array, 0);

    /* Get pointers to the data as C-types. */
    DOUBLE *X1 = (DOUBLE *)PyArray_DATA((PyArrayObject *) x1_array);
    DOUBLE *Y1 = (DOUBLE *)PyArray_DATA((PyArrayObject *) y1_array);
    DOUBLE *Z1 = (DOUBLE *)PyArray_DATA((PyArrayObject *) z1_array);

    DOUBLE *X2 = (DOUBLE *)PyArray_DATA((PyArrayObject *) x2_array);
    DOUBLE *Y2 = (DOUBLE *)PyArray_DATA((PyArrayObject *) y2_array);
    DOUBLE *Z2 = (DOUBLE *)PyArray_DATA((PyArrayObject *) z2_array);


    results_countpairs *results = countpairs(ND1,X1,Y1,Z1,
                                             ND2,X2,Y2,Z2,
#if defined(USE_OMP) && defined(_OPENMP)
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
    PyObject *x1_array = PyArray_FROM_OTF(x1_obj, ELEMENT_SIZE, NPY_ARRAY_IN_ARRAY);
    PyObject *y1_array = PyArray_FROM_OTF(y1_obj, ELEMENT_SIZE, NPY_ARRAY_IN_ARRAY);
    PyObject *z1_array = PyArray_FROM_OTF(z1_obj, ELEMENT_SIZE, NPY_ARRAY_IN_ARRAY);

    if (x1_array == NULL || y1_array == NULL || z1_array == NULL) {
        Py_XDECREF(x1_array);
        Py_XDECREF(y1_array);
        Py_XDECREF(z1_array);
        return NULL;
    }


    /* Interpret the input objects as numpy arrays. */
    PyObject *x2_array = PyArray_FROM_OTF(x2_obj, ELEMENT_SIZE, NPY_ARRAY_IN_ARRAY);
    PyObject *y2_array = PyArray_FROM_OTF(y2_obj, ELEMENT_SIZE, NPY_ARRAY_IN_ARRAY);
    PyObject *z2_array = PyArray_FROM_OTF(z2_obj, ELEMENT_SIZE, NPY_ARRAY_IN_ARRAY);

    if (x2_array == NULL || y2_array == NULL || z2_array == NULL) {
        Py_XDECREF(x2_array);
        Py_XDECREF(y2_array);
        Py_XDECREF(z2_array);
        return NULL;
    }


    /* How many data points are there? */
    const int64_t ND1 = (int64_t)PyArray_DIM((const PyArrayObject *) x1_array, 0);
    const int64_t ND2 = (int64_t)PyArray_DIM((const PyArrayObject *) x2_array, 0);

    /* Get pointers to the data as C-types. */
    DOUBLE *X1 = (DOUBLE *)PyArray_DATA((PyArrayObject *) x1_array);
    DOUBLE *Y1 = (DOUBLE *)PyArray_DATA((PyArrayObject *) y1_array);
    DOUBLE *Z1 = (DOUBLE *)PyArray_DATA((PyArrayObject *) z1_array);

    DOUBLE *X2 = (DOUBLE *)PyArray_DATA((PyArrayObject *) x2_array);
    DOUBLE *Y2 = (DOUBLE *)PyArray_DATA((PyArrayObject *) y2_array);
    DOUBLE *Z2 = (DOUBLE *)PyArray_DATA((PyArrayObject *) z2_array);

    results_countpairs_rp_pi *results = countpairs_rp_pi(ND1,X1,Y1,Z1,
                                                         ND2,X2,Y2,Z2,
#if defined(USE_OMP) && defined(_OPENMP)
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
    PyObject *x1_array = PyArray_FROM_OTF(x1_obj, ELEMENT_SIZE, NPY_ARRAY_IN_ARRAY);
    PyObject *y1_array = PyArray_FROM_OTF(y1_obj, ELEMENT_SIZE, NPY_ARRAY_IN_ARRAY);
    PyObject *z1_array = PyArray_FROM_OTF(z1_obj, ELEMENT_SIZE, NPY_ARRAY_IN_ARRAY);

    if (x1_array == NULL || y1_array == NULL || z1_array == NULL) {
        Py_XDECREF(x1_array);
        Py_XDECREF(y1_array);
        Py_XDECREF(z1_array);
        return NULL;
    }


    /* How many data points are there? */
    const int64_t ND1 = (int64_t)PyArray_DIM((const PyArrayObject *) x1_array, 0);

    /* Get pointers to the data as C-types. */
    DOUBLE *X1 = (DOUBLE *)PyArray_DATA((PyArrayObject *) x1_array);
    DOUBLE *Y1 = (DOUBLE *)PyArray_DATA((PyArrayObject *) y1_array);
    DOUBLE *Z1 = (DOUBLE *)PyArray_DATA((PyArrayObject *) z1_array);

    results_countpairs_wp *results = countpairs_wp(ND1,X1,Y1,Z1,
                                                   boxsize,
#if defined(USE_OMP) && defined(_OPENMP)
                                                   nthreads,
#endif
                                                   binfile,
                                                   pimax);

    /* Clean up. */
    Py_DECREF(x1_array);Py_DECREF(y1_array);Py_DECREF(z1_array);

#if 0
    for(int i=1;i<results->nbin;i++) {
        const DOUBLE rpavg = results->rpavg[i];
        fprintf(stderr,"%lf %lf %lf %lf %"PRIu64"\n",results->rupp[i-1],results->rupp[i],rpavg,results->wp[i],results->npairs[i]);
    }
#endif

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
    PyObject *x1_array = PyArray_FROM_OTF(x1_obj, ELEMENT_SIZE, NPY_ARRAY_IN_ARRAY);
    PyObject *y1_array = PyArray_FROM_OTF(y1_obj, ELEMENT_SIZE, NPY_ARRAY_IN_ARRAY);
    PyObject *z1_array = PyArray_FROM_OTF(z1_obj, ELEMENT_SIZE, NPY_ARRAY_IN_ARRAY);

    if (x1_array == NULL || y1_array == NULL || z1_array == NULL) {
        Py_XDECREF(x1_array);
        Py_XDECREF(y1_array);
        Py_XDECREF(z1_array);
        return NULL;
    }


    /* How many data points are there? */
    const int64_t ND1 = (int64_t)PyArray_DIM((const PyArrayObject *) x1_array, 0);

    /* Get pointers to the data as C-types. */
    DOUBLE *X1 = (DOUBLE *)PyArray_DATA((PyArrayObject *) x1_array);
    DOUBLE *Y1 = (DOUBLE *)PyArray_DATA((PyArrayObject *) y1_array);
    DOUBLE *Z1 = (DOUBLE *)PyArray_DATA((PyArrayObject *) z1_array);

    results_countpairs_xi *results = countpairs_xi(ND1,X1,Y1,Z1,
                                                   boxsize,
#if defined(USE_OMP) && defined(_OPENMP)
                                                   nthreads,
#endif
                                                   binfile);

    /* Clean up. */
    Py_DECREF(x1_array);Py_DECREF(y1_array);Py_DECREF(z1_array);

#if 0
    for(int i=1;i<results->nbin;i++) {
        const DOUBLE rpavg = results->rpavg[i];
        fprintf(stderr,"%lf %lf %lf %lf %"PRIu64"\n",results->rupp[i-1],results->rupp[i],rpavg,results->xi[i],results->npairs[i]);
    }
#endif

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
    PyObject *x1_array = PyArray_FROM_OTF(x1_obj, ELEMENT_SIZE, NPY_ARRAY_IN_ARRAY);
    PyObject *y1_array = PyArray_FROM_OTF(y1_obj, ELEMENT_SIZE, NPY_ARRAY_IN_ARRAY);
    PyObject *z1_array = PyArray_FROM_OTF(z1_obj, ELEMENT_SIZE, NPY_ARRAY_IN_ARRAY);

    if (x1_array == NULL || y1_array == NULL || z1_array == NULL) {
        Py_XDECREF(x1_array);
        Py_XDECREF(y1_array);
        Py_XDECREF(z1_array);
        return NULL;
    }

    /* How many data points are there? */
    const int64_t ND1 = (int64_t)PyArray_DIM((const PyArrayObject *) x1_array, 0);

    /* Get pointers to the data as C-types. */
    DOUBLE *X1 = (DOUBLE *)PyArray_DATA((PyArrayObject *) x1_array);
    DOUBLE *Y1 = (DOUBLE *)PyArray_DATA((PyArrayObject *) y1_array);
    DOUBLE *Z1 = (DOUBLE *)PyArray_DATA((PyArrayObject *) z1_array);

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
