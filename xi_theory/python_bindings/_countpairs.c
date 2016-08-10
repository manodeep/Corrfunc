/* File: _countpairs.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/
#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <stdio.h>

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
#define ELEMENT_TYPE NPY_DOUBLE
#else
#define ELEMENT_TYPE NPY_FLOAT
#endif

#define ELEMENT_DESCR    (PyArray_DescrFromType(ELEMENT_TYPE))
#define NOTYPE_DESCR     (PyArray_DescrFromType(NPY_NOTYPE))

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
static PyObject *countpairs_error_out(PyObject *module, const char *msg);


/* Inline documentation for the methods so that help(function) has something reasonably useful*/
static PyMethodDef module_methods[] = {
    {"countpairs_error_out"  ,(PyCFunction) countpairs_error_out        ,METH_VARARGS, error_out_docstring},
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

static PyObject *countpairs_error_out(PyObject *module, const char *msg)
{
#if PY_MAJOR_VERSION < 3
    (void) module;//to avoid unused warning with python2
#endif    

    struct module_state *st = GETSTATE(module);
    PyErr_SetString(st->error, msg);
    Py_RETURN_NONE;
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
//Python 2
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

    
static int64_t check_dims_and_datatype(PyObject *module, const PyArrayObject *x1_obj, const PyArrayObject *y1_obj, const PyArrayObject *z1_obj, size_t *element_size)
{
    char msg[1024];

    /* All the arrays should be 1-D*/
    const int nxdims = PyArray_NDIM(x1_obj);
    const int nydims = PyArray_NDIM(y1_obj);
    const int nzdims = PyArray_NDIM(z1_obj);

    if(nxdims != 1 || nydims != 1 || nzdims != 1) {
        snprintf(msg, 1024, "ERROR: Expected 1-D numpy arrays.\nFound (nxdims, nydims, nzdims) = (%d, %d, %d) instead",
                 nxdims, nydims, nzdims);
        countpairs_error_out(module, msg);
        return -1;
    }

    /* All the arrays should be floating point (only float32 and float64 are allowed) */
    const int x_type = PyArray_TYPE(x1_obj);
    const int y_type = PyArray_TYPE(y1_obj);
    const int z_type = PyArray_TYPE(z1_obj);
    if( ! ((x_type == NPY_FLOAT || x_type == NPY_DOUBLE) &&
           (y_type == NPY_FLOAT || y_type == NPY_DOUBLE) &&
           (z_type == NPY_FLOAT || z_type == NPY_DOUBLE))
        ) {
        PyArray_Descr *x_descr = PyArray_DescrFromType(x_type);
        PyArray_Descr *y_descr = PyArray_DescrFromType(y_type);
        PyArray_Descr *z_descr = PyArray_DescrFromType(z_type);
        if(x_descr == NULL || y_descr == NULL || z_descr == NULL) {
            /* Generating the dtype descriptor failed somehow. At least provide some information */
            snprintf(msg, 1024, "TypeError: Expected 3 floating point arrays (allowed types = %d or %d). Instead found type-nums (%d, %d, %d)\n",
                     NPY_FLOAT, NPY_DOUBLE, x_type, y_type, z_type);
        } else {
            snprintf(msg, 1024, "TypeError: Expected 3 floating point arrays (allowed types = %d or %d). Instead found type-nums (%d, %d, %d) "
                     "with type-names = (%s, %s, %s)\n",
                     NPY_FLOAT, NPY_DOUBLE, x_type, y_type, z_type, x_descr->typeobj->tp_name, y_descr->typeobj->tp_name, z_descr->typeobj->tp_name);
        }
        Py_XDECREF(x_descr);Py_XDECREF(y_descr);Py_XDECREF(z_descr);
        countpairs_error_out(module, msg);
        return -1;
    }

    if( x_type != y_type || y_type != z_type) {
        PyArray_Descr *x_descr = PyArray_DescrFromType(x_type);
        PyArray_Descr *y_descr = PyArray_DescrFromType(y_type);
        PyArray_Descr *z_descr = PyArray_DescrFromType(z_type);
        if(x_descr == NULL || y_descr == NULL || z_descr == NULL) {
          /* Generating the dtype descriptor failed somehow. At least provide some information */
          snprintf(msg, 1024, "TypeError: Expected *ALL* 3 floating point arrays to be the same type (allowed types = %d or %d). Instead found type-nums (%d, %d, %d)\n",
                   NPY_FLOAT, NPY_DOUBLE, x_type, y_type, z_type);
        } else {
          snprintf(msg, 1024, "TypeError: Expected *ALL* 3 floating point arrays to be the same type (allowed types = %d or %d). Instead found type-nums (%d, %d, %d)\n",
                   "with type-names = (%s, %s, %s)\n",
                   NPY_FLOAT, NPY_DOUBLE, x_type, y_type, z_type, x_descr->typeobj->tp_name, y_descr->typeobj->tp_name, z_descr->typeobj->tp_name);
        }
        Py_XDECREF(x_descr);Py_XDECREF(y_descr);Py_XDECREF(z_descr);
        countpairs_error_out(module, msg);
        return -1;
    }
    
    /* Check if the number of elements in the 3 Python arrays are identical */
    const int64_t nx1 = (int64_t)PyArray_SIZE(x1_obj);
    const int64_t ny1 = (int64_t)PyArray_SIZE(y1_obj);
    const int64_t nz1 = (int64_t)PyArray_SIZE(z1_obj);

    if(nx1 != ny1 || ny1 != nz1) {
      snprintf(msg, 1024, "ERROR: Expected arrays to have the same number of elements in all 3-dimensions.\nFound (nx, ny, nz) = (%"PRId64", %"PRId64", %"PRId64") instead",
               nx1, ny1, nz1);
      countpairs_error_out(module, msg);
      return -1;
    }


    /* Return the size of each element of the data object */
    if(x_type == NPY_FLOAT) {
      *element_size = sizeof(float);
    } else {
      *element_size = sizeof(double);
    }
    
    return nx1;
}


static PyObject *countpairs_countpairs(PyObject *self, PyObject *args)
{
    //Error-handling is global in python2 -> stored in struct module_state _struct declared at the top of this file
#if PY_MAJOR_VERSION < 3
    (void) self;
    PyObject *module = NULL;//should not be used -> setting to NULL so any attempts to dereference will result in a crash. 
#else
    //In python3, self is simply the module object that was returned earlier by init
    PyObject *module = self;
#endif    
    PyArrayObject *x1_obj, *y1_obj, *z1_obj, *x2_obj,*y2_obj,*z2_obj;
    int autocorr=0;
    int nthreads=4;
    char *binfile;

    /* Get the input args. */

    /* 
       Will throw a TypeError if the object inputs are not
       compatible with numpy ArrayType.
    */
    if (! PyArg_ParseTuple(args, "iisO!O!O!O!O!O!",&autocorr,&nthreads,&binfile,
                           &PyArray_Type,&x1_obj,
                           &PyArray_Type,&y1_obj,
                           &PyArray_Type,&z1_obj,
                           &PyArray_Type,&x2_obj,
                           &PyArray_Type,&y2_obj,
                           &PyArray_Type,&z2_obj)
        ) {
        Py_RETURN_NONE;
    }

    /* We have numpy arrays and all the required inputs*/
    /* How many data points are there? And are they all of floating point type */
    const int64_t ND1 = check_dims_and_datatype(module, x1_obj, y1_obj, z1_obj);
    if(ND1 == -1) {
        //Error has already been set -> simply return 
        Py_RETURN_NONE;
    }

    const int64_t ND2 = check_dims_and_datatype(module, x2_obj, y2_obj, z2_obj);
    if(ND2 == -1) {
        //Error has already been set -> simply return 
        Py_RETURN_NONE;
    }

    /* 
       Interpret the input objects as numpy arrays (of whatever the input type the python object has). 
       NULL initialization is necessary since we might be calling XDECREF.
       The input objects can be converted into the required DOUBLE array.
    */
    const int requirements = NPY_ARRAY_IN_ARRAY | NPY_ARRAY_FORCECAST;
    PyObject *x1_array = NULL, *y1_array = NULL, *z1_array = NULL;
    x1_array = PyArray_FromArray(x1_obj, ELEMENT_DESCR, requirements);
    y1_array = PyArray_FromArray(y1_obj, ELEMENT_DESCR, requirements);
    z1_array = PyArray_FromArray(z1_obj, ELEMENT_DESCR, requirements);
    
    /* NULL initialization is necessary since we might be calling XDECREF*/
    PyObject *x2_array = NULL, *y2_array = NULL, *z2_array = NULL;
    x2_array = PyArray_FromArray(x2_obj, ELEMENT_DESCR, requirements);
    y2_array = PyArray_FromArray(y2_obj, ELEMENT_DESCR, requirements);
    z2_array = PyArray_FromArray(z2_obj, ELEMENT_DESCR, requirements);
    
    if (x1_array == NULL || y1_array == NULL || z1_array == NULL ||
        x2_array == NULL || y2_array == NULL || z2_array == NULL) {
        Py_XDECREF(x1_array);
        Py_XDECREF(y1_array);
        Py_XDECREF(z1_array);

        Py_XDECREF(x2_array);
        Py_XDECREF(y2_array);
        Py_XDECREF(z2_array);
        char msg[1024];
        snprintf(msg, 1024, "TypeError: In %s: Could not convert to array of correct floating point type (need arrays of %s). Are you passing numpy arrays?",
                 __FUNCTION__, sizeof(DOUBLE) == 4 ? "floats":"doubles");
        countpairs_error_out(module, msg);
        Py_RETURN_NONE;
    }


    /* Get pointers to the data as DOUBLE C-types. */    
    DOUBLE *X1 = (DOUBLE *)PyArray_DATA((PyArrayObject *) x1_array); 
    DOUBLE *Y1 = (DOUBLE *)PyArray_DATA((PyArrayObject *) y1_array);
    DOUBLE *Z1 = (DOUBLE *)PyArray_DATA((PyArrayObject *) z1_array);
    
    DOUBLE *X2 = (DOUBLE *)PyArray_DATA((PyArrayObject *) x2_array);
    DOUBLE *Y2 = (DOUBLE *)PyArray_DATA((PyArrayObject *) y2_array);
    DOUBLE *Z2 = (DOUBLE *)PyArray_DATA((PyArrayObject *) z2_array);

    NPY_BEGIN_THREADS_DEF;
    NPY_BEGIN_THREADS;

    results_countpairs results = countpairs(ND1,X1,Y1,Z1,
                                            ND2,X2,Y2,Z2,
#if defined(USE_OMP) && defined(_OPENMP)
                                            nthreads,
#endif
                                            autocorr,
                                            binfile);
    NPY_END_THREADS;

    /* Clean up. */
    Py_DECREF(x1_array);Py_DECREF(y1_array);Py_DECREF(z1_array);
    Py_DECREF(x2_array);Py_DECREF(y2_array);Py_DECREF(z2_array);

    /* Build the output list */
    PyObject *ret = PyList_New(0);
    DOUBLE rlow=results.rupp[0];
    for(int i=1;i<results.nbin;i++) {
        PyObject *item = NULL;
        const DOUBLE rpavg = results.rpavg[i];

#ifdef DOUBLE_PREC
        item = Py_BuildValue("(dddk)", rlow,results.rupp[i],rpavg,results.npairs[i]);
#else
        item = Py_BuildValue("(fffk)", rlow,results.rupp[i],rpavg,results.npairs[i]);
#endif
        PyList_Append(ret, item);
        Py_XDECREF(item);
        rlow=results.rupp[i];
    }

    free_results(&results);
    return ret;
}


static PyObject *countpairs_countpairs_rp_pi(PyObject *self, PyObject *args)
{
#if PY_MAJOR_VERSION < 3
    (void) self;
    PyObject *module = NULL;//should not be used -> setting to NULL so any attempts to dereference will result in a crash. 
#else
    //In python3, self is simply the module object that was returned earlier by init
    PyObject *module = self;
#endif    
    PyArrayObject *x1_obj, *y1_obj, *z1_obj, *x2_obj,*y2_obj,*z2_obj;
    int autocorr=0;
    int nthreads=4;
    double pimax;
    char *binfile;

    if ( ! PyArg_ParseTuple(args, "iidsO!O!O!O!O!O!",
                            &autocorr,&nthreads,&pimax,&binfile,
                            &PyArray_Type,&x1_obj,
                            &PyArray_Type,&y1_obj,
                            &PyArray_Type,&z1_obj,
                            &PyArray_Type,&x2_obj,
                            &PyArray_Type,&y2_obj,
                            &PyArray_Type,&z2_obj)
         ) {
        Py_RETURN_NONE;
    }

    /* How many data points are there? And are they all of floating point type */
    const int64_t ND1 = check_dims_and_datatype(module, x1_obj, y1_obj, z1_obj);
    if(ND1 == -1) {
        //Error has already been set -> simply return 
        Py_RETURN_NONE;
    }

    const int64_t ND2 = check_dims_and_datatype(module, x2_obj, y2_obj, z2_obj);
    if(ND2 == -1) {
        //Error has already been set -> simply return 
        Py_RETURN_NONE;
    }

    /* Interpret the input objects as numpy arrays. */
    const int requirements = NPY_ARRAY_IN_ARRAY | NPY_ARRAY_FORCECAST ;
    PyObject *x1_array = PyArray_FromArray(x1_obj, ELEMENT_DESCR, requirements);
    PyObject *y1_array = PyArray_FromArray(y1_obj, ELEMENT_DESCR, requirements);
    PyObject *z1_array = PyArray_FromArray(z1_obj, ELEMENT_DESCR, requirements);

    PyObject *x2_array = PyArray_FromArray(x2_obj, ELEMENT_DESCR, requirements);
    PyObject *y2_array = PyArray_FromArray(y2_obj, ELEMENT_DESCR, requirements);
    PyObject *z2_array = PyArray_FromArray(z2_obj, ELEMENT_DESCR, requirements);

    if (x1_array == NULL || y1_array == NULL || z1_array == NULL ||
        x2_array == NULL || y2_array == NULL || z2_array == NULL) {
        Py_XDECREF(x1_array);
        Py_XDECREF(y1_array);
        Py_XDECREF(z1_array);

        Py_XDECREF(x2_array);
        Py_XDECREF(y2_array);
        Py_XDECREF(z2_array);
        char msg[1024];
        snprintf(msg, 1024, "TypeError: In %s: Could not convert to array of correct floating point type (need arrays of %s). Are you passing numpy arrays?",
                 __FUNCTION__, sizeof(DOUBLE) == 4 ? "floats":"doubles");
        countpairs_error_out(module, msg);
        Py_RETURN_NONE;
    }


    /* Get pointers to the data as C-types. */
    DOUBLE *X1 = (DOUBLE *)PyArray_DATA((PyArrayObject *) x1_array);
    DOUBLE *Y1 = (DOUBLE *)PyArray_DATA((PyArrayObject *) y1_array);
    DOUBLE *Z1 = (DOUBLE *)PyArray_DATA((PyArrayObject *) z1_array);

    DOUBLE *X2 = (DOUBLE *)PyArray_DATA((PyArrayObject *) x2_array);
    DOUBLE *Y2 = (DOUBLE *)PyArray_DATA((PyArrayObject *) y2_array);
    DOUBLE *Z2 = (DOUBLE *)PyArray_DATA((PyArrayObject *) z2_array);

    NPY_BEGIN_THREADS_DEF;
    NPY_BEGIN_THREADS;
    
    results_countpairs_rp_pi results = countpairs_rp_pi(ND1,X1,Y1,Z1,
                                                         ND2,X2,Y2,Z2,
#if defined(USE_OMP) && defined(_OPENMP)
                                                         nthreads,
#endif
                                                         autocorr,
                                                         binfile,
                                                         pimax);

    NPY_END_THREADS;
    
    /* Clean up. */
    Py_DECREF(x1_array);Py_DECREF(y1_array);Py_DECREF(z1_array);
    Py_DECREF(x2_array);Py_DECREF(y2_array);Py_DECREF(z2_array);

    /* Build the output list */
    PyObject *ret = PyList_New(0);//create an empty list
    DOUBLE rlow=results.rupp[0];
    const DOUBLE dpi = pimax/(DOUBLE)results.npibin ;

    for(int i=1;i<results.nbin;i++) {
        for(int j=0;j<results.npibin;j++) {
            const int bin_index = i*(results.npibin + 1) + j;
            PyObject *item = NULL;
            const DOUBLE rpavg = results.rpavg[bin_index];
#ifdef DOUBLE_PREC
            item = Py_BuildValue("(ddddk)", rlow,results.rupp[i],rpavg,(j+1)*dpi,results.npairs[bin_index]);
#else
            item = Py_BuildValue("(ffffk)", rlow,results.rupp[i],rpavg,(j+1)*dpi,results.npairs[bin_index]);
#endif
            PyList_Append(ret, item);
            Py_XDECREF(item);
        }
        rlow=results.rupp[i];
    }
    free_results_rp_pi(&results);
    return ret;
}

static PyObject *countpairs_countpairs_wp(PyObject *self, PyObject *args)
{
#if PY_MAJOR_VERSION < 3
    (void) self;//to suppress the unused variable warning. Terrible hack
    PyObject *module = NULL;//need not be used -> setting to NULL so any attempts to dereference will result in a crash. 
#else
    //In python3, self is simply the module object that was returned earlier by init
    PyObject *module = self;
#endif    
    PyArrayObject *x1_obj=NULL, *y1_obj=NULL, *z1_obj=NULL;
    double boxsize,pimax;
    int nthreads=1;
    char *binfile;
    size_t element_size;
    
    if( ! PyArg_ParseTuple(args, "ddisO!O!O!",&boxsize,&pimax,&nthreads,&binfile,
                           &PyArray_Type,&x1_obj,
                           &PyArray_Type,&y1_obj,
                           &PyArray_Type,&z1_obj)
        ){
        Py_RETURN_NONE;
    }
    
    /* How many data points are there? And are they all of floating point type */
    const int64_t ND1 = check_dims_and_datatype(module, x1_obj, y1_obj, z1_obj, &element_size);
    if(ND1 == -1) {
        //Error has already been set -> simply return 
        Py_RETURN_NONE;
    }
    
    /* Interpret the input objects as numpy arrays. */
    const int requirements = NPY_ARRAY_IN_ARRAY;
    PyObject *x1_array = PyArray_FromArray(x1_obj, NOTYPE_DESCR, requirements);
    PyObject *y1_array = PyArray_FromArray(y1_obj, NOTYPE_DESCR, requirements);
    PyObject *z1_array = PyArray_FromArray(z1_obj, NOTYPE_DESCR, requirements);
    
    if (x1_array == NULL || y1_array == NULL || z1_array == NULL) {
        Py_XDECREF(x1_array);
        Py_XDECREF(y1_array);
        Py_XDECREF(z1_array);
        char msg[1024];
        snprintf(msg, 1024, "TypeError: In %s: Could not convert to array of correct floating point type (need arrays of %s). Are you passing numpy arrays?",
                 __FUNCTION__, sizeof(DOUBLE) == 4 ? "floats":"doubles");
        perror(NULL);
        countpairs_error_out(module, msg);
        Py_RETURN_NONE;
    }


    /* Get pointers to the data as C-types. */
    void *X1 = PyArray_DATA(x1_array);
    void *Y1 = PyArray_DATA(y1_array);
    void *Z1 = PyArray_DATA(z1_array);

    NPY_BEGIN_THREADS_DEF;
    NPY_BEGIN_THREADS;

    
    results_countpairs_wp results;
    struct config_options options;
    options.need_avg_sep = 1;
    options.periodic = 1;
    options.float_type = element_size;
    int status = countpairs_wp(ND1,X1,Y1,Z1,
                               boxsize,
                               nthreads,
                               binfile,
                               pimax,
                               &results,
                               &options);
    
    NPY_END_THREADS;

    /* Clean up. */
    Py_DECREF(x1_array);Py_DECREF(y1_array);Py_DECREF(z1_array);

    if(status != EXIT_SUCCESS) {
        Py_RETURN_NONE;
    }

    
#if 0
    for(int i=1;i<results.nbin;i++) {
        const DOUBLE rpavg = results.rpavg[i];
        fprintf(stderr,"%lf %lf %lf %lf %"PRIu64"\n",results.rupp[i-1],results.rupp[i],rpavg,results.wp[i],results.npairs[i]);
    }
#endif

    /* Build the output list */
    PyObject *ret = PyList_New(0);
    DOUBLE rlow=results.rupp[0];
    for(int i=1;i<results.nbin;i++) {
        PyObject *item = NULL;
        const DOUBLE rpavg = results.rpavg[i];
        item = Py_BuildValue("(ddddk)", rlow,results.rupp[i],rpavg,results.wp[i],results.npairs[i]);
        PyList_Append(ret, item);
        Py_XDECREF(item);
        rlow=results.rupp[i];
    }
    free_results_wp(&results);
    return ret;
}


static PyObject *countpairs_countpairs_xi(PyObject *self, PyObject *args)
{
#if PY_MAJOR_VERSION < 3
    (void) self;//to suppress the unused variable warning. Terrible hack
    PyObject *module = NULL;//should not be used -> setting to NULL so any attempts to dereference will result in a crash. 
#else
    //In python3, self is simply the module object that was returned earlier by init
    PyObject *module = self;
#endif    

    PyArrayObject *x1_obj, *y1_obj, *z1_obj;
    double boxsize;
    int nthreads=4;
    char *binfile;

    if( ! PyArg_ParseTuple(args, "disO!O!O!",&boxsize,&nthreads,&binfile,
                           &PyArray_Type,&x1_obj,
                           &PyArray_Type,&y1_obj,
                           &PyArray_Type,&z1_obj)
        ) {
        Py_RETURN_NONE;
    }

    /* How many data points are there? And are they all of floating point type */
    const int64_t ND1 = check_dims_and_datatype(module, x1_obj, y1_obj, z1_obj);
    if(ND1 == -1) {
        //Error has already been set -> simply return 
        Py_RETURN_NONE;
    }

    /* Interpret the input objects as numpy arrays. */
    const int requirements = NPY_ARRAY_IN_ARRAY | NPY_ARRAY_FORCECAST ;    
    PyObject *x1_array = PyArray_FromArray(x1_obj, ELEMENT_DESCR, requirements);
    PyObject *y1_array = PyArray_FromArray(y1_obj, ELEMENT_DESCR, requirements);
    PyObject *z1_array = PyArray_FromArray(z1_obj, ELEMENT_DESCR, requirements);

    if (x1_array == NULL || y1_array == NULL || z1_array == NULL) {
        Py_XDECREF(x1_array);
        Py_XDECREF(y1_array);
        Py_XDECREF(z1_array);
        char msg[1024];
        snprintf(msg, 1024, "TypeError: In %s: Could not convert to array of correct floating point type (need arrays of %s). Are you passing numpy arrays?",
                 __FUNCTION__, sizeof(DOUBLE) == 4 ? "floats":"doubles");
        countpairs_error_out(module, msg);
        Py_RETURN_NONE;
    }

    /* Get pointers to the data as C-types. */
    DOUBLE *X1 = (DOUBLE *)PyArray_DATA((PyArrayObject *) x1_array);
    DOUBLE *Y1 = (DOUBLE *)PyArray_DATA((PyArrayObject *) y1_array);
    DOUBLE *Z1 = (DOUBLE *)PyArray_DATA((PyArrayObject *) z1_array);

    NPY_BEGIN_THREADS_DEF;
    NPY_BEGIN_THREADS;

    results_countpairs_xi results = countpairs_xi(ND1,X1,Y1,Z1,
                                                  boxsize,
#if defined(USE_OMP) && defined(_OPENMP)
                                                  nthreads,
#endif
                                                  binfile);
    NPY_END_THREADS;

    /* Clean up. */
    Py_DECREF(x1_array);Py_DECREF(y1_array);Py_DECREF(z1_array);

#if 0
    for(int i=1;i<results.nbin;i++) {
        const DOUBLE rpavg = results.rpavg[i];
        fprintf(stderr,"%lf %lf %lf %lf %"PRIu64"\n",results.rupp[i-1],results.rupp[i],rpavg,results.xi[i],results.npairs[i]);
    }
#endif

    /* Build the output list */
    PyObject *ret = PyList_New(0);
    DOUBLE rlow=results.rupp[0];
    for(int i=1;i<results.nbin;i++) {
        PyObject *item = NULL;
        const DOUBLE rpavg = results.rpavg[i];

#ifdef DOUBLE_PREC
        item = Py_BuildValue("(ddddk)", rlow,results.rupp[i],rpavg,results.xi[i],results.npairs[i]);
#else
        item = Py_BuildValue("(ffffk)", rlow,results.rupp[i],rpavg,results.xi[i],results.npairs[i]);
#endif//DOUBLE_PREC

        PyList_Append(ret, item);
        Py_XDECREF(item);
        rlow=results.rupp[i];
    }
    free_results_xi(&results);
    return ret;
}

static PyObject *countpairs_countspheres_vpf(PyObject *self, PyObject *args)
{
#if PY_MAJOR_VERSION < 3
    (void) self;//to suppress the unused variable warning. Terrible hack
    PyObject *module = NULL;//should not be used -> setting to NULL so any attempts to dereference will result in a crash. 
#else
    //In python3, self is simply the module object that was returned earlier by init
    PyObject *module = self;
#endif    

    PyArrayObject *x1_obj, *y1_obj, *z1_obj;
    double rmax;
    int nbin,nc,num_pN;
    unsigned long seed=-1;

    if( ! PyArg_ParseTuple(args, "diiikO!O!O!",&rmax,&nbin,&nc,&num_pN,&seed,
                           &PyArray_Type,&x1_obj,
                           &PyArray_Type,&y1_obj,
                           &PyArray_Type,&z1_obj)
        ) {
        Py_RETURN_NONE;
    }
    
    /* How many data points are there? And are they all of floating point type */
    const int64_t ND1 = check_dims_and_datatype(module, x1_obj, y1_obj, z1_obj);
    if(ND1 == -1) {
        //Error has already been set -> simply return 
        Py_RETURN_NONE;
    }

    /* Interpret the input objects as numpy arrays. */
    const int requirements = NPY_ARRAY_IN_ARRAY | NPY_ARRAY_FORCECAST ;        
    PyObject *x1_array = PyArray_FromArray(x1_obj, ELEMENT_DESCR, requirements);
    PyObject *y1_array = PyArray_FromArray(y1_obj, ELEMENT_DESCR, requirements);
    PyObject *z1_array = PyArray_FromArray(z1_obj, ELEMENT_DESCR, requirements);

    if (x1_array == NULL || y1_array == NULL || z1_array == NULL) {
        Py_XDECREF(x1_array);
        Py_XDECREF(y1_array);
        Py_XDECREF(z1_array);
        char msg[1024];
        snprintf(msg, 1024, "TypeError: In %s: Could not convert to array of correct floating point type (need arrays of %s). Are you passing numpy arrays?",
                 __FUNCTION__, sizeof(DOUBLE) == 4 ? "floats":"doubles");
        countpairs_error_out(module, msg);
        Py_RETURN_NONE;
    }

    /* Get pointers to the data as C-types. */
    DOUBLE *X1 = (DOUBLE *)PyArray_DATA((PyArrayObject *) x1_array);
    DOUBLE *Y1 = (DOUBLE *)PyArray_DATA((PyArrayObject *) y1_array);
    DOUBLE *Z1 = (DOUBLE *)PyArray_DATA((PyArrayObject *) z1_array);

    /* Do the VPF calculation */
    results_countspheres results = countspheres(ND1, X1, Y1, Z1,
                                                rmax, nbin, nc,
                                                num_pN,
                                                seed);

    /* Clean up. */
    Py_DECREF(x1_array);Py_DECREF(y1_array);Py_DECREF(z1_array);

    /* Build the output list (of lists, since num_pN is determined at runtime) */
    PyObject *ret = PyList_New(0);
    const DOUBLE rstep = rmax/(DOUBLE)nbin ;
    for(int ibin=0;ibin<results.nbin;ibin++) {
        const double r=(ibin+1)*rstep;
        PyObject *item = PyList_New(0);
        PyObject *this_val = Py_BuildValue("d",r);
        PyList_Append(item, this_val);
        Py_XDECREF(this_val);
        for(int i=0;i<num_pN;i++) {
#ifdef DOUBLE_PREC
            this_val = Py_BuildValue("d",(results.pN)[ibin][i]);
#else
            this_val = Py_BuildValue("d",(results.pN)[ibin][i]);
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
