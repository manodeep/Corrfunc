/* File: _countpairs_mocks.c */
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
#include "countpairs_rp_pi_mocks.h"
#include "countpairs_theta_mocks.h"

//for the vpf
#include "countspheres_mocks.h"

struct module_state {
    PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
//python3 follows
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#define INITERROR return NULL
PyObject *PyInit__countpairs_mocks(void);

#else
//python2 follows
#define GETSTATE(m) (&_state)
static struct module_state _state;
#define INITERROR return
PyMODINIT_FUNC init_countpairs_mocks(void);

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
/* #define THETA_UNICODE "\u03B8" */
#define THETA_UNICODE "theta"
/* #define RP_UNICODE    "r\u209a" */
/* #endif */


//Docstrings for the methods
static char module_docstring[] =    "This module provides an interface for calculating correlation functions on MOCKS (spherical geometry) using C.\n"
    "\n"
    "countpairs_rp_pi_mocks: Calculate the 2-D DD(rp,pi) auto/cross-correlation function given two sets of ra/dec/cz and ra/dec/cz arrays.\n"
    "countpairs_theta_mocks: Calculate DD(theta) auto/cross-correlation function given two sets of ra/dec/cz and ra/dec/cz arrays.\n"
    "countspheres_vpf_mocks: Calculate the counts-in-spheres given one set of ra/dec/cz.\n"
    "\n\n"
    "See `Corrfunc/call_correlation_functions_mocks.py` for example calls to each function.\n";

static char error_out_docstring[]          =  "Error-handler for the module.";

/* function proto-type*/
static PyObject *countpairs_countpairs_rp_pi_mocks(PyObject *self, PyObject *args);
static PyObject *countpairs_countpairs_theta_mocks(PyObject *self, PyObject *args);
static PyObject *countpairs_countspheres_vpf_mocks(PyObject *self, PyObject *args);
static PyObject *countpairs_mocks_error_out(PyObject *module);

static PyMethodDef module_methods[] = {
    {"countpairs_mocks_error_out"   ,(PyCFunction) countpairs_mocks_error_out        ,METH_NOARGS, error_out_docstring},
    {"countpairs_rp_pi_mocks"       ,(PyCFunction) countpairs_countpairs_rp_pi_mocks ,METH_VARARGS,
     "countpairs_rp_pi_mocks(autocorr, cosmology, nthreads, "PIMAX_UNICODE", binfile, RA1, DEC1, CZ1, RA2, DEC2, CZ2)\n"
     "\n"
     "Calculate the 2-D pair-counts, "XI_UNICODE"(rp, "PI_UNICODE"), auto/cross-correlation function given two\n"
     "sets of RA1/DEC1/CZ1 and RA2/DEC2/CZ2 arrays. This module is suitable for mock catalogs that have been\n"
     "created by carving out a survey footprint from simulated data. The module can also be used for actual\n"
     "observed galaxies, but you probably want to attach weights to the points to account for completeness etc.\n"
     "(i.e., seriously consider using some other code for real survey data). Default bins in "PI_UNICODE"\n"
     "are set to 1.0 Mpc/h.\n"
     "\n"
     "Note, that this module only returns pair counts and not the actual correlation function\n"
     ""XI_UNICODE"(rp, "PI_UNICODE"). See the xi_mocks/DDrppi/wprp_mocks.c for computing wp(rp) from DD(rp, "PI_UNICODE").\n"
     "\n"
     "parameters (all parameters are required)\n"
     "----------------------------------------\n"
     "autocorr: boolean\n"
     "    Flag for auto/cross-correlation. If autocorr is not 0, the RA2/DEC2/CZ2 arrays\n"
     "    are not used (but must still be passed, as RA1/DEC1/CZ1).\n"
     "cosmology: integer\n"
     "    Integer to select cosmology. Pre-set values for (1,2) \n"
     "    1 -> LasDamas cosmology. Om=0.25,  Ol=0.75  (other values are not used)\n"
     "    2 -> Planck   cosmology. Om=0.302, Ol=0.698 \n"
     "    To setup a new cosmology, add an entry to the function, `init_cosmology` in \n"
     "    `ROOT/utils/cosmology_params.c` and recompile the package.\n"
     "nthreads: integer\n"
     "    The number of OpenMP threads to use. Has no effect if OpenMP was not used\n"
     "    during library compilation. \n"
     ""PIMAX_UNICODE": double (Mpc/h)\n"
     "    The max. integration distance along the "PI_UNICODE" direction in Mpc/h. Typical\n"
     "    values are in the range 40-100 Mpc/h.\n"
     "binfile: filename\n"
     "    Filename containing the radial bins for the correlation function. The file\n"
     "    is expected to contain white-space separated ``rmin  rmax`` with the bin\n"
     "    edges.  Units must be Mpc/h (see the ``bins`` file in the tests directory\n"
     "    for a sample). For usual logarithmic bins, ``logbins``in the root directory\n"
     "    of this package will create a compatible ``binfile``.\n"
     "RA1: float/double (default double)\n"
     "    The right-ascension of the galaxy, in the range [0, 360]. If there are\n"
     "    negative RA's in the supplied array (input RA in the range [-180, 180]),\n"
     "    then the code will shift the entire array by 180 to put RA's in the\n"
     "    [0, 360] range.\n"
     "    Float/double must be consistent with the library otherwise there will\n"
     "    be a RunTimeError from numpy. Specifically, if ``DOUBLE_PREC`` is not\n"
     "    enabled, then the arrays should be compatible with ``numpy.float32``;\n"
     "    otherwise, use ``numpy.float``.\n"
     "DEC1: float/double (default double)\n"
     "    The declination of the galaxy, in the range [-90, 90]. If there are\n"
     "    declinations > 90 in the supplied array (input dec in the range [0, 180]),\n"
     "    then the code will shift the entire array by -90 to put declinations in\n"
     "    the [-90, 90] range. If the code finds declinations more than 180, then\n"
     "    it assumes RA and DEC have been swapped and aborts with that message.\n"
     "    Float/double must be consistent with the library otherwise there will\n"
     "    be a RunTimeError from numpy. Specifically, if ``DOUBLE_PREC`` is not\n"
     "    enabled, then the arrays should be compatible with ``numpy.float32``;\n"
     "    otherwise, use ``numpy.float``.\n"
     "CZ1: float/double (default double)\n"
     "    The redshift multiplied by speed of light for the galaxies. The code will\n"
     "    checks that cz has been supplied by comparing with a threshold (currently\n"
     "    set to 10, defined in function check_ra_dec_cz in file\n"
     "    `DDrppi/countpairs_rp_pi_mocks.c`) and multiplies by the speed of light if\n"
     "    max z is less than that threshold. If you really want to change the speed\n"
     "    of light, then edit the macro in `ROOT/utils/set_cosmo_dist.h`.\n"
     "RA2/DEC2/CZ2: float/double (default double)\n"
     "    Same as for RA1/DEC1/CZ1\n"
     "\n"
     "returns\n"
     "-------\n"
     "a Python list containing [rmin, rmax, ravg, "PI_UNICODE", npairs] \n"
     "for each "PI_UNICODE"-bin (up to "PIMAX_UNICODE") for each radial bin specified in\n"
     "the ``binfile``. For instance, for a ``"PIMAX_UNICODE"`` of 40.0 Mpc/h, each radial\n"
     "bin will be split into 40 "PI_UNICODE" bins (default "PI_UNICODE" bin is 1.0). Thus, the\n"
     "total number of items in the list is {(int) ``"PIMAX_UNICODE"`` * number of rp bins}.\n"
     "If ``OUTPUT_RPAVG`` is not defined in ``mocks.options`` \n"
     "then ``ravg`` will be set to 0.0 for all bins. "PI_UNICODE" for each bin\n"
     "is the upper limit of the "PI_UNICODE" values that were considered in that (rp, "PI_UNICODE") bin.\n"
     "``npairs``contains the number of pairs in that bin and can be used to compute the\n"
     "actual "XI_UNICODE"(rp, "PI_UNICODE") by combining with RR counts.\n"
     "(The C struct is identical to the one in `xi_theory/xi_rp_pi/countpairs_rp_pi.h`)\n"
     "\n"
     "example\n"
     "-------\n"
     "import numpy as np\n"
     "from Corrfunc._countpairs_mocks import countpairs_rp_pi_mocks\n"
     "ra,dec,cz = np.genfromtxt('../xi_mocks/tests/data/Mr19_mock_northonly.rdcz.dat',dtype=np.float,unpack=True)\n"
     "cosmology=1\n"
     "autocorr=1\n"
     "nthreads=4\n"
     "binfile='../xi_mocks/tests/bins'\n"
     "pimax=40.0\n"
     "DD=countpairs_rp_pi_mocks(autocorr, cosmology,nthreads,pimax,binfile,ra,dec,cz,ra,dec,cz)\n"
     "\n"
     "See `Corrfunc/call_correlation_functions_mocks.py`\n"
    },
    {"countpairs_theta_mocks"       ,(PyCFunction) countpairs_countpairs_theta_mocks ,METH_VARARGS,
     "countpairs_theta_mocks(autocorr, cosmology, nthreads, binfile, RA1, DEC1, RA2, DEC2)\n"
     "\n"
     "Calculate the angular pair-counts, "XI_UNICODE"(rp, "THETA_UNICODE"), auto/cross-correlation function given two\n"
     "sets of RA1/DEC1 and RA2/DEC2 arrays. This module is suitable for mock catalogs that have been\n"
     "created by carving out a survey footprint from simulated data. The module can also be used for actual\n"
     "observed galaxies, but you probably want to attach weights to the points to account for completeness etc.\n"
     "(i.e., seriously consider using some other code for real survey data).\n"
     "\n"
     "Note, that this module only returns pair counts and not the actual correlation function\n"
     ""XI_UNICODE"("THETA_UNICODE"). See the xi_mocks/wtheta/wtheta.c for computing w(theta) from DD("THETA_UNICODE").\n"
     "\n"
     "parameters (all parameters are required)\n"
     "----------------------------------------\n"
     "autocorr: boolean\n"
     "    Flag for auto/cross-correlation. If autocorr is not 0, the RA2/DEC2 arrays\n"
     "    are not used (but must still be passed as valid arrays).\n"
     "cosmology: integer\n"
     "    Integer to select cosmology. Pre-set values for (1,2) \n"
     "    1 -> LasDamas cosmology. Om=0.25,  Ol=0.75  (other values are not used)\n"
     "    2 -> Planck   cosmology. Om=0.302, Ol=0.698 \n"
     "    To setup a new cosmology, add an entry to the function, `init_cosmology` in \n"
     "    `ROOT/utils/cosmology_params.c` and recompile the package.\n"
     "nthreads: integer\n"
     "    The number of OpenMP threads to use. Has no effect if OpenMP was not used\n"
     "    during library compilation. \n"
     "binfile: filename\n"
     "    Filename containing the radial bins for the correlation function. The file\n"
     "    is expected to contain white-space separated ``thetamin  thetamax`` with the bin\n"
     "    edges.  Units must be degrees (see the ``angular_bins`` file in the tests directory\n"
     "    for a sample). For usual logarithmic bins, ``logbins``in the root directory\n"
     "    of this package will create a compatible ``binfile``.\n"
     "RA1: float/double (default double)\n"
     "    The right-ascension of the galaxy, in the range [0, 360]. If there are\n"
     "    negative RA's in the supplied array (input RA in the range [-180, 180]),\n"
     "    then the code will shift the entire array by 180 to put RA's in the\n"
     "    [0, 360] range.\n"
     "    Float/double must be consistent with the library otherwise there will\n"
     "    be a RunTimeError from numpy. Specifically, if ``DOUBLE_PREC`` is not\n"
     "    enabled, then the arrays should be compatible with ``numpy.float32``;\n"
     "    otherwise, use ``numpy.float``.\n"
     "DEC1: float/double (default double)\n"
     "    The declination of the galaxy, in the range [-90, 90]. If there are\n"
     "    declinations > 90 in the supplied array (input dec in the range [0, 180]),\n"
     "    then the code will shift the entire array by -90 to put declinations in\n"
     "    the [-90, 90] range. If the code finds declinations more than 180, then\n"
     "    it assumes RA and DEC have been swapped and aborts with that message.\n"
     "    Float/double must be consistent with the library otherwise there will\n"
     "    be a RunTimeError from numpy. Specifically, if ``DOUBLE_PREC`` is not\n"
     "    enabled, then the arrays should be compatible with ``numpy.float32``;\n"
     "    otherwise, use ``numpy.float``.\n"
     "RA2/DEC2: float/double (default double)\n"
     "    Same as for RA1/DEC1\n"
     "\n"
     "returns\n"
     "-------\n"
     "a Python list containing [thetamin, thetamax, thetaavg, npairs] for each angular\n"
     "bin specified in the ``binfile``. If ``OUTPUT_THETAAVG`` is not defined in\n"
     "``mocks.options`` then ``thetaravg`` will be set to 0.0 for all bins. ``npairs``\n"
     "contains the number of pairs in that bin and can be used to compute the\n"
     "actual "XI_UNICODE"("THETA_UNICODE") by combining  DD, (DR) and RR counts.\n"
     "\n"
     "example\n"
     "-------\n"
     "import numpy as np\n"
     "from Corrfunc._countpairs_mocks import countpairs_theta_mocks\n"
     "ra,dec,cz = np.genfromtxt('../xi_mocks/tests/data/Mr19_mock_northonly.rdcz.dat',dtype=np.float,unpack=True)\n"
     "cosmology=1\n"
     "autocorr=1\n"
     "nthreads=4\n"
     "binfile='../xi_mocks/tests/angular_bins'\n"
     "pimax=40.0\n"
     "DD=countpairs_theta_mocks(autocorr, cosmology, nthreads, binfile, ra,dec, ra,dec)\n"
     "\n"
     "See `Corrfunc/call_correlation_functions_mocks.py`\n"
    },
    {"countspheres_vpf_mocks"       ,(PyCFunction) countpairs_countspheres_vpf_mocks ,METH_VARARGS,
     /* if (!PyArg_ParseTuple(args, "diiiisiOOOOOO",&rmax,&nbin,&num_spheres,&num_pN,&threshold_neighbors,&centers_file,&cosmology,&x1_obj,&y1_obj,&z1_obj,&x2_obj,&y2_obj,&z2_obj)) */
     "countspheres_vpf_mocks(rmax, nbin, ncenters, num_pN, threshold_ngb, centers_file, cosmology, RA1, DEC1, CZ1, RA2, DEC2, CZ2)\n"
     "\n"
     "Calculates the fraction of random spheres that contain exactly *N* points, pN(r).\n"
     "\n"
     "----------------------------------------\n"
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
     "threshold_ngb: integer\n"
     "    Minimum number of random points needed in a `rmax` sphere such that it is\n"
     "    considered to be entirely within the mock footprint. `xi_mocks/vpf/vpf_mocks.c`\n"
     "    assumes that the minimum number of randoms can be at most a 1-sigma deviation from\n"
     "    the expected random number density.\n"
     "centers_file: filename\n"
     "    A file containing random sphere centers. If the file does not exist, then\n"
     "    a list of random centers will be written out. In that case, the randoms arrays,\n"
     "    RA2/DEC2/CZ2 are used to check that the sphere is entirely within the footprint.\n"
     "    If the file does exist but either `rmax' is too small or there are not enough centers\n"
     "    then the file will be overwritten. If the centers file has to be written, the code\n"
     "    will take significantly longer to finish (however, subsequent runs can re-use that\n"
     "    centers file and will be faster).\n"
     "cosmology: integer\n"
     "    Integer to select cosmology. Pre-set values for (1,2) \n"
     "    1 -> LasDamas cosmology. Om=0.25,  Ol=0.75  (other values are not used)\n"
     "    2 -> Planck   cosmology. Om=0.302, Ol=0.698 \n"
     "    To setup a new cosmology, add an entry to the function, `init_cosmology` in \n"
     "    `ROOT/utils/cosmology_params.c` and recompile the package.\n"
     "RA1: float/double (default double)\n"
     "    The right-ascension of the galaxy, in the range [0, 360]. If there are\n"
     "    negative RA's in the supplied array (input RA in the range [-180, 180]),\n"
     "    then the code will shift the entire array by 180 to put RA's in the\n"
     "    [0, 360] range.\n"
     "    Float/double must be consistent with the library otherwise there will\n"
     "    be a RunTimeError from numpy. Specifically, if ``DOUBLE_PREC`` is not\n"
     "    enabled, then the arrays should be compatible with ``numpy.float32``;\n"
     "    otherwise, use ``numpy.float``.\n"
     "DEC1: float/double (default double)\n"
     "    The declination of the galaxy, in the range [-90, 90]. If there are\n"
     "    declinations > 90 in the supplied array (input dec in the range [0, 180]),\n"
     "    then the code will shift the entire array by -90 to put declinations in\n"
     "    the [-90, 90] range. If the code finds declinations more than 180, then\n"
     "    it assumes RA and DEC have been swapped and aborts with that message.\n"
     "    Float/double must be consistent with the library otherwise there will\n"
     "    be a RunTimeError from numpy. Specifically, if ``DOUBLE_PREC`` is not\n"
     "    enabled, then the arrays should be compatible with ``numpy.float32``;\n"
     "    otherwise, use ``numpy.float``.\n"
     "CZ1: float/double (default double)\n"
     "    The redshift multiplied by speed of light for the galaxies. The code will\n"
     "    checks that cz has been supplied by comparing with a threshold (currently\n"
     "    set to 10, defined in function check_ra_dec_cz in file\n"
     "    `DDrppi/countpairs_rp_pi_mocks.c`) and multiplies by the speed of light if\n"
     "    max z is less than that threshold. If you really want to change the speed\n"
     "    of light, then edit the macro in `ROOT/utils/set_cosmo_dist.h`.\n"
     "RA2/DEC2/CZ2: float/double (default double)\n"
     "    Same as for RA1/DEC1/CZ1 but must be a random distribution of points on the\n"
     "    same footprint as RA1/DEC1/CZ1. If the centers file contains enough random\n"
     "    centre (at least `ncenters`), then these arrays are not accessed. This is a\n"
     "    very important distinction from all the other modules in this package.\n"
     "\n"
     "returns\n"
     "-------\n"
     "a Python list containing [rmax, p0, p1,..., p(num_pN-1)] for each radial bin.\n"
    },
    {NULL, NULL, 0, NULL}
};


static PyObject *countpairs_mocks_error_out(PyObject *module) {
    (void) module;//to avoid unused warning
    struct module_state *st = GETSTATE(module);
    PyErr_SetString(st->error, "ERROR: (generic) something bad happened");
    return NULL;
}


#if PY_MAJOR_VERSION >= 3

static int _countpairs_mocks_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int _countpairs_mocks_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}


static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_countpairs_mocks",
    module_docstring,
    sizeof(struct module_state),
    module_methods,
    NULL,
    _countpairs_mocks_traverse,
    _countpairs_mocks_clear,
    NULL
};


PyObject *PyInit__countpairs_mocks(void)

#else
    PyMODINIT_FUNC init_countpairs_mocks(void)
#endif
{

#if PY_MAJOR_VERSION >= 3
    PyObject *module = PyModule_Create(&moduledef);
#else
    PyObject *module = Py_InitModule3("_countpairs_mocks", module_methods, module_docstring);
#endif

    if (module == NULL) {
        INITERROR;
    }

    struct module_state *st = GETSTATE(module);

    st->error = PyErr_NewException("_countpairs_mocks.error", NULL, NULL);
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
    DOUBLE *phiD1   = (DOUBLE *)PyArray_DATA((PyArrayObject *) x1_array);
    DOUBLE *thetaD1 = (DOUBLE *)PyArray_DATA((PyArrayObject *) y1_array);
    DOUBLE *czD1    = (DOUBLE *)PyArray_DATA((PyArrayObject *) z1_array);

    DOUBLE *phiD2   = (DOUBLE *)PyArray_DATA((PyArrayObject *) x2_array);
    DOUBLE *thetaD2 = (DOUBLE *)PyArray_DATA((PyArrayObject *) y2_array);
    DOUBLE *czD2    = (DOUBLE *)PyArray_DATA((PyArrayObject *) z2_array);


    results_countpairs_mocks *results  = countpairs_mocks(ND1,phiD1,thetaD1,czD1,
                                                          ND2,phiD2,thetaD2,czD2,
#if defined(USE_OMP) && defined(_OPENMP)
                                                          nthreads,
#endif
                                                          autocorr,
                                                          binfile,
                                                          pimax,
                                                          cosmology);

    /* Clean up. */
    Py_DECREF(x1_array);Py_DECREF(y1_array);Py_DECREF(z1_array);
    Py_DECREF(x2_array);Py_DECREF(y2_array);Py_DECREF(z2_array);

#if 0
    /* Output pairs*/
    for(int i=1;i<results->nbin;i++) {
        const double logrp = LOG10(results->rupp[i]);
        for(int j=0;j<npibin;j++) {
            const int index = i*(npibin+1) + j;
            fprintf(stdout,"%10"PRIu64" %20.8lf %20.8lf  %20.8lf \n",results->npairs[index],results->rpavg[index],logrp,(j+1)*dpi);
        }
    }
#endif


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
    PyObject *x1_array = PyArray_FROM_OTF(x1_obj, ELEMENT_SIZE, NPY_ARRAY_IN_ARRAY);
    PyObject *y1_array = PyArray_FROM_OTF(y1_obj, ELEMENT_SIZE, NPY_ARRAY_IN_ARRAY);

    if (x1_array == NULL || y1_array == NULL) {
        Py_XDECREF(x1_array);
        Py_XDECREF(y1_array);
        return NULL;
    }

    /* Interpret the input objects as numpy arrays. */
    PyObject *x2_array = PyArray_FROM_OTF(x2_obj, ELEMENT_SIZE, NPY_ARRAY_IN_ARRAY);
    PyObject *y2_array = PyArray_FROM_OTF(y2_obj, ELEMENT_SIZE, NPY_ARRAY_IN_ARRAY);

    if (x2_array == NULL || y2_array == NULL) {
        Py_XDECREF(x2_array);
        Py_XDECREF(y2_array);
        return NULL;
    }


    /* How many data points are there? */
    const int64_t ND1 = (int64_t)PyArray_DIM((const PyArrayObject *) x1_array, 0);
    const int64_t ND2 = (int64_t)PyArray_DIM((const PyArrayObject *) x2_array, 0);

    /* Get pointers to the data as C-types. */
    DOUBLE *phiD1   = (DOUBLE *)PyArray_DATA((PyArrayObject *) x1_array);
    DOUBLE *thetaD1 = (DOUBLE *)PyArray_DATA((PyArrayObject *) y1_array);

    DOUBLE *phiD2   = (DOUBLE *)PyArray_DATA((PyArrayObject *) x2_array);
    DOUBLE *thetaD2 = (DOUBLE *)PyArray_DATA((PyArrayObject *) y2_array);


    results_countpairs_theta *results = countpairs_theta_mocks(ND1,phiD1,thetaD1,
                                                               ND2,phiD2,thetaD2,
#if defined(USE_OMP) && defined(_OPENMP)
                                                               nthreads,
#endif
                                                               autocorr,
                                                               binfile) ;


    /* Clean up. */
    Py_DECREF(x1_array);Py_DECREF(y1_array);
    Py_DECREF(x2_array);Py_DECREF(y2_array);

#if 0
    /*---Output-Pairs-------------------------------------*/
    DOUBLE theta_low = results->theta_upp[0];
    for(int i=1;i<results->nbin;i++) {
        fprintf(stdout,"%10"PRIu64" %20.8lf %20.8lf %20.8lf \n",results->npairs[i],results->theta_avg[i],theta_low,results->theta_upp[i]);
        theta_low=results->theta_upp[i];
    }
#endif

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
    DOUBLE *phiD1   = (DOUBLE *)PyArray_DATA((PyArrayObject *) x1_array);
    DOUBLE *thetaD1 = (DOUBLE *)PyArray_DATA((PyArrayObject *) y1_array);
    DOUBLE *czD1    = (DOUBLE *)PyArray_DATA((PyArrayObject *) z1_array);

    DOUBLE *phiD2   = (DOUBLE *)PyArray_DATA((PyArrayObject *) x2_array);
    DOUBLE *thetaD2 = (DOUBLE *)PyArray_DATA((PyArrayObject *) y2_array);
    DOUBLE *czD2    = (DOUBLE *)PyArray_DATA((PyArrayObject *) z2_array);

    results_countspheres_mocks *results = countspheres_mocks(ND1, phiD1,thetaD1, czD1,
                                                             ND2, phiD2,thetaD2, czD2,
                                                             threshold_neighbors,
                                                             rmax, nbin, num_spheres,
                                                             num_pN,
                                                             centers_file,
                                                             cosmology);
#if 0
    // Output the results
    const DOUBLE rstep = rmax/(DOUBLE)nbin ;
    for(int ibin=0;ibin<results->nbin;ibin++) {
        const double r=(ibin+1)*rstep;
        fprintf(stdout,"%10.2"DOUBLE_FORMAT" ", r);
        for(int i=0;i<num_pN;i++) {
            fprintf(stdout," %10.4e", (results->pN)[ibin][i]);
        }
        fprintf(stdout,"\n");
    }
#endif

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
