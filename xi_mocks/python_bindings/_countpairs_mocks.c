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

#define NOTYPE_DESCR     (PyArray_DescrFromType(NPY_NOTYPE))

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
static PyObject *countpairs_countpairs_rp_pi_mocks(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *countpairs_countpairs_theta_mocks(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *countpairs_countspheres_vpf_mocks(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *countpairs_mocks_error_out(PyObject *module, const char *msg);

static PyMethodDef module_methods[] = {
    {"countpairs_mocks_error_out"   ,(PyCFunction) countpairs_mocks_error_out        ,METH_VARARGS, error_out_docstring},
    {"countpairs_rp_pi_mocks"       ,(PyCFunction) countpairs_countpairs_rp_pi_mocks ,METH_VARARGS | METH_KEYWORDS,
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
    {"countpairs_theta_mocks"       ,(PyCFunction) countpairs_countpairs_theta_mocks ,METH_VARARGS | METH_KEYWORDS,
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
    {"countspheres_vpf_mocks"       ,(PyCFunction) countpairs_countspheres_vpf_mocks ,METH_VARARGS | METH_KEYWORDS,
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


static PyObject *countpairs_mocks_error_out(PyObject *module, const char *msg)
{
#if PY_MAJOR_VERSION < 3
    (void) module;//to avoid unused warning with python2
#endif
    
    struct module_state *st = GETSTATE(module);
    PyErr_SetString(st->error, msg);
    Py_RETURN_NONE;
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

static int64_t check_dims_and_datatype(PyObject *module, PyArrayObject *x1_obj, PyArrayObject *y1_obj, PyArrayObject *z1_obj, size_t *element_size)
{
    char msg[1024];

    /* All the arrays should be 1-D*/
    const int nxdims = PyArray_NDIM(x1_obj);
    const int nydims = PyArray_NDIM(y1_obj);
    const int nzdims = PyArray_NDIM(z1_obj);

    if(nxdims != 1 || nydims != 1 || nzdims != 1) {
        snprintf(msg, 1024, "ERROR: Expected 1-D numpy arrays.\nFound (nxdims, nydims, nzdims) = (%d, %d, %d) instead",
                 nxdims, nydims, nzdims);
        countpairs_mocks_error_out(module, msg);
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
        countpairs_mocks_error_out(module, msg);
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
          snprintf(msg, 1024, "TypeError: Expected *ALL* 3 floating point arrays to be the same type (allowed types = %d or %d). Instead found type-nums (%d, %d, %d) "
                   "with type-names = (%s, %s, %s)\n",
                   NPY_FLOAT, NPY_DOUBLE, x_type, y_type, z_type, x_descr->typeobj->tp_name, y_descr->typeobj->tp_name, z_descr->typeobj->tp_name);
        }
        Py_XDECREF(x_descr);Py_XDECREF(y_descr);Py_XDECREF(z_descr);
        countpairs_mocks_error_out(module, msg);
        return -1;
    }
    
    /* Check if the number of elements in the 3 Python arrays are identical */
    const int64_t nx1 = (int64_t)PyArray_SIZE(x1_obj);
    const int64_t ny1 = (int64_t)PyArray_SIZE(y1_obj);
    const int64_t nz1 = (int64_t)PyArray_SIZE(z1_obj);

    if(nx1 != ny1 || ny1 != nz1) {
      snprintf(msg, 1024, "ERROR: Expected arrays to have the same number of elements in all 3-dimensions.\nFound (nx, ny, nz) = (%"PRId64", %"PRId64", %"PRId64") instead",
               nx1, ny1, nz1);
      countpairs_mocks_error_out(module, msg);
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


static int64_t check_dims_and_datatype_ra_dec(PyObject *module, PyArrayObject *x1_obj, PyArrayObject *y1_obj, size_t *element_size)
{
    char msg[1024];

    /* All the arrays should be 1-D*/
    const int nxdims = PyArray_NDIM(x1_obj);
    const int nydims = PyArray_NDIM(y1_obj);

    if(nxdims != 1 || nydims != 1) {
        snprintf(msg, 1024, "ERROR: Expected 1-D numpy arrays.\nFound (nxdims, nydims) = (%d, %d) instead",
                 nxdims, nydims);
        countpairs_mocks_error_out(module, msg);
        return -1;
    }

    /* All the arrays should be floating point (only float32 and float64 are allowed) */
    const int x_type = PyArray_TYPE(x1_obj);
    const int y_type = PyArray_TYPE(y1_obj);
    if( ! ((x_type == NPY_FLOAT || x_type == NPY_DOUBLE) &&
           (y_type == NPY_FLOAT || y_type == NPY_DOUBLE))
        ) {
        PyArray_Descr *x_descr = PyArray_DescrFromType(x_type);
        PyArray_Descr *y_descr = PyArray_DescrFromType(y_type);
        if(x_descr == NULL || y_descr == NULL) {
            /* Generating the dtype descriptor failed somehow. At least provide some information */
            snprintf(msg, 1024, "TypeError: Expected 2 floating point arrays (allowed types = %d or %d). Instead found type-nums (%d, %d)\n",
                     NPY_FLOAT, NPY_DOUBLE, x_type, y_type);
        } else {
            snprintf(msg, 1024, "TypeError: Expected 2 floating point arrays (allowed types = %d or %d). Instead found type-nums (%d, %d) "
                     "with type-names = (%s, %s)\n",
                     NPY_FLOAT, NPY_DOUBLE, x_type, y_type, x_descr->typeobj->tp_name, y_descr->typeobj->tp_name);
        }
        Py_XDECREF(x_descr);Py_XDECREF(y_descr);
        countpairs_mocks_error_out(module, msg);
        return -1;
    }

    if( x_type != y_type) {
        PyArray_Descr *x_descr = PyArray_DescrFromType(x_type);
        PyArray_Descr *y_descr = PyArray_DescrFromType(y_type);
        if(x_descr == NULL || y_descr == NULL) {
          /* Generating the dtype descriptor failed somehow. At least provide some information */
            snprintf(msg, 1024, "TypeError: Expected *BOTH* (RA, DEC) floating point arrays to be the same type (allowed types = %d or %d). Instead found type-nums (%d, %d)\n",
                   NPY_FLOAT, NPY_DOUBLE, x_type, y_type);
        } else {
          snprintf(msg, 1024, "TypeError: Expected *BOTH* (RA, DEC) floating point arrays to be the same type (allowed types = %d or %d). Instead found type-nums (%d, %d) "
                   "with type-names = (%s, %s)\n",
                   NPY_FLOAT, NPY_DOUBLE, x_type, y_type, x_descr->typeobj->tp_name, y_descr->typeobj->tp_name);
        }
        Py_XDECREF(x_descr);Py_XDECREF(y_descr);
        countpairs_mocks_error_out(module, msg);
        return -1;
    }
    
    /* Check if the number of elements in the 3 Python arrays are identical */
    const int64_t nx1 = (int64_t)PyArray_SIZE(x1_obj);
    const int64_t ny1 = (int64_t)PyArray_SIZE(y1_obj);

    if(nx1 != ny1) {
        snprintf(msg, 1024, "ERROR: Expected *BOTH* (RA, DEC) arrays to have the same number of elements in all 2-dimensions.\nFound (nra, ndec) = "
                 "(%"PRId64", %"PRId64") instead", nx1, ny1);
      countpairs_mocks_error_out(module, msg);
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


static PyObject *countpairs_countpairs_rp_pi_mocks(PyObject *self, PyObject *args, PyObject *kwargs)
{
    //Error-handling is global in python2 -> stored in struct module_state _struct declared at the top of this file
#if PY_MAJOR_VERSION < 3
    (void) self;
    PyObject *module = NULL;//should not be used -> setting to NULL so any attempts to dereference will result in a crash. 
#else
    //In python3, self is simply the module object that was returned earlier by init
    PyObject *module = self;
#endif    

    //x1->ra (phi), y1-> declination (theta1), z1->cz (cz1)
    //x2->ra (ph2), y2-> declination (theta2), z2->cz (cz2)
    PyArrayObject *x1_obj=NULL, *y1_obj=NULL, *z1_obj=NULL, *x2_obj=NULL,*y2_obj=NULL,*z2_obj=NULL;

    struct config_options options;
    memset(&options, 0, sizeof(struct config_options));
    options.verbose = 0;
    options.instruction_set = AVX;
    options.periodic = 0;
    options.fast_divide=0;
    
    int autocorr=0;
    int nthreads=4;
    int cosmology=1;
    double pimax;
    char *binfile;

    static char *kwlist[] = {
        "autocorr",
        "cosmology",
        "nthreads",
        "pimax",
        "binfile",
        "X1",
        "Y1",
        "Z1",
        "X2",
        "Y2",
        "Z2",
        "verbose", /* keyword verbose -> print extra info at runtime + progressbar */
        "output_rpavg",
        "fast_divide",
        "isa",/* instruction set to use of type enum isa; valid values are AVX, SSE, FALLBACK */
        NULL
    };


    if ( ! PyArg_ParseTupleAndKeywords(args, kwargs, "iiidsO!O!O!|O!O!O!iii", kwlist,
                                       &autocorr,&cosmology,&nthreads,&pimax,&binfile,
                                       &PyArray_Type,&x1_obj,
                                       &PyArray_Type,&y1_obj,
                                       &PyArray_Type,&z1_obj,
                                       &PyArray_Type,&x2_obj,//optional parameters -> if autocorr == 1, not checked; required if autocorr=0
                                       &PyArray_Type,&y2_obj,
                                       &PyArray_Type,&z2_obj,
                                       &(options.verbose),
                                       &(options.need_avg_sep),
                                       &(options.fast_divide),
                                       &(options.instruction_set))

         ) {
        Py_RETURN_NONE;
    }

    
    /* We have numpy arrays and all the required inputs*/
    /* How many data points are there? And are they all of floating point type */
    size_t element_size;
    const int64_t ND1 = check_dims_and_datatype(module, x1_obj, y1_obj, z1_obj, &element_size);
    if(ND1 == -1) {
        //Error has already been set -> simply return 
        Py_RETURN_NONE;
    }

    int64_t ND2 = ND1;
    if(autocorr == 0) {
        char msg[1024];
        if(x2_obj == NULL || y2_obj == NULL || z2_obj == NULL) {
            snprintf(msg, 1024, "ValueError: In %s: If autocorr is 0, need to pass the second set of positions (X2=numpy array, Y2=numpy array, Z2=numpy array).\n",
                     __FUNCTION__);
            countpairs_mocks_error_out(module, msg);
            Py_RETURN_NONE;
        }
        size_t element_size2;
        ND2 = check_dims_and_datatype(module, x2_obj, y2_obj, z2_obj, &element_size2);
        if(ND2 == -1) {
            //Error has already been set -> simply return 
            Py_RETURN_NONE;
        }
        if(element_size != element_size2) {
            snprintf(msg, 1024, "TypeError: In %s: The two arrays must have the same data-type. First array is of type %s while second array is of type %s\n",
                     __FUNCTION__, element_size == 4 ? "floats":"doubles", element_size2 == 4 ? "floats":"doubles");
            countpairs_mocks_error_out(module, msg);
            Py_RETURN_NONE;
        }
    }


    /* Interpret the input objects as numpy arrays. */
    const int requirements = NPY_ARRAY_IN_ARRAY;
    PyObject *x1_array = NULL, *y1_array = NULL, *z1_array = NULL;
    PyObject *x2_array = NULL, *y2_array = NULL, *z2_array = NULL;
    x1_array = PyArray_FromArray(x1_obj, NOTYPE_DESCR, requirements);
    y1_array = PyArray_FromArray(y1_obj, NOTYPE_DESCR, requirements);
    z1_array = PyArray_FromArray(z1_obj, NOTYPE_DESCR, requirements);

    if(autocorr == 0) {
        x2_array = PyArray_FromArray(x2_obj, NOTYPE_DESCR, requirements);
        y2_array = PyArray_FromArray(y2_obj, NOTYPE_DESCR, requirements);
        z2_array = PyArray_FromArray(z2_obj, NOTYPE_DESCR, requirements);
    }

    if (x1_array == NULL || y1_array == NULL || z1_array == NULL ||
        (autocorr == 0 && (x2_array == NULL || y2_array == NULL || z2_array == NULL))) {
        Py_XDECREF(x1_array);
        Py_XDECREF(y1_array);
        Py_XDECREF(z1_array);

        Py_XDECREF(x2_array);
        Py_XDECREF(y2_array);
        Py_XDECREF(z2_array);
        char msg[1024];
        snprintf(msg, 1024, "TypeError: In %s: Could not convert input to arrays of allowed floating point types (doubles or floats). Are you passing numpy arrays?",
                 __FUNCTION__);
        countpairs_mocks_error_out(module, msg);
        Py_RETURN_NONE;
    }

    /* Get pointers to the data as C-types. */
    void *phiD1=NULL, *thetaD1=NULL, *czD1=NULL;
    void *phiD2=NULL, *thetaD2=NULL, *czD2=NULL;

    phiD1   = PyArray_DATA((PyArrayObject *)x1_array);
    thetaD1 = PyArray_DATA((PyArrayObject *)y1_array);
    czD1    = PyArray_DATA((PyArrayObject *)z1_array);

    if(autocorr == 0) {
        phiD2   = PyArray_DATA((PyArrayObject *) x2_array);
        thetaD2 = PyArray_DATA((PyArrayObject *) y2_array);
        czD2    = PyArray_DATA((PyArrayObject *) z2_array);
    }
    options.float_type = element_size;
    
    NPY_BEGIN_THREADS_DEF;
    NPY_BEGIN_THREADS;

    results_countpairs_mocks results;
    int status = countpairs_mocks(ND1,phiD1,thetaD1,czD1,
                                  ND2,phiD2,thetaD2,czD2,
                                  nthreads,
                                  autocorr,
                                  binfile,
                                  pimax,
                                  cosmology,
                                  &results,
                                  &options);

    NPY_END_THREADS;
    
    /* Clean up. */
    Py_DECREF(x1_array);Py_DECREF(y1_array);Py_DECREF(z1_array);//x1 should absolutely not be NULL
    Py_XDECREF(x2_array);Py_XDECREF(y2_array);Py_XDECREF(z2_array);//x2 might be NULL depending on value of autocorr

    if(status != EXIT_SUCCESS) {
        Py_RETURN_NONE;
    }
    
    
#if 0
    /* Output pairs*/
    for(int i=1;i<results.nbin;i++) {
        const double logrp = LOG10(results.rupp[i]);
        for(int j=0;j<npibin;j++) {
            const int index = i*(npibin+1) + j;
            fprintf(stdout,"%10"PRIu64" %20.8lf %20.8lf  %20.8lf \n",results.npairs[index],results.rpavg[index],logrp,(j+1)*dpi);
        }
    }
#endif


    /* Build the output list */
    PyObject *ret = PyList_New(0);//create an empty list
    double rlow=results.rupp[0];
    const double dpi = pimax/(double)results.npibin ;

    for(int i=1;i<results.nbin;i++) {
        for(int j=0;j<results.npibin;j++) {
            const int bin_index = i*(results.npibin + 1) + j;
            PyObject *item = NULL;
            const double rpavg = results.rpavg[bin_index];
            item = Py_BuildValue("(ddddk)", rlow,results.rupp[i],rpavg,(j+1)*dpi,results.npairs[bin_index]);
            PyList_Append(ret, item);
            Py_XDECREF(item);
        }
        rlow=results.rupp[i];
    }
    free_results_mocks(&results);
    return ret;
}

static PyObject *countpairs_countpairs_theta_mocks(PyObject *self, PyObject *args, PyObject *kwargs)
{
    //Error-handling is global in python2 -> stored in struct module_state _struct declared at the top of this file
#if PY_MAJOR_VERSION < 3
    (void) self;
    PyObject *module = NULL;//should not be used -> setting to NULL so any attempts to dereference will result in a crash. 
#else
    //In python3, self is simply the module object that was returned earlier by init
    PyObject *module = self;
#endif    

    PyArrayObject *x1_obj=NULL, *y1_obj=NULL, *x2_obj=NULL,*y2_obj=NULL;
    int nthreads=1;
    char *binfile;
    int autocorr=0;
    struct config_options options;
    memset(&options, 0, sizeof(struct config_options));
    options.verbose=0;
    options.instruction_set=AVX;
    options.link_in_dec=1;
    options.link_in_ra=1;
    options.fast_acos=0;
    static char *kwlist[] = {
        "autocorr",
        "nthreads",
        "binfile",
        "X1",
        "Y1",
        "X2",
        "Y2",
        "link_in_dec",
        "link_in_ra",
        "verbose", /* keyword verbose -> print extra info at runtime + progressbar */
        "output_thetaavg",
        "fast_acos",
        "isa",/* instruction set to use of type enum isa; valid values are AVX, SSE, FALLBACK */
        NULL
    };


    if ( ! PyArg_ParseTupleAndKeywords(args, kwargs, "iisO!O!|O!O!iiiiii", kwlist,
                                       &autocorr,&nthreads,&binfile,
                                       &PyArray_Type,&x1_obj,
                                       &PyArray_Type,&y1_obj,
                                       &PyArray_Type,&x2_obj,//optional parameters -> if autocorr == 1, not checked; required if autocorr=0
                                       &PyArray_Type,&y2_obj,
                                       &(options.link_in_ra),
                                       &(options.link_in_dec),
                                       &(options.verbose),
                                       &(options.need_avg_sep),
                                       &(options.fast_acos),
                                       &(options.instruction_set))

         ) {
        Py_RETURN_NONE;
    }
    options.autocorr=autocorr;
    options.periodic=0;//doesn't matter but noting intent by setting it to 0


    size_t element_size;
    /* We have numpy arrays and all the required inputs*/
    /* How many data points are there? And are they all of floating point type */
    const int64_t ND1 = check_dims_and_datatype_ra_dec(module, x1_obj, y1_obj, &element_size);
    if(ND1 == -1) {
        //Error has already been set -> simply return 
        Py_RETURN_NONE;
    }

    int64_t ND2 = ND1;
    if(autocorr==0) {
        char msg[1024];
        if(x2_obj == NULL || y2_obj == NULL) {
            snprintf(msg, 1024, "ValueError: In %s: If autocorr is 0, need to pass the second set of positions (RA2=numpy array, DEC2=numpy array).\n",
                     __FUNCTION__);
            countpairs_mocks_error_out(module, msg);
            Py_RETURN_NONE;
        }
        size_t element_size2;
        ND2 = check_dims_and_datatype_ra_dec(module, x2_obj, y2_obj,&element_size2);
        if(ND2 == -1) {
            //Error has already been set -> simply return 
            Py_RETURN_NONE;
        }
        if(element_size != element_size2) {
            snprintf(msg, 1024, "TypeError: In %s: The two arrays must have the same data-type. First array is of type %s while second array is of type %s\n",
                     __FUNCTION__, element_size == 4 ? "floats":"doubles", element_size2 == 4 ? "floats":"doubles");
            countpairs_mocks_error_out(module, msg);
            Py_RETURN_NONE;
        }
    }

    /* Interpret the input objects as numpy arrays. */
    const int requirements = NPY_ARRAY_IN_ARRAY;
    PyObject *x1_array = NULL, *y1_array = NULL;
    PyObject *x2_array = NULL, *y2_array = NULL;
    x1_array = PyArray_FromArray(x1_obj, NOTYPE_DESCR, requirements);
    y1_array = PyArray_FromArray(y1_obj, NOTYPE_DESCR, requirements);

    if(autocorr == 0) {
        x2_array = PyArray_FromArray(x2_obj, NOTYPE_DESCR, requirements);
        y2_array = PyArray_FromArray(y2_obj, NOTYPE_DESCR, requirements);
    }

    if (x1_array == NULL || y1_array == NULL || 
        (autocorr == 0 && (x2_array == NULL || y2_array == NULL))) {
        Py_XDECREF(x1_array);
        Py_XDECREF(y1_array);

        Py_XDECREF(x2_array);
        Py_XDECREF(y2_array);
        char msg[1024];
        snprintf(msg, 1024, "TypeError: In %s: Could not convert input to arrays of allowed floating point types (doubles or floats). Are you passing numpy arrays?",
                 __FUNCTION__);
        countpairs_mocks_error_out(module, msg);
        Py_RETURN_NONE;
    }


    /* Get pointers to the data as C-types. */
    void *phiD1 = NULL, *thetaD1 = NULL;
    void *phiD2 = NULL, *thetaD2 = NULL;
    phiD1   = PyArray_DATA((PyArrayObject *) x1_array); 
    thetaD1 = PyArray_DATA((PyArrayObject *) y1_array);

    if(autocorr == 0) {
        phiD2   = PyArray_DATA((PyArrayObject *) x2_array);
        thetaD2 = PyArray_DATA((PyArrayObject *) y2_array);
    }

    NPY_BEGIN_THREADS_DEF;
    NPY_BEGIN_THREADS;

    results_countpairs_theta results;
    options.float_type = element_size;
    int status = countpairs_theta_mocks(ND1,phiD1,thetaD1,
                                        ND2,phiD2,thetaD2,
                                        nthreads,
                                        autocorr,
                                        binfile,
                                        &results,
                                        &options);
    NPY_END_THREADS;
    
    /* Clean up. */
    Py_DECREF(x1_array);Py_DECREF(y1_array);//x1/y1 (representing ra1,dec1) should not be NULL
    Py_XDECREF(x2_array);Py_XDECREF(y2_array);//x2/y2 may be NULL (in case of autocorr)

    if(status != EXIT_SUCCESS) {
        Py_RETURN_NONE;
    }
    
#if 0
    /*---Output-Pairs-------------------------------------*/
    double theta_low = results.theta_upp[0];
    for(int i=1;i<results.nbin;i++) {
        fprintf(stdout,"%10"PRIu64" %20.8lf %20.8lf %20.8lf \n",results.npairs[i],results.theta_avg[i],theta_low,results.theta_upp[i]);
        theta_low=results.theta_upp[i];
    }
#endif

    /* Build the output list */
    PyObject *ret = PyList_New(0);
    double rlow=results.theta_upp[0];
    for(int i=1;i<results.nbin;i++) {
        PyObject *item = NULL;
        const double theta_avg = results.theta_avg[i];
        item = Py_BuildValue("(dddk)", rlow,results.theta_upp[i],theta_avg,results.npairs[i]);
        PyList_Append(ret, item);
        Py_XDECREF(item);
        rlow=results.theta_upp[i];
    }
    free_results_countpairs_theta(&results);
    return ret;
}



static PyObject *countpairs_countspheres_vpf_mocks(PyObject *self, PyObject *args, PyObject *kwargs)
{
    //Error-handling is global in python2 -> stored in struct module_state _struct declared at the top of this file
#if PY_MAJOR_VERSION < 3
    (void) self;
    PyObject *module = NULL;//should not be used -> setting to NULL so any attempts to dereference will result in a crash. 
#else
    //In python3, self is simply the module object that was returned earlier by init
    PyObject *module = self;
#endif    

    //x1->ra (phi), y1-> declination (theta1), z1->cz (cz1)
    //x2->ra (ph2), y2-> declination (theta2), z2->cz (cz2)
    PyArrayObject *x1_obj=NULL, *y1_obj=NULL, *z1_obj=NULL, *x2_obj=NULL,*y2_obj=NULL,*z2_obj=NULL;
    int cosmology=1;
    double rmax;
    int nbin,num_spheres, num_pN;
    char *centers_file;
    int threshold_neighbors;
    struct config_options options;
    memset(&options, 0, sizeof(struct config_options));
    options.verbose=0;
    options.instruction_set=AVX;
    static char *kwlist[] = {
        "rmax",
        "nbins",
        "numSpheres",
        "numpN",
        "threshold",
        "centers_file",
        "cosmology",
        "X1",
        "Y1",
        "Z1",
        "X2",
        "Y2",
        "Z2",
        "verbose", /* keyword verbose -> print extra info at runtime + progressbar */
        "isa",/* instruction set to use of type enum isa; valid values are AVX, SSE, FALLBACK */
        NULL
    };


    if ( ! PyArg_ParseTupleAndKeywords(args, kwargs, "diiiisiO!O!O!O!O!O!|ii", kwlist,
                                       &rmax,&nbin,&num_spheres,&num_pN,&threshold_neighbors,&centers_file,&cosmology,
                                       &PyArray_Type,&x1_obj,
                                       &PyArray_Type,&y1_obj,
                                       &PyArray_Type,&z1_obj,
                                       &PyArray_Type,&x2_obj,//optional parameters -> if autocorr == 1, not checked; required if autocorr=0
                                       &PyArray_Type,&y2_obj,
                                       &PyArray_Type,&z2_obj,
                                       &(options.verbose),
                                       &(options.instruction_set))

         ) {
        Py_RETURN_NONE;
    }
    
    size_t element_size;
    /* We have numpy arrays and all the required inputs*/
    /* How many data points are there? And are they all of floating point type */
    const int64_t ND1 = check_dims_and_datatype(module, x1_obj, y1_obj, z1_obj, &element_size);
    if(ND1 == -1) {
        //Error has already been set -> simply return 
        Py_RETURN_NONE;
    }

    size_t element_size2;
    const int64_t ND2 = check_dims_and_datatype(module, x2_obj, y2_obj, z2_obj, &element_size2);
    if(ND2 == -1) {
        //Error has already been set -> simply return 
        Py_RETURN_NONE;
    }

    if(element_size != element_size2) {
        char msg[1024];
        snprintf(msg, 1024, "TypeError: In %s: The two arrays must have the same data-type. First array is of type %s while second array is of type %s\n",
                 __FUNCTION__, element_size == 4 ? "floats":"doubles", element_size2 == 4 ? "floats":"doubles");
        countpairs_mocks_error_out(module, msg);
        Py_RETURN_NONE;
    }

    
    /* Interpret the input objects as numpy arrays. */
    const int requirements = NPY_ARRAY_IN_ARRAY;
    PyObject *x1_array = NULL, *y1_array = NULL, *z1_array = NULL;
    PyObject *x2_array = NULL, *y2_array = NULL, *z2_array = NULL;
    x1_array = PyArray_FromArray(x1_obj, NOTYPE_DESCR, requirements);
    y1_array = PyArray_FromArray(y1_obj, NOTYPE_DESCR, requirements);
    z1_array = PyArray_FromArray(z1_obj, NOTYPE_DESCR, requirements);
    x2_array = PyArray_FromArray(x2_obj, NOTYPE_DESCR, requirements);
    y2_array = PyArray_FromArray(y2_obj, NOTYPE_DESCR, requirements);
    z2_array = PyArray_FromArray(z2_obj, NOTYPE_DESCR, requirements);


    if (x1_array == NULL || y1_array == NULL || z1_array == NULL ||
        x2_array == NULL || y2_array == NULL || z2_array == NULL) {
        Py_XDECREF(x1_array);
        Py_XDECREF(y1_array);
        Py_XDECREF(z1_array);

        Py_XDECREF(x2_array);
        Py_XDECREF(y2_array);
        Py_XDECREF(z2_array);
        char msg[1024];
        snprintf(msg, 1024, "TypeError: In %s: Could not convert input to arrays of allowed floating point types (doubles or floats). Are you passing numpy arrays?",
                 __FUNCTION__);
        countpairs_mocks_error_out(module, msg);
        Py_RETURN_NONE;
    }

    /* Get pointers to the data as C-types. */
    void *phiD1=NULL, *thetaD1=NULL,*czD1=NULL;
    void *phiD2=NULL, *thetaD2=NULL,*czD2=NULL;
    
    phiD1   = PyArray_DATA((PyArrayObject *) x1_array);
    thetaD1 = PyArray_DATA((PyArrayObject *) y1_array);
    czD1    = PyArray_DATA((PyArrayObject *) z1_array);

    phiD2   = PyArray_DATA((PyArrayObject *) x2_array);
    thetaD2 = PyArray_DATA((PyArrayObject *) y2_array);
    czD2    = PyArray_DATA((PyArrayObject *) z2_array);

    NPY_BEGIN_THREADS_DEF;
    NPY_BEGIN_THREADS;
    
    results_countspheres_mocks results;
    options.float_type = element_size;
    int status = countspheres_mocks(ND1, phiD1,thetaD1, czD1,
                                    ND2, phiD2,thetaD2, czD2,
                                    threshold_neighbors,
                                    rmax, nbin, num_spheres,
                                    num_pN,
                                    centers_file,
                                    cosmology,
                                    &results,
                                    &options);
    
    NPY_END_THREADS;

    /* Clean up. */
    Py_DECREF(x1_array);Py_DECREF(y1_array);Py_DECREF(z1_array);
    Py_DECREF(x2_array);Py_DECREF(y2_array);Py_DECREF(z2_array);

    if(status != EXIT_SUCCESS) {
        Py_RETURN_NONE;
    }
    
#if 0
    // Output the results
    const double rstep = rmax/(double)nbin ;
    for(int ibin=0;ibin<results.nbin;ibin++) {
        const double r=(ibin+1)*rstep;
        fprintf(stdout,"%10.2"REAL_FORMAT" ", r);
        for(int i=0;i<num_pN;i++) {
            fprintf(stdout," %10.4e", (results.pN)[ibin][i]);
        }
        fprintf(stdout,"\n");
    }
#endif

    /* Build the output list (of lists, since num_pN is determined at runtime) */
    PyObject *ret = PyList_New(0);
    const double rstep = rmax/(double)nbin ;
    for(int ibin=0;ibin<results.nbin;ibin++) {
        const double r=(ibin+1)*rstep;
        PyObject *item = PyList_New(0);
        PyObject *this_val = Py_BuildValue("d",r);
        PyList_Append(item, this_val);
        Py_XDECREF(this_val);
        for(int i=0;i<num_pN;i++) {
            this_val = Py_BuildValue("d",(results.pN)[ibin][i]);
            PyList_Append(item, this_val);
            Py_XDECREF(this_val);
        }
        PyList_Append(ret, item);
        Py_XDECREF(item);
    }
    free_results_countspheres_mocks(&results);

    return ret;
}
