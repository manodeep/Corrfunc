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
#include "countpairs_s_mu_mocks.h"
#include "countpairs_theta_mocks.h"

//for the vpf
#include "countspheres_mocks.h"

//for the instruction set detection
#include "cpu_features.h"

//for unicode characters
#include "macros.h"


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



static int highest_isa_mocks;

//Docstrings for the methods
static char module_docstring[] =    "Python extensions for calculating clustering statistics on MOCKS (spherical geometry).\n"
    "\n"
    "countpairs_rp_pi_mocks: Calculate the 2-D DD("RP_CHAR","PI_CHAR") auto/cross-correlation function given two sets of ra/dec/cz and ra/dec/cz arrays.\n"
    "countpairs_theta_mocks: Calculate DD(theta) auto/cross-correlation function given two sets of ra/dec/cz and ra/dec/cz arrays.\n"
    "countspheres_vpf_mocks: Calculate the counts-in-spheres given one set of ra/dec/cz.\n"
    "\n\n"
    "See `Corrfunc/call_correlation_functions_mocks.py` for example calls to each function.\n";

/* function proto-type*/
static PyObject *countpairs_countpairs_rp_pi_mocks(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *countpairs_countpairs_s_mu_mocks(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *countpairs_countpairs_theta_mocks(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *countpairs_countspheres_vpf_mocks(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *countpairs_mocks_error_out(PyObject *module, const char *msg);

static PyMethodDef module_methods[] = {
    {"countpairs_rp_pi_mocks"       ,(PyCFunction)(void(*)(void)) countpairs_countpairs_rp_pi_mocks ,METH_VARARGS | METH_KEYWORDS,
     "countpairs_rp_pi_mocks(autocorr, cosmology, nthreads, pimax, binfile,\n"
     "                       RA1, DEC1, CZ1, weights1=None, weight_type=None,\n"
     "                       RA2=None, DEC2=None, CZ2=None, weights2=None,\n"
     "                       is_comoving_dist=False,\n"
     "                       verbose=False, output_rpavg=False,\n"
     "                       fast_divide_and_NR_steps=0, xbin_refine_factor=2, \n"
     "                       ybin_refine_factor=2, zbin_refine_factor=1, \n"
     "                       max_cells_per_dim=100, copy_particles=True,\n"
     "                       enable_min_sep_opt=True, c_api_timer=False, isa=-1)\n\n"

     "Calculate the 2-D pair-counts, "XI_CHAR"("RP_CHAR", "PI_CHAR"), auto/cross-correlation function given two\n"
     "sets of RA1/DEC1/CZ1 and RA2/DEC2/CZ2 arrays. This module is suitable for mock catalogs that have been\n"
     "created by carving out a survey footprint from simulated data. The module can also be used for actual\n"
     "observed galaxies, but you probably want to attach weights to the points to account for completeness etc.\n"
     "Default bins in "PI_CHAR" are set to 1.0 Mpc/h.\n\n"

     "Note, that this module only returns pair counts and not the actual correlation function\n"
     ""XI_CHAR"("RP_CHAR", "PI_CHAR"). See the `~Corrfunc.utils.convert_3d_counts_to_cf.py and \n"
     "`~Corrfunc.utils.convert_rp_pi_counts_to_wp.py for computing "XI_CHAR"("RP_CHAR","PI_CHAR") \n"
     "and wp("RP_CHAR") from DD("RP_CHAR", "PI_CHAR").\n"
     UNICODE_WARNING"\n\n"

     "Parameters\n"
     "----------\n"
     "Every parameter can be passed as a keyword of the corresponding name.\n\n"

     "autocorr: boolean\n"
     "    Flag for auto/cross-correlation. If autocorr is not 0, the RA2/DEC2/CZ2 arrays\n"
     "    are not used (but must still be passed, as RA1/DEC1/CZ1).\n\n"

     "cosmology: integer\n"
     "    Integer to select cosmology. Pre-set values for (1,2) \n"
     "    1 -> LasDamas cosmology. Om=0.25,  Ol=0.75  (other values are not used)\n"
     "    2 -> Planck   cosmology. Om=0.302, Ol=0.698 \n"
     "    To setup a new cosmology, add an entry to the function, `init_cosmology` in \n"
     "    `ROOT/utils/cosmology_params.c` and recompile the package.\n\n"

     "nthreads: integer\n"
     "    The number of OpenMP threads to use. Has no effect if OpenMP was not used\n"
     "    during library compilation. \n\n"

     ""PIMAX_CHAR": double (Mpc/h)\n"
     "    The max. integration distance along the "PI_CHAR" direction in Mpc/h. Typical\n"
     "    values are in the range 40-100 Mpc/h.\n\n"

     "binfile: filename\n"
     "    Filename containing the radial bins for the correlation function. The file\n"
     "    is expected to contain white-space separated ``rpmin  rpmax`` with the bin\n"
     "    edges.  Units must be Mpc/h (see the ``bins`` file in the tests directory\n"
     "    for a sample). For usual logarithmic bins, ``logbins``in the root directory\n"
     "    of this package will create a compatible ``binfile``.\n\n"

     "RA1: array-like, float/double (default double)\n"
     "    The right-ascension of the galaxy, in the range [0, 360]. If there are\n"
     "    negative RA's in the supplied array (input RA in the range [-180, 180]),\n"
     "    then the code will shift the entire array by 180 to put RA's in the\n"
     "    [0, 360] range.\n\n"

     "DEC1: array-like, float/double (default double)\n"
     "    The declination of the galaxy, in the range [-90, 90]. If there are\n"
     "    declinations > 90 in the supplied array (input dec in the range [0, 180]),\n"
     "    then the code will shift the entire array by -90 to put declinations in\n"
     "    the [-90, 90] range. If the code finds declinations more than 180, then\n"
     "    it assumes RA and DEC have been swapped and aborts with that message.\n\n"

     "CZ1: array-like, float/double (default double)\n"
     "    The redshift multiplied by speed of light for the galaxies. The code will\n"
     "    checks that cz has been supplied by comparing with a threshold (currently\n"
     "    set to 10, defined in function check_ra_dec_cz in file\n"
     "    `DDrppi_mocks/countpairs_rp_pi_mocks_impl.c.src`) and multiplies by the speed of light if\n"
     "    max z is less than that threshold. If you really want to change the speed\n"
     "    of light, then edit the macro in `ROOT/utils/set_cosmo_dist.h`.\n\n"

     "weights1 : array-like, real (float/double), shape (n_weights_per_particle,n_particles), optional\n"
     "    Weights for computing a weighted pair count.\n\n"

     "weight_type : str, optional\n"
     "    The type of pair weighting to apply.\n"
     "    Options: \"pair_product\", None\n"
     "    Default: None.\n\n"

     "RA2/DEC2/CZ2: float/double (default double)\n"
     "    Same as for RA1/DEC1/CZ1\n\n"

     "weights2\n : array-like, real (float/double), shape (n_weights_per_particle,n_particles), optional\n"
     "    Weights for computing a weighted pair count.\n\n"

     "is_comoving_dist: boolean (default false)\n"
     "    Boolean flag to indicate that ``cz`` values have already been\n"
     "    converted into co-moving distances. This flag allows arbitrary\n"
     "    cosmologies to be used in ``Corrfunc``.\n\n"

     "verbose : boolean (default false)\n"
     "    Boolean flag to control output of informational messages\n\n"

     "output_rpavg : boolean (default false)\n"
     "    Boolean flag to output the average ``rp`` for each bin. Code will\n"
     "    run slightly slower if you set this flag. Also, note, if you are calculating\n"
     "    in single-precision, ``rpavg`` will suffer from numerical loss of\n"
     "    precision and can not be trusted. If you need accurate ``rpavg``\n"
     "    values, then pass in double precision arrays for the particle positions.\n\n"

     "fast_divide_and_NR_steps: integer (default 0)\n"
     "    Replaces the division in ``AVX512F`` and ``AVX`` implementation with an\n"
     "    approximate reciprocal, followed by ``fast_divide_and_NR_steps`` "
     "    Newton-Raphson step. Can improve \n"
     "    runtime by ~15-20%. Value of 0 keeps the standard division.\n\n"

     "(xyz)bin_refine_factor: integer (default (2,2,1) typical values in [1-3]) \n"
     "    Controls the refinement on the cell sizes. Can have up to a 20% impact \n"
     "    on runtime. \n\n"

     "max_cells_per_dim: integer (default 100, typical values in [50-300]) \n"
     "    Controls the maximum number of cells per dimension. Total number of cells \n"
     "    can be up to (max_cells_per_dim)^3. Only increase if ``rmax`` is too small \n"
     "    relative to the boxsize (and increasing helps the runtime). \n\n"

     "copy_particles: boolean (default True)\n"
     "    Boolean flag to make a copy of the particle positions\n"
     "    If set to False, the particles will be re-ordered in-place\n\n"

     "enable_min_sep_opt: boolean (default true)\n"
     "    Boolean flag to allow optimizations based on min. separation between pairs \n"
     "    of cells. Here to allow for comparison studies.\n\n"

     "c_api_timer : boolean (default false)\n"
     "    Boolean flag to measure actual time spent in the C libraries. Here\n"
     "    to allow for benchmarking and scaling studies.\n\n"

     "isa : integer (default -1)\n"
     "    Controls the runtime dispatch for the instruction set to use. Possible\n"
     "    options are: [-1, AVX512F, AVX, SSE42, FALLBACK]\n\n"
     "    Setting isa to -1 will pick the fastest available instruction\n"
     "    set on the current computer. However, if you set ``isa`` to, say,\n"
     "    ``AVX`` and ``AVX`` is not available on the computer, then the code will\n"
     "    revert to using ``FALLBACK`` (even though ``SSE42`` might be available).\n"
     "\n"
     "    Unless you are benchmarking the different instruction sets, you should\n"
     "    always leave ``isa`` to the default value. And if you *are* benchmarking,\n"
     "    then the integer values correspond to the ``enum`` for the instruction set\n"
     "    defined in ``utils/defs.h``.\n\n"

     "Returns\n"
     "--------\n\n"

     "A tuple (results, time) \n\n"

     "results : Python list \n"
     "    Contains [rpmin, rpmax, rpavg, "PI_CHAR", npairs, weightavg] \n"
     "    for each "PI_CHAR"-bin (up to "PIMAX_CHAR") for each radial bin specified in\n"
     "    the ``binfile``. For instance, for a ``"PIMAX_CHAR"`` of 40.0 Mpc/h, each radial\n"
     "    bin will be split into 40 "PI_CHAR" bins (default "PI_CHAR" bin is 1.0). Thus, the\n"
     "    total number of items in the list is {(int) ``"PIMAX_CHAR"`` * number of rp bins}.\n"
     "    If ``output_rpavg`` is not set then ``rpavg`` will be set to 0.0 for all bins; similarly for ``weight_avg``. \n"
     "    "PI_CHAR" for each bin is the upper limit of the "PI_CHAR" values that were \n"
     "    considered in that ("RP_CHAR", "PI_CHAR") bin. ``npairs``contains the number of pairs\n"
     "    in that bin and can be used to compute the actual "XI_CHAR"("RP_CHAR", "PI_CHAR") by\n"
     "    combining with RR counts. combining  DD, (DR) and RR counts. \n"
     "    See ``~Corrfunc.utils.convert_3d_counts_to_cf``\n\n"

     "time : double\n"
     "    if ``c_api_timer`` is set, then the return value contains the time spent\n"
     "    in the API; otherwise time is set to 0.0\n\n"

     "Example\n"
     "-------\n"
     ">>> import numpy as np\n"
     ">>> from Corrfunc._countpairs_mocks import countpairs_rp_pi_mocks\n"
     ">>> ra,dec,cz = np.genfromtxt('../mocks/tests/data/Mr19_mock_northonly.rdcz.dat',dtype=np.float,unpack=True)\n"
     ">>> cosmology=1\n"
     ">>> autocorr=1\n"
     ">>> nthreads=4\n"
     ">>> binfile='../mocks/tests/bins'\n"
     ">>> pimax=40.0\n"
     ">>> (DDrppi, time) = countpairs_rp_pi_mocks(autocorr, cosmology, nthreads, pimax, binfile,\n"
     "                                            ra,dec,cz,ra,dec,cz,\n"
     "                                            verbose=True)\n"
     "\n"
    },
    {"countpairs_s_mu_mocks"       ,(PyCFunction)(void(*)(void)) countpairs_countpairs_s_mu_mocks ,METH_VARARGS | METH_KEYWORDS,
     "countpairs_s_mu_mocks(autocorr, cosmology, nthreads, mu_max, nmu_bins, binfile,\n"
     "                       RA1, DEC1, CZ1, weights1=None, weight_type=None,\n"
     "                       RA2=None, DEC2=None, CZ2=None, weights2=None,\n"
     "                       is_comoving_dist=False,\n"
     "                       verbose=False, output_savg=False,\n"
     "                       fast_divide_and_NR_steps=0, xbin_refine_factor=2, \n"
     "                       ybin_refine_factor=2, zbin_refine_factor=1, \n"
     "                       max_cells_per_dim=100, copy_particles=True, \n"
     "                       enable_min_sep_opt=True, c_api_timer=False, isa=-1)\n\n"

     "Calculate the 2-D pair-counts, "XI_CHAR"(s, "MU_CHAR"), auto/cross-correlation function given two\n"
     "sets of RA1/DEC1/CZ1 and RA2/DEC2/CZ2 arrays. This module is suitable for mock catalogs that have been\n"
     "created by carving out a survey footprint from simulated data. The module can also be used for actual\n"
     "observed galaxies, but you probably want to attach weights to the points to account for completeness etc.\n"
     "\n"
     UNICODE_WARNING
     "\n"
     "Parameters\n"
     "----------\n"
     "Every parameter can be passed as a keyword of the corresponding name.\n\n"

     "autocorr: boolean\n"
     "    Flag for auto/cross-correlation. If autocorr is not 0, the RA2/DEC2/CZ2 arrays\n"
     "    are not used (but must still be passed, as RA1/DEC1/CZ1).\n\n"

     "cosmology: integer\n"
     "    Integer to select cosmology. Pre-set values for (1,2) \n"
     "    1 -> LasDamas cosmology. Om=0.25,  Ol=0.75  (other values are not used)\n"
     "    2 -> Planck   cosmology. Om=0.302, Ol=0.698 \n"
     "    To setup a new cosmology, add an entry to the function, `init_cosmology` in \n"
     "    `ROOT/utils/cosmology_params.c` and recompile the package.\n\n"

     "nthreads: integer\n"
     "    The number of OpenMP threads to use. Has no effect if OpenMP was not used\n"
     "    during library compilation. \n\n"

     "mu_max: double \n"
     "    The maximum mu value to use; must be > 0 and <= 1.0\n\n"

     "nmu_bins: int \n"
     "    The number of "MU_CHAR" bins to use, binning from [0.0, mu_max)\n\n"

     "binfile: filename\n"
     "    Filename containing the radial bins for the correlation function. The file\n"
     "    is expected to contain white-space separated ``smin  smax`` with the bin\n"
     "    edges.  Units must be Mpc/h (see the ``bins`` file in the tests directory\n"
     "    for a sample). For usual logarithmic bins, ``logbins``in the root directory\n"
     "    of this package will create a compatible ``binfile``.\n\n"

     "RA1: array-like, float/double (default double)\n"
     "    The right-ascension of the galaxy, in the range [0, 360]. If there are\n"
     "    negative RA's in the supplied array (input RA in the range [-180, 180]),\n"
     "    then the code will shift the entire array by 180 to put RA's in the\n"
     "    [0, 360] range.\n\n"

     "DEC1: array-like, float/double (default double)\n"
     "    The declination of the galaxy, in the range [-90, 90]. If there are\n"
     "    declinations > 90 in the supplied array (input dec in the range [0, 180]),\n"
     "    then the code will shift the entire array by -90 to put declinations in\n"
     "    the [-90, 90] range. If the code finds declinations more than 180, then\n"
     "    it assumes RA and DEC have been swapped and aborts with that message.\n\n"

     "CZ1: array-like, float/double (default double)\n"
     "    The redshift multiplied by speed of light for the galaxies. The code will\n"
     "    checks that cz has been supplied by comparing with a threshold (currently\n"
     "    set to 10, defined in function check_ra_dec_cz in file\n"
     "    `DDrppi/countpairs_rp_pi_mocks_impl.c.src`) and multiplies by the speed of light if\n"
     "    max z is less than that threshold. If you really want to change the speed\n"
     "    of light, then edit the macro in `ROOT/utils/set_cosmo_dist.h`.\n\n"

     "weights1 : array-like, real (float/double), shape (n_weights_per_particle,n_particles), optional\n"
     "    Weights for computing a weighted pair count.\n\n"

     "weight_type : str, optional\n"
     "    The type of pair weighting to apply.\n"
     "    Options: \"pair_product\", None\n"
     "    Default: None.\n\n"

     "RA2/DEC2/CZ2: float/double (default double)\n"
     "    Same as for RA1/DEC1/CZ1\n\n"

     "weights2\n : array-like, real (float/double), shape (n_weights_per_particle,n_particles), optional\n"
     "    Weights for computing a weighted pair count.\n\n"

     "is_comoving_dist: boolean (default false)\n"
     "    Boolean flag to indicate that ``cz`` values have already been\n"
     "    converted into co-moving distances. This flag allows arbitrary\n"
     "    cosmologies to be used in ``Corrfunc``.\n\n"

     "verbose : boolean (default false)\n"
     "    Boolean flag to control output of informational messages\n\n"

     "output_savg : boolean (default false)\n"
     "    Boolean flag to output the average ``s`` for each bin. Code will\n"
     "    run slightly slower if you set this flag. Also, note, if you are calculating\n"
     "    in single-precision, ``savg`` will suffer from numerical loss of\n"
     "    precision and can not be trusted. If you need accurate ``savg``\n"
     "    values, then pass in double precision arrays for the particle positions.\n\n"

     "fast_divide_and_NR_steps: integer (default 0)\n"
     "    Replaces the division in ``AVX512F`` and ``AVX`` implementation with an\n"
     "    approximate reciprocal, followed by ``fast_divide_and_NR_steps`` "
     "    Newton-Raphson step. Can improve \n"
     "    runtime by ~15-20%. Value of 0 keeps the standard division.\n\n"

     "(xyz)bin_refine_factor: integer (default (2,2,1) typical values in [1-3]) \n"
     "    Controls the refinement on the cell sizes. Can have up to a 20% impact \n"
     "    on runtime. \n\n"

     "max_cells_per_dim: integer (default 100, typical values in [50-300]) \n"
     "    Controls the maximum number of cells per dimension. Total number of cells \n"
     "    can be up to (max_cells_per_dim)^3. Only increase if ``rmax`` is too small \n"
     "    relative to the boxsize (and increasing helps the runtime).\n\n"

     "copy_particles: boolean (default True)\n"
     "    Boolean flag to make a copy of the particle positions\n"
     "    If set to False, the particles will be re-ordered in-place\n\n"

     "enable_min_sep_opt: boolean (default true)\n"
     "    Boolean flag to allow optimizations based on min. separation between pairs \n"
     "    of cells. Here to allow for comparison studies.\n\n"

     "c_api_timer : boolean (default false)\n"
     "    Boolean flag to measure actual time spent in the C libraries. Here\n"
     "    to allow for benchmarking and scaling studies.\n\n"

     "isa : integer (default -1)\n"
     "    Controls the runtime dispatch for the instruction set to use. Possible\n"
     "    options are: [-1, AVX512F, AVX, SSE42, FALLBACK]\n\n"
     "    Setting isa to -1 will pick the fastest available instruction\n"
     "    set on the current computer. However, if you set ``isa`` to, say,\n"
     "    ``AVX`` and ``AVX`` is not available on the computer, then the code will\n"
     "    revert to using ``FALLBACK`` (even though ``SSE42`` might be available).\n"
     "\n"
     "    Unless you are benchmarking the different instruction sets, you should\n"
     "    always leave ``isa`` to the default value. And if you *are* benchmarking,\n"
     "    then the integer values correspond to the ``enum`` for the instruction set\n"
     "    defined in ``utils/defs.h``.\n\n"

     "Returns\n"
     "--------\n\n"

     "A tuple (results, time) \n\n"

     "results : Python list \n"
     "    Contains [smin, smax, savg, "MU_CHAR", npairs, weightavg] \n"
     "    for each "MU_CHAR"-bin (up to 1.0) for each radial bin specified in\n"
     "    the ``binfile``.  If ``output_savg`` is not set, then\n"
     "    ``savg`` will be set to 0.0 for all bins; similarly for ``weightavg``.\n"
     "    ``npairs`` contains the number of pairs in that bin and can be used to \n"
     "    compute the actual "XI_CHAR"(s, "MU_CHAR") by combining  DD, (DR) and RR counts. \n"
     "    See ``~Corrfunc.utils.convert_3d_counts_to_cf``\n\n"

     "time : double\n"
     "    if ``c_api_timer`` is set, then the return value contains the time spent\n"
     "    in the API; otherwise time is set to 0.0\n\n"

     "Example\n"
     "-------\n"
     ">>> import numpy as np\n"
     ">>> from Corrfunc._countpairs_mocks import countpairs_s_mu_mocks\n"
     ">>> ra,dec,cz = np.genfromtxt('../mocks/tests/data/Mr19_mock_northonly.rdcz.dat',dtype=np.float,unpack=True)\n"
     ">>> cosmology=1\n"
     ">>> autocorr=1\n"
     ">>> nthreads=4\n"
     ">>> binfile='../mocks/tests/bins'\n"
     ">>> nmu_bins=10\n"
     ">>> mu_max=1.0\n"
     ">>> (DDsmu, time) = countpairs_s_mu_mocks(autocorr, cosmology, nthreads, mu_max, nmu_bins, binfile,\n"
     "                                            ra,dec,cz,ra,dec,cz,\n"
     "                                            verbose=True)\n"
     "\n"
    },
    {"countpairs_theta_mocks"       ,(PyCFunction)(void(*)(void)) countpairs_countpairs_theta_mocks ,METH_VARARGS | METH_KEYWORDS,
     "countpairs_theta_mocks(autocorr, nthreads, binfile,\n"
     "                       RA1, DEC1, weights1=None, weight_type=None,\n"
     "                       RA2=None, DEC2=None, weights2=None,\n"
     "                       link_in_dec=True, link_in_ra=True,\n"
     "                       verbose=False, output_thetaavg=False,\n"
     "                       fast_acos=False, ra_refine_factor=2,\n"
     "                       dec_refine_factor=2, max_cells_per_dim=100, \n"
     "                       copy_particles=True, enable_min_sep_opt=True, \n"
     "                       c_api_timer=False, isa='fastest')\n"
     "\n"
     "Calculate the angular pair-counts, required for "OMEGA_CHAR"("THETA_CHAR"), auto/cross-correlation function given two\n"
     "sets of RA1/DEC1 and RA2/DEC2 arrays. This module is suitable for mock catalogs that have been\n"
     "created by carving out a survey footprint from simulated data. The module can also be used for actual\n"
     "observed galaxies, but you probably want to attach weights to the points to account for completeness etc.\n"
     "(i.e., seriously consider using some other code for real survey data).\n"
     "\n"
     "Note, that this module only returns pair counts and not the actual correlation function\n"
     ""OMEGA_CHAR"("THETA_CHAR"). See  `~Corrfunc.utils.convert_3d_counts_to_cf.py` for computing\n"
     ""OMEGA_CHAR"("THETA_CHAR") from DD("THETA_CHAR").\n"
     UNICODE_WARNING
     "\n"
     "Parameters\n"
     "-----------\n"
     "Every parameter can be passed as a keyword of the corresponding name.\n\n"

     "autocorr: boolean\n"
     "    Flag for auto/cross-correlation. If autocorr is not 0, the RA2/DEC2 arrays\n"
     "    are not used (but must still be passed as valid arrays).\n\n"

     "nthreads: integer\n"
     "    The number of OpenMP threads to use. Has no effect if OpenMP was not used\n"
     "    during library compilation. \n\n"

     "binfile: filename\n"
     "    Filename containing the radial bins for the correlation function. The file\n"
     "    is expected to contain white-space separated ``thetamin  thetamax`` with the bin\n"
     "    edges.  Units must be degrees (see the ``angular_bins`` file in the tests directory\n"
     "    for a sample). For usual logarithmic bins, ``logbins``in the root directory\n"
     "    of this package will create a compatible ``binfile``.\n\n"

     "RA1: float/double (default double)\n"
     "    The right-ascension of the galaxy, in the range [0, 360]. If there are\n"
     "    negative RA's in the supplied array (input RA in the range [-180, 180]),\n"
     "    then the code will shift the entire array by 180 to put RA's in the\n"
     "    [0, 360] range.\n\n"

     "DEC1: float/double (default double)\n"
     "    The declination of the galaxy, in the range [-90, 90]. If there are\n"
     "    declinations > 90 in the supplied array (input dec in the range [0, 180]),\n"
     "    then the code will shift the entire array by -90 to put declinations in\n"
     "    the [-90, 90] range. If the code finds declinations more than 180, then\n"
     "    it assumes RA and DEC have been swapped and aborts with that message.\n\n"

     "weights1 : array-like, real (float/double), shape (n_weights_per_particle,n_particles), optional\n"
     "    Weights for computing a weighted pair count.\n\n"

     "weight_type : str, optional\n"
     "    The type of pair weighting to apply.\n"
     "    Options: \"pair_product\", None\n"
     "    Default: None.\n\n"

     "RA2/DEC2: float/double (default double)\n"
     "    Same as for RA1/DEC1\n\n"

     "weights2\n : array-like, real (float/double), shape (n_weights_per_particle,n_particles), optional\n"
     "    Weights for computing a weighted pair count.\n\n"

     "verbose : boolean (default false)\n"
     "    Boolean flag to control output of informational messages\n\n"

     "output_thetaavg : boolean (default false)\n"
     "    Boolean flag to output the average ``\theta`` for each bin. Code will\n"
     "    run slower if you set this flag. Also, note, if you are calculating\n"
     "    in single-precision, ``thetaavg`` will suffer from numerical loss of\n"
     "    precision and can not be trusted. If you need accurate ``thetaavg``\n"
     "    values, then pass in double precision arrays for ``RA/DEC``.\n"
     "    **Note** Code will run significantly slower if you enable this option.\n"
     "    Use ``fast_acos`` if you can tolerate some loss of precision.\n\n"

     "fast_acos: boolean (default false)\n"
     "    Flag to use numerical approximation for the ``arccos`` - gives better\n"
     "    performance at the expense of some precision. Relevant only if\n"
     "    ``output_thetaavg==True``.\n"
     "    Developers: Two versions already coded up in ``utils/fast_acos.h``,\n"
     "    so you can choose the version you want. There are also notes on how\n"
     "    to implement faster (and less accurate) functions, particularly relevant\n"
     "    if you know your ``theta`` range is limited. If you implement a new\n"
     "    version, then you will have to reinstall the entire Corrfunc package.\n\n"
     "    **Note** that tests will fail if you run the tests with``fast_acos=True``.\n\n"

     "(radec)_refine_factor: integer (default (2,2) typical values in [1-3]) \n"
     "    Controls the refinement on the cell sizes. Can have up to a 20% impact \n"
     "    on runtime. Note, only two refine factors are to be specified and these \n"
     "    correspond to ``ra`` and ``dec`` (rather, than the usual three of \n"
     "    ``(xyz)bin_refine_factor`` for all other correlation functions).\n\n"

     "max_cells_per_dim: integer (default 100, typical values in [50-300]) \n"
     "    Controls the maximum number of cells per dimension. Total number of cells \n"
     "    can be up to (max_cells_per_dim)^3. Only increase if ``rmax`` is too small \n"
     "    relative to the boxsize (and increasing helps the runtime).\n\n"

     "copy_particles: boolean (default True)\n"
     "    Boolean flag to make a copy of the particle positions\n"
     "    If set to False, the particles will be re-ordered in-place\n\n"

     "enable_min_sep_opt: boolean (default true)\n"
     "    Boolean flag to allow optimizations based on min. separation between pairs \n"
     "    of cells. Here to allow for comparison studies.\n\n"

     "c_api_timer : boolean (default false)\n"
     "    Boolean flag to measure actual time spent in the C libraries. Here\n"
     "    to allow for benchmarking and scaling studies.\n\n"

     "isa : integer (default -1)\n"
     "    Controls the runtime dispatch for the instruction set to use. Possible\n"
     "    options are: [-1, AVX512F, AVX, SSE42, FALLBACK]\n\n"
     "    Setting isa to -1 will pick the fastest available instruction\n"
     "    set on the current computer. However, if you set ``isa`` to, say,\n"
     "    ``AVX`` and ``AVX`` is not available on the computer, then the code will\n"
     "    revert to using ``FALLBACK`` (even though ``SSE42`` might be available).\n"
     "\n"
     "    Unless you are benchmarking the different instruction sets, you should\n"
     "    always leave ``isa`` to the default value. And if you *are* benchmarking,\n"
     "    then the integer values correspond to the ``enum`` for the instruction set\n"
     "    defined in ``utils/defs.h``.\n\n"

     "Returns\n"
     "--------\n"
     "A tuple (results, time) \n\n"

     "results : Python list \n"
     "    Contains [thetamin, thetamax, thetaavg, npairs, weightavg] for each angular bin \n"
     "    specified in the ``binfile``. If ``output_thetaavg`` is not set, then\n"
     "    ``thetaravg`` will be set to 0.0 for all bins; similarly for ``weightavg``. ``npairs`` contains the number of\n"
     "    pairs in that bin and can be used to compute the actual "OMEGA_CHAR"("THETA_CHAR")\n"
     "    by combining  DD, (DR) and RR counts. See ``~Corrfunc.utils.convert_3d_counts_to_cf``\n\n"

     "time : double\n"
     "    if ``c_api_timer`` is set, then the return value contains the time spent\n"
     "    in the API; otherwise time is set to 0.0\n\n"

     "Example\n"
     "-------\n\n"

     ">>> import numpy as np\n"
     ">>> from Corrfunc._countpairs_mocks import countpairs_theta_mocks\n"
     ">>> ra, dec, _ = np.genfromtxt('../mocks/tests/data/Mr19_mock_northonly.rdcz.dat',dtype=np.float,unpack=True)\n"
     ">>> cosmology=1\n"
     ">>> autocorr=1\n"
     ">>> nthreads=4\n"
     ">>> binfile='../mocks/tests/angular_bins'\n"
     ">>> DD=countpairs_theta_mocks(autocorr, cosmology, nthreads, binfile, \n"
     "                              ra,dec,ra,dec,\n"
     "                              verbose=True)\n"
     "\n"
    },
    {"countspheres_vpf_mocks"       ,(PyCFunction)(void(*)(void)) countpairs_countspheres_vpf_mocks ,METH_VARARGS | METH_KEYWORDS,
     "countspheres_vpf_mocks(rmax, nbins, nspheres, numpN,\n"
     "                       threshold_ngb, centers_file, cosmology,\n"
     "                       RA, DEC, CZ,\n"
     "                       RAND_RA, RAND_DEC, RAND_CZ,\n"
     "                       verbose=False, is_comoving_dist=False,\n"
     "                       xbin_refine_factor=2, ybin_refine_factor=2, \n"
     "                       zbin_refine_factor=1, max_cells_per_dim=100, \n"
     "                       copy_particles=1, c_api_timer=False, isa=-1)\n"
     "\n"
     "Calculates the fraction of random spheres that contain exactly *N* points, pN(r).\n"
     "Returns a numpy structured array containing the probability of a sphere of radius\n"
     "up to ``rmax`` containing ``0--numpN-1`` galaxies.\n"
     UNICODE_WARNING
     "\n"
     "Parameters\n"
     "------------\n"
     "Every parameter can be passed as a keyword of the corresponding name.\n\n"

     "rmax : double\n"
     "    Maximum radius of the sphere to place on the particles\n\n"

     "nbins : integer\n"
     "    Number of bins in the counts-in-cells. Radius of first shell\n"
     "    is rmax/nbins\n\n"

     "nspheres: integer (>= 0)\n"
     "    Number of random spheres to place within the particle distribution.\n"
     "    For a small number of spheres, the error is larger in the measured\n"
     "    pN's.\n\n"

     "numpN: integer (>= 1)\n"
     "    Governs how many unique pN's are to returned. If ``numpN`` is set to 1,\n"
     "    then only the vpf (p0) is returned. For ``numpN=2``, p0 and p1 are\n"
     "    returned.\n\n"

     "    More explicitly, the columns in the results look like the following:\n"
     "      numpN = 1 -> p0\n"
     "      numpN = 2 -> p0 p1\n"
     "      numpN = 3 -> p0 p1 p2\n"
     "      and so on...(note that p0 is the vpf).\n\n"

     "threshold_ngb: integer\n"
     "    Minimum number of random points needed in a ``rmax`` sphere such that it\n"
     "    is considered to be entirely within the mock footprint. The\n"
     "    command-line version, ``mocks/vpf/vpf_mocks.c``, assumes that the\n"
     "    minimum number of randoms can be at most a 1-sigma deviation from\n"
     "    the expected random number density.\n\n"

     "centers_file: string, filename\n"
     "    A file containing random sphere centers. If the file does not exist,\n"
     "    then a list of random centers will be written out. In that case, the\n"
     "    randoms arrays, ``RAND_RA``, ``RAND_DEC`` and ``RAND_CZ`` are used to\n"
     "    check that the sphere is entirely within the footprint. If the file does\n"
     "    exist but either ``rmax`` is too small or there are not enough centers\n"
     "    then the file will be overwritten.\n\n"
     "    **Note** If the centers file has to be written, the code will take\n"
     "    significantly longer to finish. However, subsequent runs can re-use\n"
     "    that centers file and will be faster.\n\n"

     "cosmology: integer, required\n"
     "    Integer choice for setting cosmology. Valid values are 1->LasDamas\n"
     "    cosmology and 2->Planck cosmology. If you need arbitrary cosmology,\n"
     "    easiest way is to convert the ``CZ`` values into co-moving distance,\n"
     "    based on your preferred cosmology. Set ``is_comoving_dist=True``, to\n"
     "    indicate that the co-moving distance conversion has already been done.\n"
     "\n"
     "    Choices: 1 -> LasDamas cosmology. Om=0.25,  Ol=0.75\n"
     "             2 -> Planck   cosmology. Om=0.302, Ol=0.698\n"
     "\n"
     "    To setup a new cosmology, add an entry to the function,\n"
     "    ``init_cosmology`` in ``ROOT/utils/cosmology_params.c`` and re-install\n"
     "    the entire package.\n\n"

     "RA: array-like, real (float/double)\n"
     "    The array of Right Ascensions for the first set of points. RA's\n"
     "    are expected to be in [0.0, 360.0], but the code will try to fix cases\n"
     "    where the RA's are in [-180, 180.0]. For peace of mind, always supply\n"
     "    RA's in [0.0, 360.0].\n"
     "\n"
     "    Calculations are done in the precision of the supplied arrays.\n\n"

     "DEC: array-like, real (float/double)\n"
     "    Array of Declinations for the first set of points. DEC's are expected\n"
     "    to be in the [-90.0, 90.0], but the code will try to fix cases where\n"
     "    the DEC's are in [0.0, 180.0]. Again, for peace of mind, always supply\n"
     "    DEC's in [-90.0, 90.0].\n"
     "\n"
     "    Must be of same precision type as RA.\n\n"

     "CZ: array-like, real (float/double)\n"
     "    Array of (Speed Of Light * Redshift) values for the first set of\n"
     "    points. Code will try to detect cases where ``redshifts`` have been\n"
     "    passed and multiply the entire array with the ``speed of light``.\n"
     "\n"
     "    If ``is_comoving_dist`` is set, then ``CZ`` is interpreted as the\n"
     "    co-moving distance, rather than (Speed Of Light * Redshift).\n\n"

     "RAND_RA: array-like, real (float/double)\n"
     "    The array of Right Ascensions for the randoms. RA's are expected to be\n"
     "    in [0.0, 360.0], but the code will try to fix cases where the RA's are\n"
     "    in [-180, 180.0]. For peace of mind, always supply RA's in\n"
     "    [0.0, 360.0].\n"
     "\n"
     "    Must be of same precision type as RA/DEC/CZ.\n\n"

     "RAND_DEC: array-like, real (float/double)\n"
     "    Array of Declinations for the randoms. DEC's are expected to be in the\n"
     "    [-90.0, 90.0], but the code will try to fix cases where the DEC's are\n"
     "    in [0.0, 180.0]. Again, for peace of mind, always supply DEC's in\n"
     "    [-90.0, 90.0].\n"
     "\n"
     "    Must be of same precision type as RA/DEC/CZ.\n\n"

     "RAND_CZ: array-like, real (float/double)\n"
     "    Array of (Speed Of Light * Redshift) values for the randoms. Code\n"
     "    will try to detect cases where ``redshifts`` have been\n"
     "    passed and multiply the entire array with the ``speed of light``.\n"
     "\n"
     "    If ``is_comoving_dist`` is set, then ``CZ2`` is interpreted as the\n"
     "    co-moving distance, rather than ``(Speed Of Light * Redshift)``.\n"
     "\n"
     "    **Note** RAND_RA, RAND_DEC and RAND_CZ are only used when the\n"
     "      ``centers_file``  needs to be written out. In that case, the\n"
     "      RAND_RA, RAND_DEC, and RAND_CZ are used as random centers.\n\n"

     "is_comoving_dist: boolean (default false)\n"
     "    Boolean flag to indicate that ``cz`` values have already been\n"
     "    converted into co-moving distances. This flag allows arbitrary\n"
     "    cosmologies to be used in ``Corrfunc``.\n\n"

     "verbose : boolean (default false)\n"
     "    Boolean flag to control output of informational messages\n\n"

     "(xyz)bin_refine_factor: integer (default (1,1,1) typical values in [1-3]) \n"
     "    Controls the refinement on the cell sizes. Can have up to a 20% impact \n"
     "    on runtime. Note that the default values are different from the \n"
     "    correlation function routines.\n\n"

     "max_cells_per_dim: integer (default 100, typical values in [50-300]) \n"
     "    Controls the maximum number of cells per dimension. Total number of cells \n"
     "    can be up to (max_cells_per_dim)^3. Only increase if ``rmax`` is too small \n"
     "    relative to the boxsize (and increasing helps the runtime).\n\n"

     "copy_particles: boolean (default True)\n"
     "    Boolean flag to make a copy of the particle positions\n"
     "    If set to False, the particles will be re-ordered in-place\n\n"

     "c_api_timer : boolean (default false)\n"
     "    Boolean flag to measure actual time spent in the C libraries. Here\n"
     "    to allow for benchmarking and scaling studies.\n\n"

     "isa : integer (default -1)\n"
     "    Controls the runtime dispatch for the instruction set to use. Possible\n"
     "    options are: [-1, AVX512F, AVX, SSE42, FALLBACK]\n\n"
     "    Setting isa to -1 will pick the fastest available instruction\n"
     "    set on the current computer. However, if you set ``isa`` to, say,\n"
     "    ``AVX`` and ``AVX`` is not available on the computer, then the code will\n"
     "    revert to using ``FALLBACK`` (even though ``SSE42`` might be available).\n"
     "\n"
     "    Unless you are benchmarking the different instruction sets, you should\n"
     "    always leave ``isa`` to the default value. And if you *are* benchmarking,\n"
     "    then the integer values correspond to the ``enum`` for the instruction set\n"
     "    defined in ``utils/defs.h``.\n\n"

     "Returns\n"
     "--------\n"
     "A tuple (results, time) \n\n"

     "results : Python list of lists\n\n"
     "    Contains [rmax, pN[numpN]] with ``nbins`` elements. Each row contains\n"
     "    the maximum radius of the sphere and the ``numpN`` elements in the \n"
     "    ``pN`` array. Each element of this array contains the probability that\n"
     "    a sphere of radius ``rmax`` contains *exactly* ``N`` galaxies. For \n"
     "    example, pN[0] (p0, the void probibility function) is the probability\n"
     "    that a sphere of radius ``rmax`` contains 0 galaxies.\n\n"

     "time : double\n"
     "    if ``c_api_timer`` is set, then the return value contains the time spent\n"
     "    in the API; otherwise time is set to 0.0\n\n"

    "Example\n"
    "--------\n\n"

    ">>> import Corrfunc\n"
    ">>> import math\n"
    ">>> from os.path import dirname, abspath, join as pjoin\n"
    ">>> rmax = 10.0\n"
    ">>> nbins = 10\n"
    ">>> numbins_to_print = nbins\n"
    ">>> nspheres = 10000\n"
    ">>> numpN = 6\n"
    ">>> threshold_ngb = 1  # does not matter since we have the centers\n"
    ">>> cosmology = 1  # LasDamas cosmology\n"
    ">>> centers_file = pjoin(dirname(abspath(Corrfunc.__file__)),\n"
    "                         '../mocks/tests/data/',\n"
    "                         'Mr19_centers_xyz_forVPF_rmax_10Mpc.txt')\n"
    ">>> N = 100000\n"
    ">>> boxsize = 420.0\n"
    ">>> X = np.random.uniform(-0.5*boxsize, 0.5*boxsize, N)\n"
    ">>> Y = np.random.uniform(-0.5*boxsize, 0.5*boxsize, N)\n"
    ">>> Z = np.random.uniform(-0.5*boxsize, 0.5*boxsize, N)\n"
    ">>> CZ = np.sqrt(X*X + Y*Y + Z*Z)\n"
    ">>> inv_cz = 1.0/CZ\n"
    ">>> X *= inv_cz\n"
    ">>> Y *= inv_cz\n"
    ">>> Z *= inv_cz\n"
    ">>> DEC = 90.0 - np.arccos(Z)*180.0/math.pi\n"
    ">>> RA = (np.arctan2(Y, X)*180.0/math.pi) + 180.0\n"
    ">>> results = vpf_mocks(rmax, nbins, nspheres, numpN,\n"
    "                        threshold_ngb, centers_file, cosmology,\n"
    "                        RA, DEC, CZ,\n"
    "                        RA, DEC, CZ,\n"
    "                        verbose=True,\n"
    "                        is_comoving_dist=True)\n"
    "\n"
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
    PyErr_Print();
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

    highest_isa_mocks = get_max_usable_isa();

#if PY_MAJOR_VERSION >= 3
    return module;
#endif

}

static int print_kwlist_into_msg(char *msg, const size_t totsize, size_t len, char *kwlist[], const size_t nitems)
{
    for(size_t i=0;i<nitems;i++) {

        if(len+strlen(kwlist[i]) >= totsize-2) {
            return EXIT_FAILURE;
        }

        memcpy(msg+len, kwlist[i], strlen(kwlist[i]));
        len += strlen(kwlist[i]);
        msg[len] = ',';
        msg[len+1] = ' ';
        len += 2;
    }

    msg[len]='\0';
    return EXIT_SUCCESS;
}


// weights1_obj may be NULL, in which case it is ignored.
// If it is not NULL, it will be checked alongside the positions
static int64_t check_dims_and_datatype(PyObject *module, PyArrayObject *x1_obj, PyArrayObject *y1_obj, PyArrayObject *z1_obj, PyArrayObject *weights1_obj, size_t *element_size)
{
    char msg[1024];

    const int check_weights = weights1_obj != NULL;

    /* All the position arrays should be 1-D*/
    const int nxdims = PyArray_NDIM(x1_obj);
    const int nydims = PyArray_NDIM(y1_obj);
    const int nzdims = PyArray_NDIM(z1_obj);

    if(nxdims != 1 || nydims != 1 || nzdims != 1) {
        snprintf(msg, 1024, "ERROR: Expected 1-D numpy arrays.\nFound (nxdims, nydims, nzdims) = (%d, %d, %d) instead",
                 nxdims, nydims, nzdims);
        countpairs_mocks_error_out(module, msg);
        return -1;
    }

    /* The weights array can be 1-D or 2-D of shape (n_weights, n_particles) */
    const int n_weight_dims = check_weights ? PyArray_NDIM(weights1_obj) : 2;

    if(n_weight_dims != 2) {
        snprintf(msg, 1024, "ERROR: Expected 2-D weight array of shpae (n_weights_per_particle,n_particles).\nFound n_weight_dims = %d instead", n_weight_dims);
        countpairs_mocks_error_out(module, msg);
        return -1;
    }

    /* All the arrays should be floating point (only float32 and float64 are allowed) */
    const int x_type = PyArray_TYPE(x1_obj);
    const int y_type = PyArray_TYPE(y1_obj);
    const int z_type = PyArray_TYPE(z1_obj);
    const int weights_type = check_weights ? PyArray_TYPE(weights1_obj) : NPY_NOTYPE;
    if( ! ((x_type == NPY_FLOAT || x_type == NPY_DOUBLE) &&
           (y_type == NPY_FLOAT || y_type == NPY_DOUBLE) &&
           (z_type == NPY_FLOAT || z_type == NPY_DOUBLE) &&
           (!check_weights || weights_type == NPY_FLOAT || weights_type == NPY_DOUBLE))
        ) {
        PyArray_Descr *x_descr = PyArray_DescrFromType(x_type);
        PyArray_Descr *y_descr = PyArray_DescrFromType(y_type);
        PyArray_Descr *z_descr = PyArray_DescrFromType(z_type);
        PyArray_Descr *weights_descr = PyArray_DescrFromType(weights_type);
        if(x_descr == NULL || y_descr == NULL || z_descr == NULL || weights_descr == NULL) {
            /* Generating the dtype descriptor failed somehow. At least provide some information */
            snprintf(msg, 1024, "TypeError: Expected floating point arrays (allowed types = %d or %d). Instead found type-nums (%d, %d, %d, %d)\n",
                     NPY_FLOAT, NPY_DOUBLE, x_type, y_type, z_type, weights_type);
        } else {
            snprintf(msg, 1024, "TypeError: Expected floating point arrays (allowed types = %d or %d). Instead found type-nums (%d, %d, %d, %d) "
                     "with type-names = (%s, %s, %s, %s)\n",
                     NPY_FLOAT, NPY_DOUBLE, x_type, y_type, z_type, weights_type, x_descr->typeobj->tp_name, y_descr->typeobj->tp_name, z_descr->typeobj->tp_name, weights_descr->typeobj->tp_name);
        }
        Py_XDECREF(x_descr);Py_XDECREF(y_descr);Py_XDECREF(z_descr);Py_XDECREF(weights_descr);
        countpairs_mocks_error_out(module, msg);
        return -1;
    }

    // Current version of the code only supports weights of the same dtype as positions
    if( x_type != y_type || y_type != z_type || (check_weights && z_type != weights_type)) {
        PyArray_Descr *x_descr = PyArray_DescrFromType(x_type);
        PyArray_Descr *y_descr = PyArray_DescrFromType(y_type);
        PyArray_Descr *z_descr = PyArray_DescrFromType(z_type);
        PyArray_Descr *weights_descr = PyArray_DescrFromType(weights_type);
        if(x_descr == NULL || y_descr == NULL || z_descr == NULL || weights_descr == NULL) {
            /* Generating the dtype descriptor failed somehow. At least provide some information */
            snprintf(msg, 1024, "TypeError: Expected *ALL* 3 floating point arrays to be the same type (allowed types = %d or %d). Instead found type-nums (%d, %d, %d, %d)\n",
                     NPY_FLOAT, NPY_DOUBLE, x_type, y_type, z_type, weights_type);
        } else {
            snprintf(msg, 1024, "TypeError: Expected *ALL* 3 floating point arrays to be the same type (allowed types = %d or %d). Instead found type-nums (%d, %d, %d, %d) "
                     "with type-names = (%s, %s, %s, %s)\n",
                     NPY_FLOAT, NPY_DOUBLE, x_type, y_type, z_type, weights_type, x_descr->typeobj->tp_name, y_descr->typeobj->tp_name, z_descr->typeobj->tp_name, weights_descr->typeobj->tp_name);
        }
        Py_XDECREF(x_descr);Py_XDECREF(y_descr);Py_XDECREF(z_descr);Py_XDECREF(weights_descr);
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

    // The last dimension of the weights array must match the number of positions
    if(check_weights){
        const int64_t n_weights1 = (int64_t) PyArray_DIMS(weights1_obj)[n_weight_dims-1];
        if(nx1 != n_weights1){
            snprintf(msg, 1024, "ERROR: the last dimension of `weights` must match the number of positions.  Instead found n_weights=%"PRId64", nx=%"PRId64,
                   n_weights1, nx1);
            countpairs_mocks_error_out(module, msg);
            return -1;
        }
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
    PyArrayObject *x1_obj=NULL, *y1_obj=NULL, *z1_obj=NULL, *weights1_obj=NULL;
    PyArrayObject *x2_obj=NULL, *y2_obj=NULL, *z2_obj=NULL, *weights2_obj=NULL;

    struct config_options options = get_config_options();
    options.is_comoving_dist = 0;
    options.verbose = 0;
    options.instruction_set = -1;
    options.periodic = 0;
    options.fast_divide_and_NR_steps=0;
    options.c_api_timer = 0;
    options.enable_min_sep_opt = 1;
    options.copy_particles = 1;

    int8_t xbin_ref=options.bin_refine_factors[0],
        ybin_ref=options.bin_refine_factors[1],
        zbin_ref=options.bin_refine_factors[2];

    int autocorr=1;
    int nthreads=4;
    int cosmology=1;
    double pimax;
    char *binfile, *weighting_method_str = NULL;

    static char *kwlist[] = {
        "autocorr",
        "cosmology",
        "nthreads",
        "pimax",
        "binfile",
        "RA1",
        "DEC1",
        "CZ1",
        "weights1",
        "RA2",
        "DEC2",
        "CZ2",
        "weights2",
        "is_comoving_dist",
        "verbose", /* keyword verbose -> print extra info at runtime + progressbar */
        "output_rpavg",
        "fast_divide_and_NR_steps",
        "xbin_refine_factor",
        "ybin_refine_factor",
        "zbin_refine_factor",
        "max_cells_per_dim",
        "copy_particles",
        "enable_min_sep_opt",
        "c_api_timer",
        "isa",/* instruction set to use of type enum isa; valid values are AVX512F, AVX, SSE, FALLBACK (enum) */
        "weight_type",
        NULL
    };

    if ( ! PyArg_ParseTupleAndKeywords(args, kwargs, "iiidsO!O!O!|O!O!O!O!O!bbbbbbbhbbbis", kwlist,
                                       &autocorr,&cosmology,&nthreads,&pimax,&binfile,
                                       &PyArray_Type,&x1_obj,
                                       &PyArray_Type,&y1_obj,
                                       &PyArray_Type,&z1_obj,
                                       &PyArray_Type,&weights1_obj,
                                       &PyArray_Type,&x2_obj,//optional parameters -> if autocorr == 1, not checked; required if autocorr=0
                                       &PyArray_Type,&y2_obj,
                                       &PyArray_Type,&z2_obj,
                                       &PyArray_Type,&weights2_obj,
                                       &(options.is_comoving_dist),
                                       &(options.verbose),
                                       &(options.need_avg_sep),
                                       &(options.fast_divide_and_NR_steps),
                                       &xbin_ref, &ybin_ref, &zbin_ref,
                                       &(options.max_cells_per_dim),
                                       &(options.copy_particles),
                                       &(options.enable_min_sep_opt),
                                       &(options.c_api_timer),
                                       &(options.instruction_set),
                                       &weighting_method_str)

         ) {

        PyObject_Print(kwargs, stdout, 0);
        fprintf(stdout, "\n");

        char msg[1024];
        int len=snprintf(msg, 1024,"ArgumentError: In DDrppi_mocks> Could not parse the arguments. Input parameters are: \n");

        /* How many keywords do we have? Subtract 1 because of the last NULL */
        const size_t nitems = sizeof(kwlist)/sizeof(*kwlist) - 1;
        int status = print_kwlist_into_msg(msg, 1024, len, kwlist, nitems);
        if(status != EXIT_SUCCESS) {
            fprintf(stderr,"Error message does not contain all of the keywords\n");
        }

        countpairs_mocks_error_out(module,msg);

        Py_RETURN_NONE;
    }

    /*This is for the fastest isa */
    if(options.instruction_set == -1) {
        options.instruction_set = highest_isa_mocks;
    }
    if(xbin_ref != options.bin_refine_factors[0] ||
       ybin_ref != options.bin_refine_factors[1] ||
       zbin_ref != options.bin_refine_factors[2]) {
        options.bin_refine_factors[0] = xbin_ref;
        options.bin_refine_factors[1] = ybin_ref;
        options.bin_refine_factors[2] = zbin_ref;
        set_bin_refine_scheme(&options, BINNING_CUST);//custom binning -> code will honor requested binning scheme
    }

    /* Validate the user's choice of weighting method */
    weight_method_t weighting_method;
    int wstatus = get_weight_method_by_name(weighting_method_str, &weighting_method);
    if(wstatus != EXIT_SUCCESS){
        char msg[1024];
        snprintf(msg, 1024, "ValueError: In %s: unknown weight_type \"%s\"!", __FUNCTION__, weighting_method_str);
        countpairs_mocks_error_out(module, msg);
        Py_RETURN_NONE;
    }
    if(weighting_method == NONE){
        // Do not attempt to validate the weights array if it will not be used!
        weights1_obj = NULL;
        weights2_obj = NULL;
    }

    /* We have numpy arrays and all the required inputs*/
    /* How many data points are there? And are they all of floating point type */
    size_t element_size;
    const int64_t ND1 = check_dims_and_datatype(module, x1_obj, y1_obj, z1_obj, weights1_obj, &element_size);
    if(ND1 == -1) {
        //Error has already been set -> simply return
        Py_RETURN_NONE;
    }

    int found_weights = weights1_obj == NULL ? 0 : PyArray_SHAPE(weights1_obj)[0];
    struct extra_options extra = get_extra_options(weighting_method);
    if(extra.weights0.num_weights > 0 && extra.weights0.num_weights != found_weights){
        char msg[1024];
        snprintf(msg, 1024, "ValueError: In %s: specified weighting method %s which requires %"PRId64" weight(s)-per-particle, but found %d weight(s) instead!\n",
                 __FUNCTION__, weighting_method_str, extra.weights0.num_weights, found_weights);
        countpairs_mocks_error_out(module, msg);
        Py_RETURN_NONE;
    }

    if(extra.weights0.num_weights > 0 && found_weights > MAX_NUM_WEIGHTS){
        char msg[1024];
        snprintf(msg, 1024, "ValueError: In %s: Provided %d weights-per-particle, but the code was compiled with MAX_NUM_WEIGHTS=%d.\n",
                 __FUNCTION__, found_weights, MAX_NUM_WEIGHTS);
        countpairs_mocks_error_out(module, msg);
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
        if((weights1_obj == NULL) != (weights2_obj == NULL)){
            snprintf(msg, 1024, "ValueError: In %s: If autocorr is 0, must pass either zero or two sets of weights.\n",
                     __FUNCTION__);
            countpairs_mocks_error_out(module, msg);
            Py_RETURN_NONE;
        }

        size_t element_size2;
        ND2 = check_dims_and_datatype(module, x2_obj, y2_obj, z2_obj, weights2_obj, &element_size2);
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
    PyObject *x1_array = PyArray_FromArray(x1_obj, NOTYPE_DESCR, requirements);
    PyObject *y1_array = PyArray_FromArray(y1_obj, NOTYPE_DESCR, requirements);
    PyObject *z1_array = PyArray_FromArray(z1_obj, NOTYPE_DESCR, requirements);
    PyObject *weights1_array = NULL;
    if(weights1_obj != NULL){
        weights1_array = PyArray_FromArray(weights1_obj, NOTYPE_DESCR, requirements);
    }

    PyObject *x2_array = NULL, *y2_array = NULL, *z2_array = NULL, *weights2_array = NULL;
    if(autocorr == 0) {
        x2_array = PyArray_FromArray(x2_obj, NOTYPE_DESCR, requirements);
        y2_array = PyArray_FromArray(y2_obj, NOTYPE_DESCR, requirements);
        z2_array = PyArray_FromArray(z2_obj, NOTYPE_DESCR, requirements);
        if(weights2_obj != NULL){
            weights2_array = PyArray_FromArray(weights2_obj, NOTYPE_DESCR, requirements);
        }
    }

    if (x1_array == NULL || y1_array == NULL || z1_array == NULL ||
        (autocorr == 0 && (x2_array == NULL || y2_array == NULL || z2_array == NULL))) {
        Py_XDECREF(x1_array);
        Py_XDECREF(y1_array);
        Py_XDECREF(z1_array);
        Py_XDECREF(weights1_array);

        Py_XDECREF(x2_array);
        Py_XDECREF(y2_array);
        Py_XDECREF(z2_array);
        Py_XDECREF(weights2_array);
        char msg[1024];
        snprintf(msg, 1024, "TypeError: In %s: Could not convert input to arrays of allowed floating point types (doubles or floats). Are you passing numpy arrays?",
                 __FUNCTION__);
        countpairs_mocks_error_out(module, msg);
        Py_RETURN_NONE;
    }

    /* Get pointers to the data as C-types. */
    void *phiD1   = PyArray_DATA((PyArrayObject *)x1_array);
    void *thetaD1 = PyArray_DATA((PyArrayObject *)y1_array);
    void *czD1    = PyArray_DATA((PyArrayObject *)z1_array);
    void *weights1=NULL;
    if(weights1_array != NULL){
        weights1 = PyArray_DATA((PyArrayObject *) weights1_array);
    }

    void *phiD2=NULL, *thetaD2=NULL, *czD2=NULL, *weights2=NULL;
    if(autocorr == 0) {
        phiD2   = PyArray_DATA((PyArrayObject *) x2_array);
        thetaD2 = PyArray_DATA((PyArrayObject *) y2_array);
        czD2    = PyArray_DATA((PyArrayObject *) z2_array);
        if(weights2_array != NULL){
            weights2 = PyArray_DATA((PyArrayObject *) weights2_array);
        }
    }
    options.float_type = element_size;

    /* Pack the weights into extra_options */
    for(int64_t w = 0; w < extra.weights0.num_weights; w++){
        extra.weights0.weights[w] = (char *) weights1 + w*ND1*element_size;
        if(autocorr == 0){
            extra.weights1.weights[w] = (char *) weights2 + w*ND2*element_size;
        }
    }

    NPY_BEGIN_THREADS_DEF;
    NPY_BEGIN_THREADS;

    results_countpairs_mocks results;
    double c_api_time = 0.0;
    int status = countpairs_mocks(ND1,phiD1,thetaD1,czD1,
                                  ND2,phiD2,thetaD2,czD2,
                                  nthreads,
                                  autocorr,
                                  binfile,
                                  pimax,
                                  cosmology,
                                  &results,
                                  &options,
                                  &extra);
    if(options.c_api_timer) {
        c_api_time = options.c_api_time;
    }
    NPY_END_THREADS;

    /* Clean up. */
    Py_DECREF(x1_array);Py_DECREF(y1_array);Py_DECREF(z1_array);Py_XDECREF(weights1_array);//x1 should absolutely not be NULL
    Py_XDECREF(x2_array);Py_XDECREF(y2_array);Py_XDECREF(z2_array);Py_XDECREF(weights2_array);//x2 might be NULL depending on value of autocorr

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
            const double rpavg = results.rpavg[bin_index];
            const double weight_avg = results.weightavg[bin_index];
            PyObject *item = Py_BuildValue("(ddddkd)", rlow,results.rupp[i],rpavg,(j+1)*dpi,results.npairs[bin_index], weight_avg);
            PyList_Append(ret, item);
            Py_XDECREF(item);
        }
        rlow=results.rupp[i];
    }
    free_results_mocks(&results);

    PyObject *rettuple = Py_BuildValue("(Od)", ret, c_api_time);
    Py_DECREF(ret);  // transfer reference ownership to the tuple
    return rettuple;
}

static PyObject *countpairs_countpairs_s_mu_mocks(PyObject *self, PyObject *args, PyObject *kwargs)
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
    PyArrayObject *x1_obj=NULL, *y1_obj=NULL, *z1_obj=NULL, *weights1_obj=NULL;
    PyArrayObject *x2_obj=NULL, *y2_obj=NULL, *z2_obj=NULL, *weights2_obj=NULL;

    struct config_options options = get_config_options();
    options.is_comoving_dist = 0;
    options.verbose = 0;
    options.instruction_set = -1;
    options.periodic = 0;
    options.fast_divide_and_NR_steps=0;
    options.enable_min_sep_opt = 1;
    options.copy_particles = 1;
    options.c_api_timer = 0;
    int8_t xbin_ref=options.bin_refine_factors[0],
        ybin_ref=options.bin_refine_factors[1],
        zbin_ref=options.bin_refine_factors[2];

    int autocorr=1;
    int nthreads=4;
    int cosmology=1;
    int nmu_bins=10;
    double mu_max=1.0;
    char *binfile, *weighting_method_str = NULL;

    static char *kwlist[] = {
        "autocorr",
        "cosmology",
        "nthreads",
        "mu_max",
        "nmu_bins",
        "binfile",
        "RA1",
        "DEC1",
        "CZ1",
        "weights1",
        "RA2",
        "DEC2",
        "CZ2",
        "weights2",
        "is_comoving_dist",
        "verbose", /* keyword verbose -> print extra info at runtime + progressbar */
        "output_savg",
        "fast_divide_and_NR_steps",
        "xbin_refine_factor",
        "ybin_refine_factor",
        "zbin_refine_factor",
        "max_cells_per_dim",
        "copy_particles",
        "enable_min_sep_opt",
        "c_api_timer",
        "isa",/* instruction set to use of type enum isa; valid values are AVX512F, AVX, SSE, FALLBACK (enum) */
        "weight_type",
        NULL
    };

    if ( ! PyArg_ParseTupleAndKeywords(args, kwargs, "iiidisO!O!O!|O!O!O!O!O!bbbbbbbhbbbis", kwlist,
                                       &autocorr,&cosmology,&nthreads,&mu_max,&nmu_bins,&binfile,
                                       &PyArray_Type,&x1_obj,
                                       &PyArray_Type,&y1_obj,
                                       &PyArray_Type,&z1_obj,
                                       &PyArray_Type,&weights1_obj,
                                       &PyArray_Type,&x2_obj,//optional parameters -> if autocorr == 1, not checked; required if autocorr=0
                                       &PyArray_Type,&y2_obj,
                                       &PyArray_Type,&z2_obj,
                                       &PyArray_Type,&weights2_obj,
                                       &(options.is_comoving_dist),
                                       &(options.verbose),
                                       &(options.need_avg_sep),
                                       &(options.fast_divide_and_NR_steps),
                                       &xbin_ref, &ybin_ref, &zbin_ref,
                                       &(options.max_cells_per_dim),
                                       &(options.copy_particles),
                                       &(options.enable_min_sep_opt),
                                       &(options.c_api_timer),
                                       &(options.instruction_set),
                                       &weighting_method_str)

         ) {

        PyObject_Print(kwargs, stdout, 0);
        fprintf(stdout, "\n");

        char msg[1024];
        int len=snprintf(msg, 1024,"ArgumentError: In DDsmu_mocks> Could not parse the arguments. Input parameters are: \n");

        /* How many keywords do we have? Subtract 1 because of the last NULL */
        const size_t nitems = sizeof(kwlist)/sizeof(*kwlist) - 1;
        int status = print_kwlist_into_msg(msg, 1024, len, kwlist, nitems);
        if(status != EXIT_SUCCESS) {
            fprintf(stderr,"Error message does not contain all of the keywords\n");
        }

        countpairs_mocks_error_out(module,msg);

        Py_RETURN_NONE;
    }

    /*This is for the fastest isa */
    if(options.instruction_set == -1) {
        options.instruction_set = highest_isa_mocks;
    }
    if(xbin_ref != options.bin_refine_factors[0] ||
       ybin_ref != options.bin_refine_factors[1] ||
       zbin_ref != options.bin_refine_factors[2]) {
        options.bin_refine_factors[0] = xbin_ref;
        options.bin_refine_factors[1] = ybin_ref;
        options.bin_refine_factors[2] = zbin_ref;
        set_bin_refine_scheme(&options, BINNING_CUST);//custom binning -> code will honor requested binning scheme
    }

    /* Validate the user's choice of weighting method */
    weight_method_t weighting_method;
    int wstatus = get_weight_method_by_name(weighting_method_str, &weighting_method);
    if(wstatus != EXIT_SUCCESS){
        char msg[1024];
        snprintf(msg, 1024, "ValueError: In %s: unknown weight_type \"%s\"!", __FUNCTION__, weighting_method_str);
        countpairs_mocks_error_out(module, msg);
        Py_RETURN_NONE;
    }
    if(weighting_method == NONE){
        // Do not attempt to validate the weights array if it will not be used!
        weights1_obj = NULL;
        weights2_obj = NULL;
    }

    /* We have numpy arrays and all the required inputs*/
    /* How many data points are there? And are they all of floating point type */
    size_t element_size;
    const int64_t ND1 = check_dims_and_datatype(module, x1_obj, y1_obj, z1_obj, weights1_obj, &element_size);
    if(ND1 == -1) {
        //Error has already been set -> simply return
        Py_RETURN_NONE;
    }

    int found_weights = weights1_obj == NULL ? 0 : PyArray_SHAPE(weights1_obj)[0];
    struct extra_options extra = get_extra_options(weighting_method);
    if(extra.weights0.num_weights > 0 && extra.weights0.num_weights != found_weights){
        char msg[1024];
        snprintf(msg, 1024, "ValueError: In %s: specified weighting method %s which requires %"PRId64" weight(s)-per-particle, but found %d weight(s) instead!\n",
                 __FUNCTION__, weighting_method_str, extra.weights0.num_weights, found_weights);
        countpairs_mocks_error_out(module, msg);
        Py_RETURN_NONE;
    }

    if(extra.weights0.num_weights > 0 && found_weights > MAX_NUM_WEIGHTS){
        char msg[1024];
        snprintf(msg, 1024, "ValueError: In %s: Provided %d weights-per-particle, but the code was compiled with MAX_NUM_WEIGHTS=%d.\n",
                 __FUNCTION__, found_weights, MAX_NUM_WEIGHTS);
        countpairs_mocks_error_out(module, msg);
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
        if((weights1_obj == NULL) != (weights2_obj == NULL)){
            snprintf(msg, 1024, "ValueError: In %s: If autocorr is 0, must pass either zero or two sets of weights.\n",
                     __FUNCTION__);
            countpairs_mocks_error_out(module, msg);
            Py_RETURN_NONE;
        }

        size_t element_size2;
        ND2 = check_dims_and_datatype(module, x2_obj, y2_obj, z2_obj, weights2_obj, &element_size2);
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
    PyObject *x1_array = PyArray_FromArray(x1_obj, NOTYPE_DESCR, requirements);
    PyObject *y1_array = PyArray_FromArray(y1_obj, NOTYPE_DESCR, requirements);
    PyObject *z1_array = PyArray_FromArray(z1_obj, NOTYPE_DESCR, requirements);
    PyObject *weights1_array = NULL;
    if(weights1_obj != NULL){
        weights1_array = PyArray_FromArray(weights1_obj, NOTYPE_DESCR, requirements);
    }

    PyObject *x2_array = NULL, *y2_array = NULL, *z2_array = NULL, *weights2_array = NULL;
    if(autocorr == 0) {
        x2_array = PyArray_FromArray(x2_obj, NOTYPE_DESCR, requirements);
        y2_array = PyArray_FromArray(y2_obj, NOTYPE_DESCR, requirements);
        z2_array = PyArray_FromArray(z2_obj, NOTYPE_DESCR, requirements);
        if(weights2_obj != NULL){
            weights2_array = PyArray_FromArray(weights2_obj, NOTYPE_DESCR, requirements);
        }
    }

    if (x1_array == NULL || y1_array == NULL || z1_array == NULL ||
        (autocorr == 0 && (x2_array == NULL || y2_array == NULL || z2_array == NULL))) {
        Py_XDECREF(x1_array);
        Py_XDECREF(y1_array);
        Py_XDECREF(z1_array);
        Py_XDECREF(weights1_array);

        Py_XDECREF(x2_array);
        Py_XDECREF(y2_array);
        Py_XDECREF(z2_array);
        Py_XDECREF(weights2_array);
        char msg[1024];
        snprintf(msg, 1024, "TypeError: In %s: Could not convert input to arrays of allowed floating point types (doubles or floats). Are you passing numpy arrays?",
                 __FUNCTION__);
        countpairs_mocks_error_out(module, msg);
        Py_RETURN_NONE;
    }

    /* Get pointers to the data as C-types. */
    void *phiD1   = PyArray_DATA((PyArrayObject *)x1_array);
    void *thetaD1 = PyArray_DATA((PyArrayObject *)y1_array);
    void *czD1    = PyArray_DATA((PyArrayObject *)z1_array);
    void *weights1=NULL;
    if(weights1_array != NULL){
        weights1 = PyArray_DATA((PyArrayObject *) weights1_array);
    }

    void *phiD2=NULL, *thetaD2=NULL, *czD2=NULL, *weights2=NULL;
    if(autocorr == 0) {
        phiD2   = PyArray_DATA((PyArrayObject *) x2_array);
        thetaD2 = PyArray_DATA((PyArrayObject *) y2_array);
        czD2    = PyArray_DATA((PyArrayObject *) z2_array);
        if(weights2_array != NULL){
            weights2 = PyArray_DATA((PyArrayObject *) weights2_array);
        }
    }
    options.float_type = element_size;

    /* Pack the weights into extra_options */
    for(int64_t w = 0; w < extra.weights0.num_weights; w++){
        extra.weights0.weights[w] = (char *) weights1 + w*ND1*element_size;
        if(autocorr == 0){
            extra.weights1.weights[w] = (char *) weights2 + w*ND2*element_size;
        }
    }

    NPY_BEGIN_THREADS_DEF;
    NPY_BEGIN_THREADS;

    results_countpairs_mocks_s_mu results;
    double c_api_time = 0.0;
    int status = countpairs_mocks_s_mu(ND1,phiD1,thetaD1,czD1,
                                       ND2,phiD2,thetaD2,czD2,
                                       nthreads,
                                       autocorr,
                                       binfile,
                                       mu_max,
                                       nmu_bins,
                                       cosmology,
                                       &results,
                                       &options,
                                       &extra);
    if(options.c_api_timer) {
        c_api_time = options.c_api_time;
    }
    NPY_END_THREADS;

    /* Clean up. */
    Py_DECREF(x1_array);Py_DECREF(y1_array);Py_DECREF(z1_array);Py_XDECREF(weights1_array);//x1 should absolutely not be NULL
    Py_XDECREF(x2_array);Py_XDECREF(y2_array);Py_XDECREF(z2_array);Py_XDECREF(weights2_array);//x2 might be NULL depending on value of autocorr

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
    double rlow=results.supp[0];
    const double dmu = mu_max/(double)results.nmu_bins ;

    for(int i=1;i<results.nsbin;i++) {
        for(int j=0;j<results.nmu_bins;j++) {
            const int bin_index = i*(results.nmu_bins + 1) + j;
            const double savg = results.savg[bin_index];
            const double weight_avg = results.weightavg[bin_index];
            PyObject *item = Py_BuildValue("(ddddkd)", rlow,results.supp[i],savg,(j+1)*dmu,results.npairs[bin_index], weight_avg);
            PyList_Append(ret, item);
            Py_XDECREF(item);
        }
        rlow=results.supp[i];
    }
    free_results_mocks_s_mu(&results);

    PyObject *rettuple = Py_BuildValue("(Od)", ret, c_api_time);
    Py_DECREF(ret);  // transfer reference ownership to the tuple
    return rettuple;
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

    PyArrayObject *x1_obj=NULL, *y1_obj=NULL, *weights1_obj=NULL;
    PyArrayObject *x2_obj=NULL, *y2_obj=NULL, *weights2_obj=NULL;
    int nthreads=1;
    char *binfile, *weighting_method_str = NULL;;
    int autocorr=0;
    struct config_options options = get_config_options();
    options.verbose=0;
    options.instruction_set=-1;
    options.link_in_dec=1;
    options.link_in_ra=1;
    options.fast_acos=0;
    options.c_api_timer=0;
    options.enable_min_sep_opt = 1;
    options.copy_particles = 1;
    int8_t ra_bin_ref=options.bin_refine_factors[0],
        dec_bin_ref=options.bin_refine_factors[1];
    static char *kwlist[] = {
        "autocorr",
        "nthreads",
        "binfile",
        "RA1",
        "DEC1",
        "weights1",
        "RA2",
        "DEC2",
        "weights2",
        "link_in_dec",
        "link_in_ra",
        "verbose", /* keyword verbose -> print extra info at runtime + progressbar */
        "output_thetaavg",
        "fast_acos",
        "ra_refine_factor",
        "dec_refine_factor",
        "max_cells_per_dim",
        "copy_particles",
        "enable_min_sep_opt",
        "c_api_timer",
        "isa",/* instruction set to use of type enum isa; valid values are AVX512F, AVX, SSE, FALLBACK */
        "weight_type",
        NULL
    };


    if ( ! PyArg_ParseTupleAndKeywords(args, kwargs, "iisO!O!|O!O!O!O!bbbbbbbhbbbis", kwlist,
                                       &autocorr,&nthreads,&binfile,
                                       &PyArray_Type,&x1_obj,
                                       &PyArray_Type,&y1_obj,
                                       &PyArray_Type,&weights1_obj,
                                       &PyArray_Type,&x2_obj,//optional parameters -> if autocorr == 1, not checked; required if autocorr=0
                                       &PyArray_Type,&y2_obj,
                                       &PyArray_Type,&weights2_obj,
                                       &(options.link_in_dec),
                                       &(options.link_in_ra),
                                       &(options.verbose),
                                       &(options.need_avg_sep),
                                       &(options.fast_acos),
                                       &ra_bin_ref, &dec_bin_ref,
                                       &(options.max_cells_per_dim),
                                       &(options.copy_particles),
                                       &(options.enable_min_sep_opt),
                                       &(options.c_api_timer),
                                       &(options.instruction_set),
                                       &weighting_method_str)

         ) {
        PyObject_Print(kwargs, stdout, 0);
        fprintf(stdout, "\n");

        char msg[1024];
        int len=snprintf(msg, 1024,"ArgumentError: In DDtheta_mocks> Could not parse the arguments. Input parameters are: \n");

        /* How many keywords do we have? Subtract 1 because of the last NULL */
        const size_t nitems = sizeof(kwlist)/sizeof(*kwlist) - 1;
        int status = print_kwlist_into_msg(msg, 1024, len, kwlist, nitems);
        if(status != EXIT_SUCCESS) {
            fprintf(stderr,"Error message does not contain all of the keywords\n");
        }
        countpairs_mocks_error_out(module,msg);
        Py_RETURN_NONE;
    }
    options.autocorr=autocorr;
    options.periodic=0;//doesn't matter but noting intent by setting it to 0

    /*This is for the fastest isa */
    if(options.instruction_set == -1) {
        options.instruction_set = highest_isa_mocks;
    }
    if(ra_bin_ref != options.bin_refine_factors[0] ||
       dec_bin_ref != options.bin_refine_factors[1]) {
        options.bin_refine_factors[0] = ra_bin_ref;
        options.bin_refine_factors[1] = dec_bin_ref;
        set_bin_refine_scheme(&options, BINNING_CUST);//custom binning -> code will honor requested binning scheme
    }

    /* Validate the user's choice of weighting method */
    weight_method_t weighting_method;
    int wstatus = get_weight_method_by_name(weighting_method_str, &weighting_method);
    if(wstatus != EXIT_SUCCESS){
        char msg[1024];
        snprintf(msg, 1024, "ValueError: In %s: unknown weight_type \"%s\"!", __FUNCTION__, weighting_method_str);
        countpairs_mocks_error_out(module, msg);
        Py_RETURN_NONE;
    }
    if(weighting_method == NONE){
        // Do not attempt to validate the weights array if it will not be used!
        weights1_obj = NULL;
        weights2_obj = NULL;
    }

    size_t element_size;
    /* We have numpy arrays and all the required inputs*/
    /* How many data points are there? And are they all of floating point type */
    const int64_t ND1 = check_dims_and_datatype_ra_dec(module, x1_obj, y1_obj, &element_size);
    if(ND1 == -1) {
        //Error has already been set -> simply return
        Py_RETURN_NONE;
    }

    int found_weights = weights1_obj == NULL ? 0 : PyArray_SHAPE(weights1_obj)[0];
    struct extra_options extra = get_extra_options(weighting_method);
    if(extra.weights0.num_weights > 0 && extra.weights0.num_weights != found_weights){
        char msg[1024];
        snprintf(msg, 1024, "ValueError: In %s: specified weighting method %s which requires %"PRId64" weight(s)-per-particle, but found %d weight(s) instead!\n",
                 __FUNCTION__, weighting_method_str, extra.weights0.num_weights, found_weights);
        countpairs_mocks_error_out(module, msg);
        Py_RETURN_NONE;
    }

    if(extra.weights0.num_weights > 0 && found_weights > MAX_NUM_WEIGHTS){
        char msg[1024];
        snprintf(msg, 1024, "ValueError: In %s: Provided %d weights-per-particle, but the code was compiled with MAX_NUM_WEIGHTS=%d.\n",
                 __FUNCTION__, found_weights, MAX_NUM_WEIGHTS);
        countpairs_mocks_error_out(module, msg);
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
        if((weights1_obj == NULL) != (weights2_obj == NULL)){
            snprintf(msg, 1024, "ValueError: In %s: If autocorr is 0, must pass either zero or two sets of weights.\n",
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
    PyObject *x1_array = PyArray_FromArray(x1_obj, NOTYPE_DESCR, requirements);
    PyObject *y1_array = PyArray_FromArray(y1_obj, NOTYPE_DESCR, requirements);
    PyObject *weights1_array = NULL;
    if(weights1_obj != NULL){
        weights1_array = PyArray_FromArray(weights1_obj, NOTYPE_DESCR, requirements);
    }

    PyObject *x2_array = NULL, *y2_array = NULL, *weights2_array = NULL;
    if(autocorr == 0) {
        x2_array = PyArray_FromArray(x2_obj, NOTYPE_DESCR, requirements);
        y2_array = PyArray_FromArray(y2_obj, NOTYPE_DESCR, requirements);
        if(weights2_obj != NULL){
            weights2_array = PyArray_FromArray(weights2_obj, NOTYPE_DESCR, requirements);
        }
    }

    if (x1_array == NULL || y1_array == NULL ||
        (autocorr == 0 && (x2_array == NULL || y2_array == NULL))) {
        Py_XDECREF(x1_array);
        Py_XDECREF(y1_array);
        Py_XDECREF(weights1_array);

        Py_XDECREF(x2_array);
        Py_XDECREF(y2_array);
        Py_XDECREF(weights2_array);
        char msg[1024];
        snprintf(msg, 1024, "TypeError: In %s: Could not convert input to arrays of allowed floating point types (doubles or floats). Are you passing numpy arrays?",
                 __FUNCTION__);
        countpairs_mocks_error_out(module, msg);
        Py_RETURN_NONE;
    }

    /* Get pointers to the data as C-types. */
    void *phiD1   = PyArray_DATA((PyArrayObject *) x1_array);
    void *thetaD1 = PyArray_DATA((PyArrayObject *) y1_array);
    void *weights1 = NULL;
    if(weights1_array != NULL){
        weights1 = PyArray_DATA((PyArrayObject *) weights1_array);
    }

    void *phiD2=NULL, *thetaD2=NULL, *weights2=NULL;
    if(autocorr == 0) {
        phiD2   = PyArray_DATA((PyArrayObject *) x2_array);
        thetaD2 = PyArray_DATA((PyArrayObject *) y2_array);
        if(weights2_array != NULL){
            weights2 = PyArray_DATA((PyArrayObject *) weights2_array);
        }
    }

    /* Pack the weights into extra_options */
    for(int64_t w = 0; w < extra.weights0.num_weights; w++){
        extra.weights0.weights[w] = (char *) weights1 + w*ND1*element_size;
        if(autocorr == 0){
            extra.weights1.weights[w] = (char *) weights2 + w*ND2*element_size;
        }
    }

    NPY_BEGIN_THREADS_DEF;
    NPY_BEGIN_THREADS;

    results_countpairs_theta results;
    options.float_type = element_size;
    double c_api_time=0.0;
    int status = countpairs_theta_mocks(ND1,phiD1,thetaD1,
                                        ND2,phiD2,thetaD2,
                                        nthreads,
                                        autocorr,
                                        binfile,
                                        &results,
                                        &options,
                                        &extra);
    if(options.c_api_timer) {
        c_api_time = options.c_api_time;
    }
    NPY_END_THREADS;


    /* Clean up. */
    Py_DECREF(x1_array);Py_DECREF(y1_array);Py_XDECREF(weights1_array);//x1/y1 (representing ra1,dec1) should not be NULL
    Py_XDECREF(x2_array);Py_XDECREF(y2_array);Py_XDECREF(weights2_array);//x2/y2 may be NULL (in case of autocorr)

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
        const double theta_avg = results.theta_avg[i];
        const double weight_avg = results.weightavg[i];
        PyObject *item = Py_BuildValue("(dddkd)", rlow,results.theta_upp[i],theta_avg,results.npairs[i], weight_avg);
        PyList_Append(ret, item);
        Py_XDECREF(item);
        rlow=results.theta_upp[i];
    }
    free_results_countpairs_theta(&results);

    PyObject *rettuple = Py_BuildValue("(Od)", ret, c_api_time);
    Py_DECREF(ret);  // transfer reference ownership to the tuple
    return rettuple;
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
    struct config_options options = get_config_options();
    options.is_comoving_dist = 0;
    options.verbose=0;
    options.instruction_set=-1;
    options.copy_particles=1;
    options.c_api_timer=0;

    /* Reset the bin refine factors default (since the VPF is symmetric in XYZ, conceptually the binning should be identical in all three directions)*/
    int bin_ref[] = {1,1,1};
    set_bin_refine_factors(&options, bin_ref);

    int8_t xbin_ref=options.bin_refine_factors[0],
        ybin_ref=options.bin_refine_factors[1],
        zbin_ref=options.bin_refine_factors[2];

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
        "is_comoving_dist",
        "verbose", /* keyword verbose -> print extra info at runtime + progressbar */
        "xbin_refine_factor",
        "ybin_refine_factor",
        "zbin_refine_factor",
        "max_cells_per_dim",
        "copy_particles",
        "c_api_timer",
        "isa",/* instruction set to use of type enum isa; valid values are AVX512F, AVX, SSE, FALLBACK */
        NULL
    };


    if ( ! PyArg_ParseTupleAndKeywords(args, kwargs, "diiiisiO!O!O!O!O!O!|bbbbbhbbi", kwlist,
                                       &rmax,&nbin,&num_spheres,&num_pN,&threshold_neighbors,&centers_file,&cosmology,
                                       &PyArray_Type,&x1_obj,
                                       &PyArray_Type,&y1_obj,
                                       &PyArray_Type,&z1_obj,
                                       &PyArray_Type,&x2_obj,
                                       &PyArray_Type,&y2_obj,
                                       &PyArray_Type,&z2_obj,
                                       &(options.is_comoving_dist),
                                       &(options.verbose),
                                       &xbin_ref, &ybin_ref, &zbin_ref,
                                       &(options.max_cells_per_dim),
                                       &(options.copy_particles),
                                       &(options.c_api_timer),
                                       &(options.instruction_set))

         ) {

        PyObject_Print(kwargs, stdout, 0);
        fprintf(stdout, "\n");

        char msg[1024];
        int len=snprintf(msg, 1024,"ArgumentError: In vpf_mocks> Could not parse the arguments. Input parameters are: \n");

        /* How many keywords do we have? Subtract 1 because of the last NULL */
        const size_t nitems = sizeof(kwlist)/sizeof(*kwlist) - 1;
        int status = print_kwlist_into_msg(msg, 1024, len, kwlist, nitems);
        if(status != EXIT_SUCCESS) {
            fprintf(stderr,"Error message does not contain all of the keywords\n");
        }
        countpairs_mocks_error_out(module,msg);

        Py_RETURN_NONE;
    }
    /*This is for the fastest isa */
    if(options.instruction_set == -1) {
        options.instruction_set = highest_isa_mocks;
    }

    if(xbin_ref != options.bin_refine_factors[0] ||
       ybin_ref != options.bin_refine_factors[1] ||
       zbin_ref != options.bin_refine_factors[2]) {
        options.bin_refine_factors[0] = xbin_ref;
        options.bin_refine_factors[1] = ybin_ref;
        options.bin_refine_factors[2] = zbin_ref;
        set_bin_refine_scheme(&options, BINNING_CUST);//custom binning -> code will honor requested binning scheme
    }

    size_t element_size;
    /* We have numpy arrays and all the required inputs*/
    /* How many data points are there? And are they all of floating point type */
    const int64_t ND1 = check_dims_and_datatype(module, x1_obj, y1_obj, z1_obj, NULL, &element_size);
    if(ND1 == -1) {
        //Error has already been set -> simply return
        Py_RETURN_NONE;
    }

    size_t element_size2;
    const int64_t ND2 = check_dims_and_datatype(module, x2_obj, y2_obj, z2_obj, NULL, &element_size2);
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
    PyObject *x1_array = PyArray_FromArray(x1_obj, NOTYPE_DESCR, requirements);
    PyObject *y1_array = PyArray_FromArray(y1_obj, NOTYPE_DESCR, requirements);
    PyObject *z1_array = PyArray_FromArray(z1_obj, NOTYPE_DESCR, requirements);
    PyObject *x2_array = PyArray_FromArray(x2_obj, NOTYPE_DESCR, requirements);
    PyObject *y2_array = PyArray_FromArray(y2_obj, NOTYPE_DESCR, requirements);
    PyObject *z2_array = PyArray_FromArray(z2_obj, NOTYPE_DESCR, requirements);

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
    void *phiD1   = PyArray_DATA((PyArrayObject *) x1_array);
    void *thetaD1 = PyArray_DATA((PyArrayObject *) y1_array);
    void *czD1    = PyArray_DATA((PyArrayObject *) z1_array);

    void *phiD2   = PyArray_DATA((PyArrayObject *) x2_array);
    void *thetaD2 = PyArray_DATA((PyArrayObject *) y2_array);
    void *czD2    = PyArray_DATA((PyArrayObject *) z2_array);

    NPY_BEGIN_THREADS_DEF;
    NPY_BEGIN_THREADS;

    results_countspheres_mocks results;
    options.float_type = element_size;
    double c_api_time = 0.0;
    int status = countspheres_mocks(ND1, phiD1,thetaD1, czD1,
                                    ND2, phiD2,thetaD2, czD2,
                                    threshold_neighbors,
                                    rmax, nbin, num_spheres,
                                    num_pN,
                                    centers_file,
                                    cosmology,
                                    &results,
                                    &options, NULL);
    if(options.c_api_timer) {
        c_api_time = options.c_api_time;
    }
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

    PyObject *rettuple = Py_BuildValue("(Od)", ret, c_api_time);
    Py_DECREF(ret);  // transfer reference ownership to the tuple
    return rettuple;
}
