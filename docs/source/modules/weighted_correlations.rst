.. _weighted_correlations:

Computing Weighted Correlation Functions
========================================

Every clustering statistic in ``Corrfunc`` accepts an array
of weights that can be used to compute weighted correlation
functions. The API reference for each clustering statistic
(:py:mod:`Corrfunc.theory.xi`, :py:mod:`Corrfunc.mocks.DDrppi_mocks`,
etc.) contains examples of how to do this.  The interface is standard across functions: the
inputs are a ``weights`` array and a ``weight_type`` string
that specifies how to use the "point weights" to compute a "pair weight".
Currently, the only supported ``weight_type`` is ``pair_product``,
in which the pair weight is the product of the point weights
(but see :ref:`custom_weighting` for how to write your own
function).

.. warning::
    The computation of the weighted result is susceptible to loss of floating
    point precision, especially in single precision.  If you are using single
    precision, make sure you test double precision as well (by casting all
    pos and weight input arrays to type ``np.float64``, for example)
    and check that the difference with the single-precision result
    is acceptable.

If ``weight_type`` and ``weights`` (or ``weights1`` and ``weights2``
for cross-correlations) are given, the mean pair weight in a
separation bin will be given in the ``weightavg`` field of the
output.  This field is 0.0 if weights are disabled.

Pair counts (i.e. the ``npairs`` field in the ``results`` array)
are never affected by weights.  For theory functions like
:py:mod:`Corrfunc.theory.xi` and :py:mod:`Corrfunc.theory.wp`
that actually return a clustering statistic, the statistic is weighted.
For ``pair_product``, the distribution used to compute the
expected bin weight from an unclustered particle set (the ``RR`` term)
is taken to be a spatially uniform particle set where every particle
has the mean weight.  See :ref:`weighted_rr` for more discussion.

Running with weights incurrs a modest performance hit (around
20%, similar to enabling ``ravg``).  Weights are supported for
all instruction sets (SSE, AVX, and fallback).

Consider the following simple example adapted from the :py:mod:`Corrfunc.theory.xi`
docstring, in which we assign a weight of 0.5 to every particle and get
the expected average pair weight of 0.25 (last column of the output).
Note that ``xi`` (fourth column) is also weighted, but the case of uniform
weights is equivalent to the unweighted case.

::

    >>> from __future__ import print_function
    >>> import numpy as np
    >>> from os.path import dirname, abspath, join as pjoin
    >>> import Corrfunc
    >>> from Corrfunc.theory.xi import xi
    >>> binfile = pjoin(dirname(abspath(Corrfunc.__file__)),
    ...                 "../theory/tests/", "bins")
    >>> N = 100000
    >>> boxsize = 420.0
    >>> nthreads = 4
    >>> seed = 42
    >>> np.random.seed(seed)
    >>> X = np.random.uniform(0, boxsize, N)
    >>> Y = np.random.uniform(0, boxsize, N)
    >>> Z = np.random.uniform(0, boxsize, N)
    >>> weights = np.full_like(X, 0.5)
    >>> results = xi(boxsize, nthreads, binfile, X, Y, Z, weights=weights, weight_type='pair_product', output_ravg=True)
    >>> for r in results: print("{0:10.6f} {1:10.6f} {2:10.6f} {3:10.6f} {4:10d} {5:10.6f}"
    ...                         .format(r['rmin'], r['rmax'],
    ...                         r['ravg'], r['xi'], r['npairs'], r['weightavg']))
    ...                   # doctest: +NORMALIZE_WHITESPACE
          0.167536   0.238755   0.226592  -0.205733          4   0.250000
          0.238755   0.340251   0.289277  -0.176729         12   0.250000
          0.340251   0.484892   0.426819  -0.051829         40   0.250000
          0.484892   0.691021   0.596187  -0.131853        106   0.250000
          0.691021   0.984777   0.850100  -0.049207        336   0.250000
          0.984777   1.403410   1.225112   0.028543       1052   0.250000
          1.403410   2.000000   1.737153   0.011403       2994   0.250000
          2.000000   2.850200   2.474588   0.005405       8614   0.250000
          2.850200   4.061840   3.532018  -0.014098      24448   0.250000
          4.061840   5.788530   5.022241  -0.010784      70996   0.250000
          5.788530   8.249250   7.160648  -0.001588     207392   0.250000
          8.249250  11.756000  10.207213  -0.000323     601002   0.250000
         11.756000  16.753600  14.541171   0.000007    1740084   0.250000
         16.753600  23.875500  20.728773  -0.001595    5028058   0.250000
