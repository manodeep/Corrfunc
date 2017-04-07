.. _computing_wp_and_xi:

Directly Computing :math:`\xi(r)` and :math:`wp(rp)`
====================================================

For a periodic cosmological box, the 3-d auto correlation, :math:`\xi(r)`, and
the projected auto correlation function, :math:`wp(rp)`, can be directly computed
using the Natural Estimator. The relevant python wrappers are present in
:py:mod:`Corrfunc.theory.xi` and :py:mod:`Corrfunc.theory.wp`.  See :ref:`rr_autocorrelations`
for details on how the Natural Estimator is computed.

.. code-block:: python

          >>> import numpy as np
          >>> from Corrfunc.theory.wp import wp
          >>> from Corrfunc.theory.xi import xi
          >>> from Corrfunc.io import read_catalog
          >>> X, Y, Z = read_catalog()
          >>> boxsize = 420.0
          >>> nthreads = 2
          >>> pimax = 40.0
          >>> nbins = 10
          >>> bins = np.linspace(0.1, 10.0, nbins + 1) # Note the + 1 to nbins
          >>> wp_counts = wp(boxsize, pimax, nthreads, bins, X, Y, Z)
          >>> xi_counts = xi(boxsize, nthreads, bins, X, Y, Z)
                

See the complete reference here :py:mod:`Corrfunc`.
