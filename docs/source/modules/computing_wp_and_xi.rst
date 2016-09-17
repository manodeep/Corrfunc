.. _computing_wp_and_xi:

Directly Computing :math:`\xi(r)` and :math:`wp(rp)`
====================================================

For a periodic cosmological box, the 3-d auto correlation, :math:`\xi(r)`, and
the projected auto correlation function, :math:`wp(rp)`, can be directly computed
using the Natural Estimator. The relevant python wrappers are present in
`Corrfunc.theory.xi` and `Corrfunc.theory.wp`

.. code:: python

          from Corrfunc.theory import wp, xi
          from Corrfunc.io import read_catalog
          X, Y, Z = read_catalog()
          boxsize = 420.0
          nthreads = 2
          pimax = 40.0
          bins = np.linspace(0.1, 10.0, 10)
          wp_counts = wp(boxsize, nthreads, pimax, bins, X, Y, Z,
                          verbose=True)
          xi_counts = xi(boxsize, nthreads, bins, X, Y, Z,
                         verbose=True)
          
See the complete reference here `api/modules`.
   
                   
