.. _computing_wp_and_xi:

Directly Computing :math:`\xi(r)` and :math:`wp(rp)`
====================================================

For a periodic cosmological box, the 3-d auto correlation, :math:`\xi(r)`, and
the projected auto correlation function, :math:`wp(rp)`, can be directly computed
using the Natural Estimator. The relevant python wrappers are present in
`Corrfunc.theory.xi` and `Corrfunc.theory.wp`


          >>> import numpy as np
          >>> from Corrfunc.theory.wp import wp
          >>> from Corrfunc.theory.xi import xi
          >>> from Corrfunc.io import read_catalog
          >>> X, Y, Z = read_catalog()
          >>> boxsize = 420.0
          >>> nthreads = 2
          >>> pimax = 40.0
          >>> bins = np.linspace(0.1, 10.0, 10)
          >>> wp_counts = wp(boxsize, pimax, nthreads, bins, X, Y, Z)
          >>> print(wp_counts) # doctest: +NORMALIZE_WHITESPACE
          [(0.1, 1.2, 0.0, 104.5339046378777, 17091638L)
          (1.2, 2.3, 0.0, 42.07376597707123, 30440694L)
          (2.3, 3.4, 0.0, 28.920472781084925, 44233218L)
          (3.4, 4.5, 0.0, 23.058255725998702, 58006150L)
          (4.5, 5.6, 0.0, 19.386905571613653, 71517892L)
          (5.6, 6.7, 0.0, 16.61351453067251, 84665632L)
          (6.7, 7.8, 0.0, 14.450311699089156, 97574326L)
          (7.8, 8.9, 0.0, 12.826732246035473, 110446942L)
          (8.9, 10.0, 0.0, 11.521964904756938, 123239890L)]
          >>> xi_counts = xi(boxsize, nthreads, bins, X, Y, Z)
          >>> print(xi_counts) # doctest: +NORMALIZE_WHITESPACE
          [(0.1, 1.2, 0.0, 39.28812933254489, 6008688L)
          (1.2, 2.3, 0.0, 6.3336689308781455, 6611354L)
          (2.3, 3.4, 0.0, 2.5742902672606, 8376486L)
          (3.4, 4.5, 0.0, 1.5565253776430241, 11441056L)
          (4.5, 5.6, 0.0, 1.085519230008439, 15217204L)
          (5.6, 6.7, 0.0, 0.807223184600961, 19531808L)
          (6.7, 7.8, 0.0, 0.6290747078679189, 24449698L)
          (7.8, 8.9, 0.0, 0.5067697952809056, 29982762L)
          (8.9, 10.0, 0.0, 0.4206338629869881, 36195954L)]
                

See the complete reference here :py:mod:`Corrfunc`.

   
                   
