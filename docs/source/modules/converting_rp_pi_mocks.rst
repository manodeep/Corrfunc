.. _converting_rp_pi_mocks:

Calculating the projected correlation function, :math:`wp(rp)`
==============================================================

2-D Pair counts can be converted into a :math:`wp(rp)`
by using the helper function :py:mod:`Corrfunc.utils.convert_rp_pi_counts_to_wp`.
First, we have to compute the relevant pair counts using the python
wrapper :py:mod:`Corrfunc.mocks.DDrppi_mocks`

.. code-block:: python

          >>> import numpy as np
          >>> from os.path import dirname, abspath, join as pjoin          
          >>> import Corrfunc
          >>> from Corrfunc.mocks.DDrppi_mocks import DDrppi_mocks
          >>> from Corrfunc.io import read_catalog
          >>> from Corrfunc.utils import convert_rp_pi_counts_to_wp

          >>> galaxy_catalog=pjoin(dirname(abspath(Corrfunc.__file__)),
          ...                      "../mocks/tests/data", "Mr19_mock_northonly.rdcz.ff")

          # Read the supplied galaxies on a periodic box
          >>> RA, DEC, CZ = read_catalog(galaxy_catalog)
          >>> N = len(RA)

          # Read the supplied randoms catalog
          >>> random_catalog=pjoin(dirname(abspath(Corrfunc.__file__)),
          ...                      "../mocks/tests/data", "Mr19_randoms_northonly.rdcz.ff")
          >>> rand_RA, rand_DEC, rand_CZ = read_catalog(random_catalog)
          >>> rand_N = len(rand_RA)
          
          # Setup the bins
          >>> nbins = 10
          >>> bins = np.linspace(0.1, 20.0, nbins + 1)
          >>> pimax = 40.0

          >>> cosmology = 1
          >>> nthreads = 2

          # Auto pair counts in DD
          >>> autocorr=1
          >>> DD_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins,
          ...                          RA, DEC, CZ)


          # Cross pair counts in DR
          >>> autocorr=0
          >>> DR_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins,
          ...                          RA, DEC, CZ, 
          ...                          RA2=rand_RA, DEC2=rand_DEC, CZ2=rand_CZ)

                         
          # Auto pairs counts in RR
          >>> autocorr=1                         
          >>> RR_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins,
          ...                          rand_RA, rand_DEC, rand_CZ)

          # All the pair counts are done, get the angular correlation function
          >>> wp = convert_rp_pi_counts_to_wp(N, N, rand_N, rand_N,
          ...                                 DD_counts, DR_counts,
          ...                                 DR_counts, RR_counts, nbins, pimax)

See the complete reference here :py:mod:`Corrfunc`.
   
                   
