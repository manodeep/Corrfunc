.. _converting_rp_pi_mocks:

Calculating the projected correlation function, :math:`wp(rp)`
==============================================================

Angular pair counts can be converted into a :math:`wp(rp)`
by using the helper function `Corrfunc.utils.convert_rp_pi_counts_to_wp`.
First, we have to compute the relevant pair counts using the python
wrapper `Corrfunc.mocks.DDrppi_mocks`

.. code:: python

          import Corrfunc
          from Corrfunc.mocks import DDrppi_mocks
          from Corrfunc.io import read_catalog
          from Corrfunc.utils import convert_rp_pi_counts_to_wp

          galaxy_catalog=pjoin(dirname(abspath(Corrfunc.__file__)),
                               "../mocks/tests/data", "Mr19_mock_northonly.rdcz.ff")
          # Read the supplied galaxies on a periodic box
          RA, DEC, CZ = read_catalog(galaxy_catalog)

          # Read the supplied randoms catalog
          random_catalog=pjoin(dirname(abspath(Corrfunc.__file__)),
                               "../mocks/tests/data", "Mr19_randoms_northonly.rdcz.ff")
          rand_RA, rand_DEC, rand_CZ = read_catalog(random_catalog)
          
          # Setup the bins
          bins = np.linspace(0.1, 20.0, 10)
          pimax = 40.0

          # Auto pair counts in DD
          autocorr=1
          DD_counts = DDrppi_mocks(autocorr, nthreads, pimax, bins,
                                    RA, DEC, CZ, 
                                    verbose=True)

          # Cross pair counts in DR
          autocorr=0
          DR_counts = DDrppi_mocks(autocorr, nthreads, pimax, bins,
                                    RA, DEC, CZ, 
                                    RA2=rand_RA, DEC2=rand_DEC, CZ2=rand_CZ, 
                                    verbose=True)
                         
          # Auto pairs counts in RR
          autocorr=1                         
          RR_counts = DDrppi_mocks(autocorr, nthreads, pimax, bins,
                                    rand_RA, rand_DEC, rand_CZ,
                                    verbose=True)

          # All the pair counts are done, get the angular correlation function
          wp = convert_rp_pi_counts_to_wp(N, N, rand_N, rand_N,
                                          DD_counts, DR_counts,
                                          DR_counts, RR_counts)

See the complete reference here :py:mod:`Corrfunc`.
   
                   
