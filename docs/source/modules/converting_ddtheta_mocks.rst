.. _converting_ddtheta_mocks:

Calculating the angular correlation function, :math:`\omega(\theta)`
====================================================================

Angular pair counts can be converted into a :math:`\omega(\theta)`
by using the helper function `Corrfunc.utils.convert_3d_counts_to_cf`.
First, we have to compute the relevant pair counts using the python
wrapper `Corrfunc.mocks.DDtheta_mocks`

.. code:: python

          import Corrfunc
          from Corrfunc.mocks import DDtheta_mocks
          from Corrfunc.io import read_catalog
          from Corrfunc.utils import convert_3d_counts_to_cf

          galaxy_catalog=pjoin(dirname(abspath(Corrfunc.__file__)),
                               "../mocks/tests/data", "Mr19_mock_northonly.rdcz.ff")
          # Read the supplied galaxies on a periodic box
          RA, DEC, _ = read_catalog(galaxy_catalog)

          # Read the supplied randoms catalog
          random_catalog=pjoin(dirname(abspath(Corrfunc.__file__)),
                               "../mocks/tests/data", "Mr19_randoms_northonly.rdcz.ff")
          rand_RA, rand_DEC, _ = read_catalog(random_catalog)
          
          # Setup the bins
          bins = np.linspace(0.1, 10.0, 10)

          # Auto pair counts in DD
          autocorr=1
          DD_counts = DDtheta_mocks(autocorr, nthreads, bins,
                                    RA, DEC,
                                    verbose=True)

          # Cross pair counts in DR
          autocorr=0
          DR_counts = DDtheta_mocks(autocorr, nthreads, bins,
                                    RA, DEC,
                                    rand_RA, rand_DEC,
                                    verbose=True)
                         
          # Auto pairs counts in RR
          autocorr=1                         
          RR_counts = DDtheta_mocks(autocorr, nthreads, bins,
                                    rand_RA, rand_DEC, 
                                    verbose=True)

          # All the pair counts are done, get the angular correlation function
          wtheta = convert_3d_counts_to_cf(N, N, rand_N, rand_N,
                                           DD_counts, DR_counts,
                                           DR_counts, RR_counts)
          
See the complete reference here `api/modules`.
   
                   
