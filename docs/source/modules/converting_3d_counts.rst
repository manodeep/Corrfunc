.. _converting_3d_counts:

Converting 3D pair counts into a correlation function
======================================================

3D pair counts can be converted into a correlation function
by using the helper function :py:mod:`Corrfunc.utils.convert_3d_counts_to_cf`.
First, we have to compute the relevant pair counts using the python
wrapper :py:mod:`Corrfunc.theory.DD`

.. code-block:: python

          >>> import numpy as np
          >>> from os.path import dirname, abspath, join as pjoin
          >>> from Corrfunc.theory.DD import DD
          >>> from Corrfunc.io import read_catalog
          >>> from Corrfunc.utils import convert_3d_counts_to_cf

          >>> # Read the supplied galaxies on a periodic box
          >>> X, Y, Z = read_catalog()
          >>> N = len(X)
          >>> boxsize = 420.0
          >>> nthreads = 2

          # Generate randoms on the box
          >>> rand_N = 3*N
          >>> rand_X = np.random.uniform(0, boxsize, rand_N)
          >>> rand_Y = np.random.uniform(0, boxsize, rand_N)
          >>> rand_Z = np.random.uniform(0, boxsize, rand_N)

          # Setup the bins
          >>> nbins = 10
          >>> bins = np.linspace(0.1, 10.0, nbins + 1) # note that +1 to nbins
              
          # Auto pair counts in DD
          >>> autocorr=1
          >>> DD_counts = DD(autocorr, nthreads, bins, X, Y, Z,
          ...               periodic=False, verbose=True)

          # Cross pair counts in DR
          >>> autocorr=0
          >>> DR_counts = DD(autocorr, nthreads, bins, X, Y, Z,
          ...               X2=rand_X, Y2=rand_Y, Z2=rand_Z,
          ...               periodic=False, verbose=True)
                         
          # Auto pairs counts in RR
          >>> autocorr=1                         
          >>> RR_counts = DD(autocorr, nthreads, bins, rand_X, rand_Y, rand_Z,
          ...                periodic=False, verbose=True)

          # All the pair counts are done, get the correlation function
          >>> cf = convert_3d_counts_to_cf(N, N, rand_N, rand_N,
          ...                             DD_counts, DR_counts,
          ...                             DR_counts, RR_counts)
          
See the complete reference here :py:mod:`Corrfunc`.
   
                   
