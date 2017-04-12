.. _staticlibrary-interface:

***********************************************
Using the static library interface in Corrfunc
***********************************************

This guide assumes that you already followed the :ref:`step_by_step_install`
section of the documentation to get the package and its dependencies set
up on your machine. This guide also assumes some familiarity with C coding.

This concepts in this guide are implemented in the files
``theory/examples/run_correlations.c`` and
``mocks/examples/run_correlations_mocks.c`` for simulations and mock
catalogs respectively.

The basic principle of using the static libraries has the following steps:

* Include the appropriate header to get the correct function signature (at
  compile time)
* In your code, include call with clustering function with appropriate parameters
* Compile your code with ``-I </path/to/Corrfunc/include>`` flags. If you have
  installed Corrfunc via ``pip``, then use
  ``os.path.join(os.path.dirname(Corrfunc.__file__), ../include/)`` as the
  include header.
* Link your code with the appropriate static library. Look in the
  ``examples/Makefile`` for the linker flags.
* Run your code


Worked out example C code for clustering statistics in simulation boxes
========================================================================

Common setup code for the simulation C routines
------------------------------------------------

In this code section, we will setup the arrays and the overall common inputs
required by the C static libraries. 

.. code-block:: c

          #include "io.h"
          
          const char file[] = {"theory/tests/data/gals_Mr19.ff"}; 
          const char fileformat[] = {"f"};  
          const char binfile[] = {"theory/tests/bins"};
          const double boxsize=420.0;
          const double pimax=40.0;
          int autocorr=1;
          const int nthreads=2;

          double *x1=NULL, *y1=NULL, *z1=NULL, *x2=NULL, *y2=NULL, *z2=NULL;

          const int64_t ND1 = read_positions(file,fileformat,sizeof(*x1),3, &x1, &y1, &z1);
          x2 = x1;
          y2 = y1;
          z2 = z1;
          const int64_t ND2 = ND1;

          struct config_options options = get_config_options();
          options.verbose = 1;        
          options.need_avg_sep = 1;   
          options.periodic = 1;       
          options.float_type = sizeof(*x1); 


Calculating 2-D projected auto-correlation (``theory/wp/libcountpairs_wp.a``)
--------------------------------------------------------------------------------

Corrfunc can directly compute the projected auto-correlation function,
:math:`w_p(r_p)`. This calculation sets periodic boundary conditions. Randoms
are calculated analytically based on the supplied boxsize. The projected
separation, :math:`r_p` is calculated in the X-Y plane while the line-of-sight
separation, :math:`\pi` is calculated in the Z plane. Only pairs with
:math:`\pi` separation less than :math:`\pi_{max}` are counted.

.. code-block:: c

          #include "countpairs_wp.h"
          
          results_countpairs_wp results;
          int status = countpairs_wp(ND1,x1,y1,z1,
                                     boxsize,
                                     nthreads,
                                     binfile,
                                     pimax,
                                     &results,
                                     &options, NULL);
                                     
          if(status != EXIT_SUCCESS) {
              fprintf(stderr,"Runtime error occurred while using wp static library\n");
              return status;
          }
          
          double rlow=results.rupp[0];
          for(int i=1;i<results.nbin;++i) {
              fprintf(stdout,"%e\t%e\t%e\t%e\t%12"PRIu64" \n",
                             results.wp[i],results.rpavg[i],rlow,results.rupp[i],results.npairs[i]);
              rlow=results.rupp[i];
         }

This is the generic pattern for using all of the correlation function. Look in
``theory/examples/run_correlations.c`` for details on how to use all of the available
static libraries.
          
Worked out example C code for clustering statistics in mock catalogs
======================================================================
Corrfunc can calculate pair counts for mock catalogs. The input positions are
expected to be ``Right Ascension``, ``Declination`` and ``CZ`` (speed of light
times redshift, in ``Mpc/h``). Cosmology has to be specified since ``CZ`` needs
to be converted into co-moving distance. If you want to calculate in arbitrary
cosmology, then you have two options:

* convert ``CZ`` into co-moving distance, and then pass the converted array while setting ``config_option.is_comoving_dist=1``.
* Add another cosmology in ``utils/cosmology_params.c`` in the function
  ``init_cosmology``. Then, recompile the ``Corrfunc.mocks`` and pass
  ``cosmology=integer_for_newcosmology`` into the relevant functions.


Common setup code for the mocks C routines
--------------------------------------------
In this code section, we will setup the arrays and the overall common inputs
required by the C static libraries. 

.. code-block:: c

   #include "io.h"   //for read_positions function
          
   const char file[] = {"mocks/tests/data/Mr19_mock_northonly.rdcz.dat"};
   const char fileformat[] = {"a"};     // ascii format
   const char binfile[] = {"mocks/tests/bins"};
   const double pimax=40.0;
   int autocorr=1;
   const int nthreads=2;
   const int cosmology=1;  // 1->LasDamas cosmology, 2->Planck

   // This computes in double-precision. Change to float for computing in float
   double *ra1=NULL, *dec1=NULL, *cz1=NULL, *ra2=NULL, *dec2=NULL, *cz2=NULL;

   //Read-in the data
   const int64_t ND1 = read_positions(file,fileformat,sizeof(*ra1),3, &ra1, &dec1, &cz1);

   ra2 = ra1;
   dec2 = dec1;
   cz2 = cz1;
   const int64_t ND2 = ND1;

   struct config_options options = get_config_options();
   options.verbose=1;
   options.periodic=0;
   options.need_avg_sep=1;
   options.float_type = sizeof(*ra1);

   
Calculating 2-D pair counts (``mocks/DDrppi_mocks/libcountpairs_rp_pi_mocks.a``)
-----------------------------------------------------------------------------------
Here is a code snippet demonstrating how to calculate :math:`DD(r_p, \pi)` for
mock catalogs. The projected separation, :math:`r_p` and line of sight
separation, :math:`\pi` are calculated using the following equations from `Zehavi et
al 2002 <http://adsabs.harvard.edu/abs/2002ApJ...571..172Z>`_:

.. math::
   
   \mathbf{s} &= \mathbf{v_1} - \mathbf{v_2}, \\
   \mathbf{l} &= \frac{1}{2}\left(\mathbf{v_1} + \mathbf{v_2}\right), \\
   \pi &= \left(\mathbf{s} \cdot \mathbf{l}\right)/\Vert\mathbf{l}\Vert, \\
   r_p^2 &= \mathbf{s} \cdot \mathbf{s} - \pi^2

where, :math:`\mathbf{v_1}` and :math:`\mathbf{v_2}` are the vectors for the
two points under consideration. Here is the C code for calling ``DDrppi_mocks``:

.. code-block:: c

                #include "countpairs_rp_pi_mocks.h"

                results_countpairs_mocks results;
                int status = countpairs_mocks(ND1,ra1,dec1,cz1,
                                              ND2,ra2,dec2,cz2,
                                              nthreads,
                                              autocorr,
                                              binfile,
                                              pimax,
                                              cosmology,
                                              &results,
                                              &options, NULL);

                const double dpi = pimax/(double)results.npibin ;
                const int npibin = results.npibin;
                for(int i=1;i<results.nbin;i++) {
                    const double logrp = LOG10(results.rupp[i]);
                    for(int j=0;j<npibin;j++) {
                        int index = i*(npibin+1) + j;
                        fprintf(stdout,"%10"PRIu64" %20.8lf %20.8lf  %20.8lf \n",results.npairs[index],results.rpavg[index],logrp,(j+1)*dpi);
                    }
                }

This is the generic pattern for using all of the correlation function. Look in
``mocks/examples/run_correlations_mocks.c`` for details on how to use all of the available
static libraries.
