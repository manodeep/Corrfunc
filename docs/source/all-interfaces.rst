.. _all-interfaces:

*****************************************************
Cheat-sheet for all available interfaces in Corrfunc
*****************************************************

This guide assumes that you already followed the :ref:`step_by_step_install`
section of the documentation to get the package and its dependencies set
up on your machine. There are three available interfaces in Corrfunc

- :ref:`python-interface`
- :ref:`staticlibrary-interface`. The static libraries
  have the form ``libcount<statistic>.a``; the corresponding header file is named
  ``count<statistic>.h``.
- :ref:`commandline-interface`



Calculating spatial clustering statistics in simulation boxes
==============================================================

Corrfunc can compute a range of spatial correlation functions and the
counts-in-cells. The easiest way to get help on the command-line is by calling
the executables without any input parameters. Here is the list of executables
associated with each type of clustering statistic:

======================    ================================  ========================================  ====================================
Clustering Statistic      Python Interface                  Static library                            Command-line  (executable name)
======================    ================================  ========================================  ====================================
:math:`\xi(r)`            :py:mod:`Corrfunc.theory.DD`       ``theory/DD/libcountpairs.a``            ``theory/DD/DD``             
:math:`\xi(r_p,\pi)`      :py:mod:`Corrfunc.theory.DDrppi`   ``theory/DDrppi/libcountpairs_rp_pi.a``  ``theory/DDrppi/DDrppi``        
:math:`\xi(s,\mu)`        :py:mod:`Corrfunc.theory.DDsmu`    ``theory/DDsmu/libcountpairs_s_mu.a``    ``theory/DDsmu/DDsmu``        
:math:`w_p(r_p)`          :py:mod:`Corrfunc.theory.wp`       ``theory/wp/libcountpairs_wp.a``         ``theory/wp/wp``         
:math:`\xi(r)`            :py:mod:`Corrfunc.theory.xi`       ``theory/xi/libcountpairs_xi.a``         ``theory/xi/xi``         
:math:`pN(n)`             :py:mod:`Corrfunc.theory.vpf`      ``theory/vpf/libcountspheres.a``         ``theory/vpf/vpf``       
======================    ================================  ========================================  ====================================
      

Calculating clustering statistics in mock catalogs
===================================================
The list of clustering statistics supported on mock catalogs and the associated
command-line executables are:

======================   ======================================  =====================================================    =====================================
Clustering Statistic     Python Interface                        Static library                                           Command-line (executable name)
======================   ======================================  =====================================================    =====================================
:math:`\xi(r_p,\pi)`     :py:mod:`Corrfunc.mocks.DDrppi_mocks`    ``mocks/DDrppi_mocks/libcountpairs_rp_pi_mocks.a``      ``mocks/DDrppi_mocks/DDrppi_mocks``  
:math:`\xi(s,\mu)`       :py:mod:`Corrfunc.mocks.DDsmu_mocks`     ``mocks/DDsmu_mocks/libcountpairs_s_mu_mocks.a``        ``mocks/DDsmu_mocks/DDsmu_mocks``  
:math:`\omega(\theta)`   :py:mod:`Corrfunc.mocks.DDtheta_mocks`   ``mocks/DDtheta_mocks/libcountpairs_theta_mocks.a``     ``mocks/DDtheta_mocks/DDtheta_mocks``
:math:`pN(n)`            :py:mod:`Corrfunc.mocks.vpf_mocks`       ``mocks/vpf_mocks/libcountspheres_mocks``               ``mocks/vpf_mocks/vpf_mocks``        
======================   ======================================  =====================================================    =====================================
                                                                  
                                                                 
If you are not sure which correlation function to use, then please also see :ref:`which_corrfunc`.
