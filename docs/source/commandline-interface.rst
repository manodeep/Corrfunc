.. _commandline-interface:

*********************************************
Using the command-line interface in Corrfunc
*********************************************

This guide assumes that you already followed the :ref:`step_by_step_install`
section of the documentation to get the package and its dependencies set
up on your machine. 

Calculating spatial clustering statistics in simulation boxes
==============================================================

Corrfunc can compute a range of spatial correlation functions and the
counts-in-cells. The easiest way to get help on the command-line is by calling
the executables without any input parameters. Here is the list of executables
associated with each type of clustering statistic:

======================  ========================
Clustering Statistic    Full path to executable
======================  ========================
:math:`DD(r)`            ``xi_theory/xi_of_r/DD``
:math:`DD(r_p,\pi)`      ``xi_theory/xi_rp_pi/DDrppi``
:math:`w_p(r_p)`         ``xi_theory/wp/wp``
:math:`\xi(r)`           ``xi_theory/xi/xi``
:math:`pN(n)`            ``xi_theory/vpf/vpf``
======================  ========================      
      

Calculating clustering statistics in mock catalogs
===================================================
The list of clustering statistics supported on mock catalogs and the associated
command-line executables are:

======================  ========================
Clustering Statistic    Full path to executable
======================  ========================
:math:`DD(r_p,\pi)`      ``xi_mocks/DDrppi/DDrppi_mocks``
:math:`DD(\theta)`       ``xi_mocks/wtheta/DDtheta_mocks``
:math:`pN(n)`            ``xi_mocks/vpf/vpf_mocks``
======================  ========================      


