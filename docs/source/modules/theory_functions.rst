.. _theory_functions:

Detailed API for Clustering Statistics on Simulations
=======================================================

All of these can be imported from :py:mod:`Corrfunc.theory`. See the complete reference here :py:mod:`Corrfunc`.

.. currentmodule:: Corrfunc.theory

Clustering in 3-D
------------------

* Pair counts for (auto or cross) correlations for :math:`\xi(r)` -- :py:mod:`Corrfunc.theory.DD`
* Auto-correlation on periodic, cosmological boxes, :math:`\xi(r)`, -- :py:mod:`Corrfunc.theory.xi`

Clustering in 2-D
------------------

* Pair counts (auto or cross) correlations for :math:`\xi(rp, \pi)` --  :py:mod:`Corrfunc.theory.DDrppi`
* Pair counts (auto or cross) correlations for :math:`\xi(s, \mu)` -- :py:mod:`Corrfunc.theory.DDsmu`     
* Projected auto-correlation function, :math:`wp(rp)` --  :py:mod:`Corrfunc.theory.wp`

Counts-in-cells
----------------

* Void Probability functions and counts-in-cells stats :math:`pN(r)` -- :py:mod:`Corrfunc.theory.vpf`
