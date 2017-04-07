.. _rr_autocorrelations:

Notes on the Random-Random Term in Autocorrelations
===================================================

The following discussion is adapted from `this notebook <http://nbviewer.jupyter.org/gist/lgarrison/1efabe4430429996733a9d29397423d2>`_ by `Lehman Garrison <https://lgarrison.github.io>`_.

.. math::
    \newcommand{\RR}{\mathrm{RR}}
    \newcommand{\DD}{\mathrm{DD}}
    
When computing a two-point correlation function estimator like

.. math::
    \xi(r) = \frac{\DD}{\RR} - 1,
    
the :math:`\RR` term can be computed analytically if the domain is a periodic box.
Often, this is done as

.. math::
    \begin{align}
    \RR_i &= N V_i \bar\rho \\
    &= N V_i \frac{N}{L^3}
    \end{align}
    
where :math:`\RR_i` is the expected number of random-random pairs in bin :math:`i`, :math:`N` is the total number of points, :math:`V_i` is the volume (or area if 2D) of bin :math:`i`, :math:`L` is the box size, and :math:`\bar\rho` is the average density in the box.

However, using :math:`\bar\rho = \frac{N}{L^3}` is only correct for continuous fields, not sets of particles.  When sitting on a particle, only :math:`N-1` particles are available to be in a bin at some non-zero distance.  The remaining particle is the particle you're sitting on, which is always at distance :math:`0`.  Thus, the correct expression is

.. math::
    \RR_i = N V_i \frac{N-1}{L^3}.

See `this notebook <http://nbviewer.jupyter.org/gist/lgarrison/1efabe4430429996733a9d29397423d2>`_ for an empirical demonstration of this effect; specifically, that computing the density with :math:`N-1` is correct, and that using :math:`N` introduces bias of order :math:`\frac{1}{N}` into the estimator.  This is a tiny correction for large :math:`N` problems, but important for small :math:`N`.

Any ``Corrfunc`` function that returns a clustering statistic (not just raw pair counts) implements this correction. 
Currently, this includes :py:mod:`Corrfunc.theory.xi` and :py:mod:`Corrfunc.theory.wp`.

Cross-correlations of two different particle sets don't suffer from this problem; the particle you're sitting on is never part of the set of particles under consideration for pair-making.

``Corrfunc`` also allows bins of zero separation, in which "self-pairs" are included in the pair counting.  :math:`\RR_i` must reflect this by simply adding :math:`N` to any such bin.

.. _weighted_rr:

RR in Weighted Clustering Statistics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We can extend the above discussion to weighted correlation functions in which
each particle is assigned a weight, and the pair weight is taken as the product
of the particle weights (see :ref:`weighted_correlations`).

Let :math:`w_j` be the weight of particle :math:`j`, and :math:`W` be the sum of the weights.
We will define the "unclustered" particle distribution to be the case of :math:`N` particles
uniformly distributed, where each is assigned the mean weight :math:`\bar w`.  We thus have

.. math::
    \begin{align}
    \RR_i &= \sum_{j=1}^N \bar w (W - \bar w) \frac{V_i}{L^3} \\
    &= (W^2 - \bar w W) \frac{V_i}{L^3} \\
    &= W^2\left(1 - \frac{1}{N}\right) \frac{V_i}{L^3}.
    \end{align}

When the particles all have :math:`w_j = 1`, then :math:`W = N` and we recover the unweighted result from above.

There are other ways to define the unclustered distribution.  If we were to redistribute
the particles uniformly but preserve their individual weights, we would find

.. math::
    \begin{align}
    \RR_i &= \sum_{j=1}^N w_j (W - w_j) \frac{V_i}{L^3} \\
    &= \left(W^2 - \sum_{j=1}^N w_j^2\right) \frac{V_i}{L^3}.
    \end{align}

This is not what we use in ``Corrfunc``, but this should help illuminate some of the considerations that
go into defining the "unclustered" case when writing a custom weight function (see :ref:`custom_weighting`).
