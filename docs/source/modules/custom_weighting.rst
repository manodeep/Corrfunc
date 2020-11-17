.. _custom_weighting:

Implementing Custom Weight Functions
====================================

``Corrfunc`` supports custom weight functions.  On this page we describe
the recommended procedure for writing your own.  When in doubt, follow
the example of ``pair_product``.

First, see :ref:`weighted_correlations` for basic usage of ``Corrfunc``'s weight features.

The steps are:

#. Add a type to the ``weight_method_t`` enum in ``utils/defs.h`` (something like ``MY_WEIGHT_SCHEME=1``).

#. Determine how many weights per particle your scheme needs, and add a case to the switch-case block in ``get_num_weights_by_method()`` in ``utils/defs.h``.  ``Corrfunc`` supports up to ``MAX_NUM_WEIGHTS=10`` weights per particle; most schemes will simply need 1.  To provide multiple weights per particle via the Python interface, simply pass a ``weights`` array of shape ``(N_WEIGHTS_PER_PARTICLE, N_PARTICLES)``.

#. Add an ``if`` statement that maps a string name (like "my_weight_scheme") to the ``weight_method_t`` (which you created above) in ``get_weight_method_by_name()`` in ``utils/defs.h``.

#. Write a function in ``utils/weight_functions.h.src`` that returns the weight for a particle pair, given the weights for the two particles.  The weights for each particle are packed in a ``const pair_struct_DOUBLE`` struct, which also contains the pair separation.  You must write one function for every instruction set you wish to support.  This can be quite easy for simple weight schemes; the four functions for ``pair_product`` are:

.. code-block:: c

    /*
     * The pair weight is the product of the particle weights
     */
    static inline DOUBLE pair_product_DOUBLE(const pair_struct_DOUBLE *pair){
        return pair->weights0[0].d*pair->weights1[0].d;
    }

    #ifdef __AVX512F__
    static inline AVX512_FLOATS avx512_pair_product_DOUBLE(const pair_struct_DOUBLE *pair){
        return AVX512_MULTIPLY_FLOATS(pair->weights0[0].a512, pair->weights1[0].a512);
    }
    #endif

    #ifdef __AVX__
    static inline AVX_FLOATS avx_pair_product_DOUBLE(const pair_struct_DOUBLE *pair){
        return AVX_MULTIPLY_FLOATS(pair->weights0[0].a, pair->weights1[0].a);
    }
    #endif

    #ifdef __SSE4_2__
    static inline SSE_FLOATS sse_pair_product_DOUBLE(const pair_struct_DOUBLE *pair){
        return SSE_MULTIPLY_FLOATS(pair->weights0[0].s, pair->weights1[0].s);
    }
    #endif

See ``utils/avx512_calls.h``, ``utils/avx_calls.h`` and ``utils/sse_calls.h`` for the lists of available vector instructions.

5. For each function you wrote in the last step, add a case to the switch-case
   block in the appropriate dispatch function in
   ``utils/weight_functions.h.src``.  If you wrote a weighting function for all
   four instruction sets, then you'll need to add the corresponding function to
   ``get_weight_func_by_method_DOUBLE()``,
   ``get_avx512_weight_func_by_method_DOUBLE``,
   ``get_avx_weight_func_by_method_DOUBLE()``,
   and  ``get_sse_weight_func_by_method_DOUBLE()``.

#. Done!  Your weight scheme should now be accessible through the Python and C interfaces via the name ("my_weight_scheme") that you specified above.  The output will be accessible in the ``weightavg`` field of the ``results`` array.

Pair counts (i.e. the ``npairs`` field in the ``results`` array)
are never affected by weights.  For theory functions like :py:mod:`Corrfunc.theory.xi` and :py:mod:`Corrfunc.theory.wp`
that actually return a clustering statistic, the statistic is weighted.
For ``pair_product``, the random distribution used to compute the
expected bin weight from an unclustered particle set (the ``RR`` term)
is taken to be a spatially uniform particle set where every particle
has the mean weight.  See :ref:`weighted_rr` for more discussion.
This behavior (automatically returning weighted clustering statistics)
is only implemented for ``pair_product``, since that is the only weighting
method for which we know the desired equivalent random distribution.
Custom weighting methods can implement similar behavior by modifying
``countpairs_xi_DOUBLE()`` in ``theory/xi/countpairs_xi_impl.c.src`` and
``countpairs_wp_DOUBLE()`` in ``theory/wp/countpairs_wp_impl.c.src``.
