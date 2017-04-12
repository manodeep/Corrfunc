.. _generate_bins:

*******************************************
Specifying the separation bins in Corrfunc
*******************************************

All of the python extensions for `Corrfunc` accept
either a filename or an array for specifying the
:math:`r_p` or :math:`\theta`. 

Manually creating a file with arbitrary bins
--------------------------------------------
This manual method lets you specify the most
generic bins. The upper-edge of one bin does
not have to be the lower-edge of the next, the
bins themselves can have arbitrary widths, the
smallest bin can start from 0.0. 

* Open a text editor with a new file
* Add two columns per bin you want, the first
  column should be low-edge of the bin while
  the second column should be the high-edge
  of the bin. Like so:

::
   
    0.10     0.15

* Now add as many such lines as the number of bins you
  want. Here is a valid example:

::
  
     0.10     0.15
     0.20     0.50
     0.50     5.00

This example specifies 3 bins, with the individual
bin limits specified on each line. Notice that the
bins need not be continuous, and the width of each
bin can be independently specified.
  
**NOTE** Make sure that the bins are in increasing
order -- smallest bin first, then the next smallest
bin and so on up to the largest bin.

Specifying bins as an array
---------------------------

You can specify the bins using ``numpy.linspace`` or
``numpy.logspace``. Note, in this case the bins are constrained
to be continuous (unlike the preceeding case). 

.. code:: python

          import numpy as np
          rmin = 0.1
          rmax = 10.0
          nbins = 20
          rbins = np.linspace(rmin, rmax, nbins + 1)
          log_rbins = np.logspace(np.log10(rmin), np.log10(rmax), nbins + 1)
