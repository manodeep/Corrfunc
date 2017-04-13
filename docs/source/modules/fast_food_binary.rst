.. _fast_food_binary:

************************
Fast-food binary format
************************

The fast-food format is a fortran binary format -- all fields are surrounded
with 4 bytes padding. These value of these padding bytes
is the number of bytes of data contained in between the padding bytes. For
example, to write out ``20 bytes of data`` in
a fast-food file format would require a total of ``4+20+4=28`` bytes. The first
and last 4 bytes of the file will contain the value 20 --
showing that 20 bytes of real data are contained in between the two paddings.

The ``fast-food`` file consists of a header:

.. code:: C
          
          int idat[5];
          float fdat[9];
          float znow;

For the purposes of these correlation function codes, the only useful quantity
is ``idat[1]`` which contains ``N`` -- the number of particles
in the  file. The rest can simply filled with `0`.

After this header, the actual ``X/Y/Z`` values are stored. The first 4
bytes after the header contains ``4*N`` for float precision or
``8*N`` for  double precision where ``N=idat[1]``, is the number
of particles in the file. After all of the ``X`` values there will 
be another 4 bytes containing ``4*N`` or ``8*N``.

.. note:: Even when the ``X/Y/Z`` arrays are written out in double-precision, the padding is still 4 bytes.
          The blocks for ``Y/Z`` similarly follow after the ``X`` block.
