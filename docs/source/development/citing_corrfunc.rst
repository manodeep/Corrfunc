:orphan:

.. _citing_corrfunc:

==============================================
License and Citation Information
==============================================

Citing Corrfunc
------------------


If you use ``Corrfunc`` for research, please cite using the MNRAS code paper with the following
bibtex entry:

::

   @ARTICLE{2020MNRAS.491.3022S,
       author = {{Sinha}, Manodeep and {Garrison}, Lehman H.},
       title = "{CORRFUNC - a suite of blazing fast correlation functions on
       the CPU}",
       journal = {\mnras},
       keywords = {methods: numerical, galaxies: general, galaxies:
       haloes, dark matter, large-scale structure of Universe, cosmology:
       theory},
       year = "2020",
       month = "Jan",
       volume = {491},
       number = {2},
       pages = {3022-3041},
       doi = {10.1093/mnras/stz3157},
       adsurl =
       {https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.3022S},
       adsnote = {Provided by the SAO/NASA
       Astrophysics Data System}
   }


The MNRAS paper (also on `arXiv:1911.03545
<https://arxiv.org/abs/1911.03545>`_) targets ``Corrfunc v2.0.0``. If you are
using ``Corrfunc v2.3.0`` or later, **and** you benefit from the
enhanced vectorised kernels, then please additionally cite this paper:

::

      @InProceedings{10.1007/978-981-13-7729-7_1,
          author="Sinha, Manodeep and Garrison, Lehman",
          editor="Majumdar, Amit and Arora, Ritu",
          title="CORRFUNC: Blazing Fast Correlation Functions with AVX512F SIMD Intrinsics",
          booktitle="Software Challenges to Exascale Computing",
          year="2019",
          publisher="Springer Singapore",
          address="Singapore",
          pages="3--20",
          isbn="978-981-13-7729-7",
          url={https://doi.org/10.1007/978-981-13-7729-7_1}
      }



Corrfunc License
---------------------

Corrfunc comes with a MIT LICENSE - see the LICENSE file.

Copyright (C) 2014 Manodeep Sinha (manodeep@gmail.com)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to
deal in the Software without restriction, including without limitation the
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
