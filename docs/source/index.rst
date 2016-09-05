=======================
Corrfunc Documentation
=======================


Corrfunc is a set of high-performance routines to measure clustering statistics.
The main features of Corrfunc are:

* **Fast** All theory pair-counting is at least an order of magnitude faster than all existing public codes. Particularly suited for MCMC.
* **Python Extensions** Python extensions allow you to do the compute-heavy bits using C while retaining all of the user-friendliness of python.
* **Modular** The code is written in a modular fashion and is easily extensible to compute arbitrary clustering statistics.
* **Future-proof** As I get access to newer instruction-sets, the codes will
  get updated to use the latest and greatest CPU features. 

The source code is publicly available at https://github.com/manodeep/Corrfunc.

******************
Getting Started
******************

.. toctree::
   :maxdepth: 1

   install
   modules/index
   development/index
              
*********
Reference
*********

.. toctree::
   :maxdepth: 1

   development/corrfunc_full_api

*********************
License and Credits
*********************

.. toctree::
   :maxdepth: 1

   development/contributors
   development/citing_corrfunc

