#pragma once

#include "function_precision.h"
#include "cellarray.h"

cellarray * gridlink(const int np,
					 const DOUBLE *x,const DOUBLE *y,const DOUBLE *z,
					 const DOUBLE xmin, const DOUBLE xmax,
					 const DOUBLE ymin, const DOUBLE ymax,
					 const DOUBLE zmin, const DOUBLE zmax,
					 const DOUBLE rmax,
					 const int xbin_refine_factor,
					 const int ybin_refine_factor,
					 const int zbin_refine_factor,
					 int *nlattice_x,
					 int *nlattice_y,
					 int *nlattice_z) __attribute__((warn_unused_result)); 
