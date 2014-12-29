#pragma once

#include "function_precision.h" //for definition of DOUBLE
#include <inttypes.h> //for uint64_t

void countpairs(const int ND1, const DOUBLE * const X1, const DOUBLE * const Y1, const DOUBLE  * const Z1,
								const int ND2, const DOUBLE * const X2, const DOUBLE * const Y2, const DOUBLE  * const Z2,
								const DOUBLE xmin, const DOUBLE xmax,
								const DOUBLE ymin, const DOUBLE ymax,
								const DOUBLE zmin, const DOUBLE zmax,
								const int autocorr,
								const double rpmax,
#ifdef USE_OMP
								const int numthreads,
#endif
								const int nrpbin, const double * restrict rupp);

