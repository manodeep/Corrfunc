#pragma once

#include "function_precision.h" //for definition of DOUBLE
#include <inttypes.h> //for uint64_t

void countpairs_rp_pi(const int ND1,
											const DOUBLE *x1, const DOUBLE *y1, const DOUBLE *z1,
											const int ND2,
											const DOUBLE *x2,  const DOUBLE *y2, const DOUBLE *z2,
											const DOUBLE xmin, const DOUBLE xmax,
											const DOUBLE ymin, const DOUBLE ymax,
											const DOUBLE zmin, const DOUBLE zmax,
											const int autocorr,
											const double rpmax,
#ifdef USE_OMP
											const int nthreads,
#endif
											const int nrpbin,const double * restrict rupp,
											const double pimax, const int npibin);
