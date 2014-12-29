#pragma once

#include "function_precision.h"
#include <inttypes.h>

void countpairs_wp(const int ND1, const DOUBLE * restrict X1, const DOUBLE * restrict Y1, const DOUBLE * restrict Z1,
									 const double boxsize, const DOUBLE pimax, 
#ifdef USE_OMP
									 const int numthreads,
#endif
									 uint64_t * restrict npair,
#ifdef OUTPUT_RPAVG
									 DOUBLE *rpavg,
#endif
									 const double * restrict rupp,const int nbin);

