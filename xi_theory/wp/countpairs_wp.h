#pragma once

#include "function_precision.h"
#include <inttypes.h>

//define the results structure
typedef struct{
	uint64_t *npairs;
	DOUBLE *wp;
	DOUBLE *rupp;
#ifdef OUTPUT_RPAVG
	DOUBLE *rpavg;
#endif	
	DOUBLE pimax;
	int nbin;
} results_countpairs_wp;

results_countpairs_wp *countpairs_wp(const int64_t ND1, const DOUBLE * restrict X1, const DOUBLE * restrict Y1, const DOUBLE * restrict Z1,
									 const double boxsize, 
#ifdef USE_OMP
									 const int numthreads,
#endif
									 const char *binfile,
									 const double pimax) __attribute__((warn_unused_result));


void free_results_wp(results_countpairs_wp **results);
