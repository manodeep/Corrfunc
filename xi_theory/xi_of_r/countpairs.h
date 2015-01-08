#pragma once

#include "function_precision.h" //for definition of DOUBLE
#include <inttypes.h> //for uint64_t


//define the results structure
typedef struct{
	uint64_t *npairs;
	DOUBLE *rupp;
#ifdef OUTPUT_RPAVG
	DOUBLE *rpavg;
#endif	
	int nbin;
} results_countpairs;

results_countpairs * countpairs(const int64_t ND1, const DOUBLE * const X1, const DOUBLE * const Y1, const DOUBLE  * const Z1,
								const int64_t ND2, const DOUBLE * const X2, const DOUBLE * const Y2, const DOUBLE  * const Z2,
#ifdef USE_OMP
								const int numthreads,
#endif
								const int autocorr,
								const char *binfile) __attribute__((warn_unused_result));


void free_results(results_countpairs **results)  __attribute__((noinline));
