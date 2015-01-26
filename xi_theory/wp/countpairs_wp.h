#pragma once


#ifdef __cplusplus
extern "C" {
#endif

#include "function_precision.h"
#include <inttypes.h>

//define the results structure
typedef struct{
	uint64_t *npairs;
	DOUBLE *wp;
	DOUBLE *rupp;
	DOUBLE *rpavg;
	DOUBLE pimax;
	int nbin;
} results_countpairs_wp;

results_countpairs_wp *countpairs_wp(const int64_t ND1, DOUBLE * restrict X1, DOUBLE * restrict Y1, DOUBLE * restrict Z1,
									 const double boxsize, 
#ifdef USE_OMP
									 const int numthreads,
#endif
									 const char *binfile,
									 const double pimax) __attribute__((warn_unused_result));


void free_results_wp(results_countpairs_wp **results);

#ifdef __cplusplus
}
#endif
	
