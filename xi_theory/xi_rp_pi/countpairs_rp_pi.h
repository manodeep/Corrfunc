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
	DOUBLE pimax;
	int nbin;
	int npibin;
} results_countpairs_rp_pi;

results_countpairs_rp_pi * countpairs_rp_pi(const int64_t ND1, const DOUBLE *X1, const DOUBLE *Y1, const DOUBLE *Z1,
											const int64_t ND2, const DOUBLE *X2, const DOUBLE *Y2, const DOUBLE *Z2,
#ifdef USE_OMP
											const int numthreads,
#endif
											const int autocorr,
											const char *binfile,
											const double pimax)  __attribute__((warn_unused_result));

void free_results_rp_pi(results_countpairs_rp_pi **results);
