/* File: countpairs.h */
/*
		This file is a part of the Corrfunc package
		Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
		License: MIT LICENSE. See LICENSE file under the top-level
		directory at https://bitbucket.org/manodeep/corrfunc/
*/

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "function_precision.h" //for definition of DOUBLE
#include <inttypes.h> //for uint64_t


//define the results structure
typedef struct{
	uint64_t *npairs;
	DOUBLE *rupp;
	DOUBLE *rpavg;
	int nbin;
} results_countpairs;

results_countpairs * countpairs(const int64_t ND1, const DOUBLE * const X1, const DOUBLE * const Y1, const DOUBLE  * const Z1,
								const int64_t ND2, const DOUBLE * const X2, const DOUBLE * const Y2, const DOUBLE  * const Z2,
#ifdef USE_OMP
								const int numthreads,
#endif
								const int autocorr,
								const char *binfile) __attribute__((warn_unused_result));


void free_results(results_countpairs **results);

#ifdef __cplusplus
}
#endif
	
