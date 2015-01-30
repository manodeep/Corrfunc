/* File: countpairs_xi.h */
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
	DOUBLE *xi;
	DOUBLE *rupp;
	DOUBLE *rpavg;
	int nbin;
} results_countpairs_xi;

results_countpairs_xi *countpairs_xi(const int64_t ND1, DOUBLE * restrict X1, DOUBLE * restrict Y1, DOUBLE * restrict Z1,
									 const double boxsize, 
#ifdef USE_OMP
									 const int numthreads,
#endif
									 const char *binfile) __attribute__((warn_unused_result));
									 

void free_results_xi(results_countpairs_xi **results);

#ifdef __cplusplus
}
#endif
	
