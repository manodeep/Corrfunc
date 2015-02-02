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
	DOUBLE **pN;
	DOUBLE rmax;
	unsigned int nbin;
	unsigned int nc;
	unsigned int num_pN;
} results_countspheres;


results_countspheres * countspheres(const int64_t np, const DOUBLE * restrict X, const DOUBLE * restrict Y, const DOUBLE * restrict Z,
																		const double rmax, const unsigned int nbin, const unsigned int nc,
																		const unsigned int num_pN,
																		unsigned long seed) __attribute__((warn_unused_result));

void free_results_countspheres(results_countspheres **results);	

#ifdef __cplusplus
}
#endif
	
