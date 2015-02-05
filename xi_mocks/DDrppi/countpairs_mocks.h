/* File: countpairs_data.h */
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
  DOUBLE pimax;
  int nbin;
  int npibin;
} results_countpairs_data;

results_countpairs_data * countpairs_data(const int64_t N1, const DOUBLE *theta1, const DOUBLE *phi1, const DOUBLE *d1,
										  const int64_t N2, const DOUBLE *theta2, const DOUBLE *phi2, const DOUBLE *d2,
#ifdef USE_OMP  
										  const int numthreads,
#endif
										  const int autocorr,
										  const char *binfile,
										  const DOUBLE pimax) __attribute__((warn_unused_result));


void free_results_data(results_countpairs_data **results);

#ifdef __cplusplus
}
#endif
	
