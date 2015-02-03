/* File: cellarray.h */
/*
		This file is a part of the Corrfunc package
		Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
		License: MIT LICENSE. See LICENSE file under the top-level
		directory at https://bitbucket.org/manodeep/corrfunc/
*/

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

#ifdef __cplusplus
extern "C" {
#endif

#include "function_precision.h" //definition of DOUBLE depending on Makefile flag -DOUBLE_PREC

#define NLATMAX   100      /* maximum grid dimension in X-Y plane */

typedef struct{
  DOUBLE *x;
  DOUBLE *y;
  DOUBLE *z;
  int64_t nelements;
} cellarray;

typedef struct{
  DOUBLE *x;
  DOUBLE *y;
  DOUBLE *z;
  DOUBLE *cz;
  int nelements;
  int nallocated;
} cellarray_mocks;


#ifdef __cplusplus
}
#endif
