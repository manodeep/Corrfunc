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

#ifdef __cplusplus
}
#endif
