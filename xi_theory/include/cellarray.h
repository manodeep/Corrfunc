#pragma once

#include <stdio.h>
#include <stdlib.h>

#include "function_precision.h" //definition of DOUBLE depending on Makefile flag -DOUBLE_PREC

#define NLATMAX   100      /* maximum grid dimension in X-Y plane */


//pointer alignment is 8 bytes. So the 4 bytes will be padded and wasted anyway ->
//changing nelements to a int64_t is fine 
typedef struct{
  DOUBLE *x;
  DOUBLE *y;
  DOUBLE *z;
  int nelements;
} cellarray;


