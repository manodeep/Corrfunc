/* File: function_precision.h */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#pragma once

#ifdef __cplusplus
extern "C" {
#endif


#ifndef M_PI
#define M_PI            3.14159265358979323846264338327950288
#endif

#define PI_OVER_180       0.017453292519943295769236907684886127134428718885417254560971
#define INV_PI_OVER_180   57.29577951308232087679815481410517033240547246656432154916024

//Define the Macros
#ifdef DOUBLE_PREC
#define COSD(X)            cos(X*PI_OVER_180)
#define SIND(X)            sin(X*PI_OVER_180)
#else
#define COSD(X)            cosf(X*PI_OVER_180)
#define SIND(X)            sinf(X*PI_OVER_180)
#endif


#include <float.h>
    
#ifdef DOUBLE_PREC
#define DOUBLE double
#define REAL_FORMAT "lf"
#define ZERO   0.0
#define SQRT   sqrt
#define LOG    log
#define LOG10  log10
#define LOG2   log2
#define FABS   fabs
#define COS    cos
#define SIN    sin
#define ACOS   acos
#define ASIN   asin
#define POW    pow
#define ABS    fabs
#define MAX_POSITIVE_FLOAT DBL_MAX
#else
#define DOUBLE float
#define REAL_FORMAT "f"
#define ZERO   0.0f
#define SQRT   sqrtf
#define LOG    logf
#define LOG10  log10f
#define LOG2   log2f
#define FABS   fabsf
#define COS    cosf
#define SIN    sinf
#define ACOS   acosf
#define ASIN   asinf
#define POW    powf
#define ABS    fabsf
#define MAX_POSITIVE_FLOAT FLT_MAX
#endif

#ifdef __cplusplus
}
#endif
