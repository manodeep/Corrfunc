#pragma once

#define PI_OVER_180       0.01745329251
#define INV_PI_OVER_180   57.2957795131

#define REGISTER_WIDTH 256  //cpu supports avx instructions
#define NVECF  8  //8 floats per ymm register
#define NVECD  4  //4 doubles per ymm register

#ifdef DOUBLE_PREC
#define DOUBLE double
#define NVEC   NVECD
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
#else
#define DOUBLE float
#define NVEC   NVECF
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
#endif
