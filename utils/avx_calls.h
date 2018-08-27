/* File: avx_calls.h */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>

#ifdef __cplusplus
extern "C" {
#endif

#include "function_precision.h" 

#define PREFETCH(mem)        asm ("prefetcht0 %0"::"m"(mem))

#if defined(__GNUC__) || defined(__GNUG__)
#define AVX_BIT_COUNT_INT(X)                __builtin_popcount(X)
#else
#define AVX_BIT_COUNT_INT(X)                _popcnt32(X)
#endif


#define AVX_BIT_COUNT_LONG(X)               _popcnt64(X)
#define AVX_BIT_COUNT_UNSIGNED_INT(X)       _mm_popcnt_u32(X)
#define AVX_BIT_COUNT_UNSIGNED_LONG(X)      _mm_popcnt_u64(X)
#define AVX_SET_INT(X)                      _mm256_set1_epi32(X)


#ifndef DOUBLE_PREC

#define DOUBLE                           float  
#define AVX_NVEC                         8    
#define AVX_INTS                         __m256i
#define AVX_FLOATS                       __m256

#define AVX_SETZERO_FLOAT()              _mm256_setzero_ps()
    
#define AVX_LOAD_FLOATS_UNALIGNED(X)     _mm256_loadu_ps(X)
#define AVX_LOAD_FLOATS_ALIGNED(X)       _mm256_load_ps(X)
#define AVX_MULTIPLY_FLOATS(X,Y)         _mm256_mul_ps(X,Y)
#define AVX_DIVIDE_FLOATS(X,Y)           _mm256_div_ps(X,Y)
#define AVX_SUBTRACT_FLOATS(X,Y)         _mm256_sub_ps(X,Y)
#define AVX_ADD_FLOATS(X,Y)              _mm256_add_ps(X,Y)
#define AVX_SQRT_FLOAT(X)                _mm256_sqrt_ps(X)
#define AVX_SVML_SQRT_FLOAT(X)           _mm256_svml_sqrt_ps(X)
#define AVX_TRUNCATE_FLOAT_TO_INT(X)     _mm256_cvttps_epi32(X)
#define AVX_STORE_FLOATS_TO_MEMORY(X,Y)  _mm256_storeu_ps(X,Y)
#define AVX_SQUARE_FLOAT(X)              _mm256_mul_ps(X,X)
#define AVX_LOG_FLOAT(X)                 _mm256_log_ps(X)
#define AVX_LOG10_FLOAT(X)                 _mm256_log10_ps(X)
#define AVX_LOG2_FLOAT(X)                _mm256_log2_ps(X)
#define AVX_RECIPROCAL_FLOATS(X)         _mm256_rcp_ps(X)

#define AVX_BROADCAST_FLOAT(X)           _mm256_broadcast_ss(X)
#define AVX_SET_FLOAT(X)                 _mm256_set1_ps(X)


    // X OP Y
#define AVX_COMPARE_FLOATS(X,Y,OP)        _mm256_cmp_ps(X,Y,OP)
#define AVX_BITWISE_AND(X,Y)              _mm256_and_ps(X,Y)
#define AVX_BITWISE_OR(X,Y)               _mm256_or_ps(X,Y)
#define AVX_XOR_FLOATS(X,Y)               _mm256_xor_ps(X,Y)
#define AVX_AND_NOT(X,Y)                  _mm256_andnot_ps(X,Y)  //~X & Y

//MoveMask
#define AVX_TEST_COMPARISON(X)            _mm256_movemask_ps(X)

#define AVX_BLEND_FLOATS_WITH_MASK(FALSEVALUE,TRUEVALUE,MASK) _mm256_blendv_ps(FALSEVALUE,TRUEVALUE,MASK)
#define AVX_MASKSTORE_FLOATS(dest, mask, source)   _mm256_maskstore_ps(dest, mask, source)

//Trig
#ifdef  __INTEL_COMPILER
#define AVX_ARC_COSINE(X, order)                 _mm256_acos_ps(X)
#else
    //Other compilers do not have the vectorized arc-cosine
#define AVX_ARC_COSINE(X, order)                  inv_cosine_avx(X, order)
#endif

    //Max
#define AVX_MAX_FLOATS(X,Y)               _mm256_max_ps(X,Y)


  //Absolute value
#define AVX_ABS_FLOAT(X)                  _mm256_max_ps(_mm256_sub_ps(_mm256_setzero_ps(), X), X)
  
    //Casting (does not actual convert between types)
#define AVX_CAST_FLOAT_TO_INT(X)          _mm256_castps_si256(X)
#define AVX_CAST_INT_TO_FLOAT(X)          _mm256_castsi256_ps(X)

    //Streaming store
#define AVX_STREAMING_STORE_FLOATS(X,Y)   _mm256_stream_ps(X,Y)
#define AVX_STREAMING_STORE_INTS(X,Y)     _mm256_stream_si256(X,Y)

#else //DOUBLE PRECISION CALCULATIONS
  
#define DOUBLE                           double
#define AVX_NVEC                         4    
#define AVX_INTS                         __m128i
#define AVX_FLOATS                       __m256d

#define AVX_SETZERO_FLOAT()              _mm256_setzero_pd()    

#define AVX_LOAD_FLOATS_UNALIGNED(X)     _mm256_loadu_pd(X)
#define AVX_LOAD_FLOATS_ALIGNED(X)       _mm256_load_pd(X)
#define AVX_MULTIPLY_FLOATS(X,Y)         _mm256_mul_pd(X,Y)
#define AVX_DIVIDE_FLOATS(X,Y)           _mm256_div_pd(X,Y)
#define AVX_SUBTRACT_FLOATS(X,Y)         _mm256_sub_pd(X,Y)
#define AVX_ADD_FLOATS(X,Y)              _mm256_add_pd(X,Y)
#define AVX_SQRT_FLOAT(X)                _mm256_sqrt_pd(X)
#define AVX_SVML_SQRT_FLOAT(X)           _mm256_svml_sqrt_pd(X)
#define AVX_TRUNCATE_FLOAT_TO_INT(X)     _mm256_cvttpd_epi32(X)
#define AVX_STORE_FLOATS_TO_MEMORY(X,Y)  _mm256_storeu_pd(X,Y)
#define AVX_SQUARE_FLOAT(X)              _mm256_mul_pd(X,X)
#define AVX_LOG_FLOAT(X)                 _mm256_log_pd(X)
#define AVX_LOG2_FLOAT(X)                _mm256_log2_pd(X)
#define AVX_LOG10_FLOAT(X)                _mm256_log10_pd(X)
#define AVX_RECIPROCAL_FLOATS(X)         _mm256_rcp_pd(X)

    // X OP Y
#define AVX_COMPARE_FLOATS(X,Y,OP)        _mm256_cmp_pd(X,Y,OP)
#define AVX_BITWISE_AND(X,Y)              _mm256_and_pd(X,Y)
#define AVX_BITWISE_OR(X,Y)               _mm256_or_pd(X,Y)
#define AVX_XOR_FLOATS(X,Y)               _mm256_xor_pd(X,Y)
#define AVX_AND_NOT(X,Y)                  _mm256_andnot_pd((X),(Y))  //~X & Y

#define AVX_BROADCAST_FLOAT(X)            _mm256_broadcast_sd(X)
#define AVX_SET_FLOAT(X)                  _mm256_set1_pd(X)
//MoveMask
#define AVX_TEST_COMPARISON(X)            _mm256_movemask_pd(X)

#define AVX_BLEND_FLOATS_WITH_MASK(FALSEVALUE,TRUEVALUE,MASK) _mm256_blendv_pd(FALSEVALUE,TRUEVALUE,MASK)
#define AVX_MASKSTORE_FLOATS(dest, mask, source)   _mm256_maskstore_pd(dest, mask, source)

//Trig
#ifdef  __INTEL_COMPILER
#define AVX_ARC_COSINE(X, order)                 _mm256_acos_pd(X)
#else
#define AVX_ARC_COSINE(X, order)                  inv_cosine_avx(X, order)
#endif

    //Max
#define AVX_MAX_FLOATS(X,Y)               _mm256_max_pd(X,Y)

  //Absolute value
#define AVX_ABS_FLOAT(X)                  _mm256_max_pd(_mm256_sub_pd(_mm256_setzero_pd(), X), X)
  
    //Casting (does not actual convert between types)
#define AVX_CAST_FLOAT_TO_INT(X)          _mm256_castpd_si256(X)
#define AVX_CAST_INT_TO_FLOAT(X)          _mm256_castsi256_pd(X)

    //Streaming store
#define AVX_STREAMING_STORE_FLOATS(X,Y)   _mm256_stream_pd(X,Y)
#define AVX_STREAMING_STORE_INTS(X,Y)     _mm_stream_si128(X,Y)

#endif //DOUBLE_PREC

#ifndef  __INTEL_COMPILER
#include "fast_acos.h"
    
static inline AVX_FLOATS inv_cosine_avx(const AVX_FLOATS X, const int order)
{
    union cos{
        AVX_FLOATS m;
        DOUBLE x[AVX_NVEC];
    };
    union cos union_costheta;
    union cos union_returnvalue;
    union_costheta.m = X;
    const DOUBLE minus_one = (DOUBLE) -1.0;
    const DOUBLE one = (DOUBLE) 1.0;

    //Force everything to be in range [0,1]
    for(int ii=0;ii<AVX_NVEC;ii++) {
        const DOUBLE costheta = union_costheta.x[ii];
        union_costheta.x[ii] = costheta <= minus_one ? minus_one:costheta;
        union_costheta.x[ii] = costheta >= one ? one:costheta;
    }
    
    if(order == 0) {
        for(int ii=0;ii<AVX_NVEC;ii++) {
            const DOUBLE costheta = union_costheta.x[ii];
            union_returnvalue.x[ii] = ACOS(costheta);
        }
    } else {
        //fast acos
        /*Taken from associated C++ code in http://www.geometrictools.com/GTEngine/Include/Mathematics/GteACosEstimate.h*/
        for(int ii=0;ii<AVX_NVEC;ii++) {
            union_returnvalue.x[ii] = FAST_ACOS(union_costheta.x[ii]);
        }
    }
    return union_returnvalue.m;
  }


#endif
  
  //The three different unions used
  //for computing rpavg and weightavg
  union int8 {
    AVX_INTS m_ibin;
    int ibin[AVX_NVEC];
  };
  union float8{
    AVX_FLOATS m_Dperp;
    DOUBLE Dperp[AVX_NVEC];
  };
  union float8_weights{
    AVX_FLOATS m_weights;
    DOUBLE weights[AVX_NVEC];
  };

#ifdef DOUBLE_PREC    
#define CHECK_AND_FAST_DIVIDE_AVX(result, numerator, denominator, fast_divide_and_NR_steps)                      { \
        /* For double precision floats */                               \
        if (fast_divide_and_NR_steps == 0) {                            \
            result = AVX_DIVIDE_FLOATS(numerator, denominator);         \
            /* The divide is the actual operation we need */            \
            /* but divides are about 10x slower than multiplies. So, I am replacing it */ \
            /* with a approximate reciprocal in floating point */       \
            /* + 2 iterations of newton-raphson in case of DOUBLE */    \
        } else {                                                        \
            /* following blocks do an approximate reciprocal followed by two iterations of Newton-Raphson */ \
            const __m128 float_tmp1 =  _mm256_cvtpd_ps(denominator);/* convert double to float -> not avx_floats := _m256d */ \
            /*(convert 4 doubles into 4 floats -> use half of available 256 bit SIMD registers) */ \
            __m128 float_inv_tmp1 = _mm_rcp_ps(float_tmp1);/* intrinsic for 128 bit float approximate reciprocal */ \
            const AVX_FLOATS rc = _mm256_cvtps_pd(float_inv_tmp1);/* convert back to double */ \
            /* We have the double->float->approx. reciprocal->double process done. */ \
            /* Now improve the accuracy of the divide with newton-raphson. */ \
            const AVX_FLOATS two = AVX_SET_FLOAT((DOUBLE) 2.0);         \
            AVX_FLOATS rc_iter = rc;                                    \
            /* Do NewtonRaphson iterations */                           \
            for(unsigned int _ii=0;_ii<fast_divide_and_NR_steps;_ii++) { \
                rc_iter = AVX_MULTIPLY_FLOATS(rc_iter,                  \
                                              AVX_SUBTRACT_FLOATS(two,  \
                                                                  AVX_MULTIPLY_FLOATS(denominator, rc_iter))); /*2.0 - l^2*rc */ \
            }                                                           \
            result = AVX_MULTIPLY_FLOATS(numerator, rc_iter);           \
        } /* end of FAST_DIVIDE */                                      \
    }
#else
#define CHECK_AND_FAST_DIVIDE_AVX(result, numerator, denominator, fast_divide_and_NR_steps)                      { \
        /* single precision floats */                                   \
        if (fast_divide_and_NR_steps == 0) {                            \
            result = AVX_DIVIDE_FLOATS(numerator, denominator);         \
            /* The divide is the actual operation we need */            \
            /* but divides are about 10x slower than multiplies. So, I am replacing it */ \
            /* with a approximate reciprocal in floating point */       \
            /* + 2 iterations of newton-raphson in case of DOUBLE */    \
        } else {                                                        \
            /* following blocks do an approximate reciprocal followed by two iterations of Newton-Raphson */ \
            const AVX_FLOATS rc  = _mm256_rcp_ps(denominator);/* intrinsic for 256 bit approximate reciprocal */ \
            /* We have the double->float->approx. reciprocal->double process done. */ \
            /* Now improve the accuracy of the divide with newton-raphson. */ \
            const AVX_FLOATS two = AVX_SET_FLOAT((DOUBLE) 2.0);         \
            AVX_FLOATS rc_iter = rc;                                    \
            /* Do NewtonRaphson iterations */                           \
            for(unsigned int _ii=0;_ii<fast_divide_and_NR_steps;_ii++) {             \
                rc_iter = AVX_MULTIPLY_FLOATS(rc_iter,                  \
                                              AVX_SUBTRACT_FLOATS(two,  \
                                                                  AVX_MULTIPLY_FLOATS(denominator, rc_iter))); /*2.0 - l^2*rc */ \
            }                                                           \
            result = AVX_MULTIPLY_FLOATS(numerator, rc_iter);           \
        } /* end of FAST_DIVIDE */                                      \
    }
#endif /* end of DOUBLE_PREC for defining check_and_fast_divide macro */


#ifdef __cplusplus
}
#endif
