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
#include <stdint.h>
#include <immintrin.h>

#ifdef __cplusplus
extern "C" {
#endif

#include "function_precision.h" 

#define PREFETCH(mem)        asm ("prefetcht0 %0"::"m"(mem))

#if defined(__GNUC__) || defined(__GNUG__)
#define AVX512_BIT_COUNT_INT(X)                __builtin_popcount(X)
#else
#define AVX512_BIT_COUNT_INT(X)                _popcnt32(X)
#endif


#define AVX512_BIT_COUNT_LONG(X)               _popcnt64(X)
#define AVX512_BIT_COUNT_UNSIGNED_INT(X)       _mm_popcnt_u32(X)
#define AVX512_BIT_COUNT_UNSIGNED_LONG(X)      _mm_popcnt_u64(X)


#define AVX512_MASK_BITWISE_AND(X,Y)              _mm512_kand(X,Y)
#define AVX512_MASK_BITWISE_OR(X,Y)               _mm512_kor(X,Y)
#define AVX512_MASK_BITWISE_XOR_FLOATS(X,Y)       _mm512_kxor(X,Y)
#define AVX512_MASK_BITWISE_AND_NOT(X,Y)          _mm512_kandn(X,Y)  //~X & Y


  /* For setting up the array that contains the number bits set in a mask */
  // Can be used to generate up to 16 bit lookup tables
#   define B2(n)  n,      n+1,      n+1,      n+2
#   define B4(n)  B2(n),  B2(n+1),  B2(n+1),  B2(n+2)
#   define B6(n)  B4(n),  B4(n+1),  B4(n+1),  B4(n+2)
#   define B8(n)  B6(n),  B6(n+1),  B6(n+1),  B6(n+2)
#   define B10(n) B8(n),  B8(n+1),  B8(n+1),  B8(n+2)
#   define B12(n) B10(n), B10(n+1), B10(n+1), B10(n+2)
#   define B14(n) B12(n), B12(n+1), B12(n+1), B12(n+2)
#   define B16(n) B14(0),B14(1), B14(1),   B14(2)

#ifndef DOUBLE_PREC

#define DOUBLE                              float  
#define AVX512_NVEC                         16    
#define AVX512_MASK                         __mmask16
#define AVX512_FLOATS                       __m512

#define AVX512_INTS                            __m512i
#define AVX512_SET_INT(X)                      _mm512_set1_epi32(X)

#define AVX512_LOAD_FLOATS_UNALIGNED(X)     _mm512_loadu_ps(X)
#define AVX512_LOAD_FLOATS_ALIGNED(X)       _mm512_load_ps(X)
#define AVX512_MASK_LOAD_FLOATS_UNALIGNED(FALSEVALS, MASK, X)     _mm512_mask_loadu_ps(FALSEVALS, MASK, X)
#define AVX512_MASK_LOAD_FLOATS_ALIGNED(FALSEVALS, MASK, X)       _mm512_mask_load_ps(FALSEVALS, MASK, X)


#define AVX512_MULTIPLY_FLOATS(X,Y)         _mm512_mul_ps(X,Y)
#define AVX512_DIVIDE_FLOATS(X,Y)           _mm512_div_ps(X,Y)
#define AVX512_SUBTRACT_FLOATS(X,Y)         _mm512_sub_ps(X,Y)
#define AVX512_ADD_FLOATS(X,Y)              _mm512_add_ps(X,Y)
#define AVX512_FMA_ADD_FLOATS(X,Y,Z)            _mm512_fmadd_ps(X,Y,Z)
#define AVX512_MASK_FMA_ADD_FLOATS(X, MASK, Y, Z) _mm512_mask_fmadd_ps(X, MASK, Y, Z)

#ifdef  __INTEL_COMPILER
#define AVX512_HORIZONTAL_SUM_FLOATS(X)     _mm512_reduce_add_ps(X)
#define AVX512_MASK_HORIZONTAL_SUM_FLOATS(MASK, X)     _mm512_mask_reduce_add_ps(MASK, X)
#else
#define AVX512_HORIZONTAL_SUM_FLOATS(X)                _horizontal_sum_floats(X)
#define AVX512_MASK_HORIZONTAL_SUM_FLOATS(MASK, X)     _horizontal_mask_sum_floats(X)
#endif

#define AVX512_SQRT_FLOAT(X)                        _mm512_sqrt_ps(X)
#define AVX512_MASK_SQRT_FLOAT(FALSEVALS, MASK, X)  _mm512_mask_sqrt_ps(FALSEVALS, MASK, X)
#define AVX512_SVML_SQRT_FLOAT(X)                   _mm512_svml_sqrt_ps(X)
#define AVX512_TRUNCATE_FLOAT_TO_INT(X)             _mm512_cvttps_epi32(X)
#define AVX512_STORE_FLOATS_TO_MEMORY(X,Y)          _mm512_storeu_ps(X,Y)
#define AVX512_SQUARE_FLOAT(X)                      _mm512_mul_ps(X,X)
#define AVX512_LOG_FLOAT(X)                         _mm512_log_ps(X)
#define AVX512_LOG10_FLOAT(X)                       _mm512_log10_ps(X)
#define AVX512_LOG2_FLOAT(X)                        _mm512_log2_ps(X)
#define AVX512_RECIPROCAL_FLOATS(X)                 _mm512_rcp_ps(X)

#define AVX512_BROADCAST_FLOAT(X)                   _mm512_broadcast_ss(X)
#define AVX512_SET_FLOAT(X)                         _mm512_set1_ps(X)


    // X OP Y
#define AVX512_COMPARE_FLOATS(X, Y, OP)             _mm512_cmp_ps_mask(X, Y, OP)

  //Mask operations (new in AVX512)
#define AVX512_MASK_COMPARE_FLOATS(M, X, Y, OP)     _mm512_mask_cmp_ps_mask(M, X, Y, OP)

//MoveMask
#define AVX512_TEST_COMPARISON(X)                   _mm512_movemask_ps(X)

#define AVX512_BLEND_FLOATS_WITH_MASK(MASK, FALSE,TRUE) _mm512_mask_blend_ps(MASK, FALSE,TRUE)

#if defined(__AVX512VL__)
/* #warning blending integers with avx512vl */
#define AVX512_BLEND_INTS_WITH_MASK(MASK, FALSE,TRUE)   _mm512_mask_blend_epi32(MASK, FALSE,TRUE)
#elif defined(__AVX2__)
/* #warning blending integers with avx2 */
#define AVX512_BLEND_INTS_WITH_MASK(MASK, FALSE,TRUE) _mm256_blend_epi32(FALSE, TRUE, MASK)//AVX2
#else
#error Do not know how to implement a integer blend intrinsics on AVX512 cpu without AVX2 support. 
#endif

#define AVX512_MASKSTORE_FLOATS(dest, mask, source)   _mm512_maskstore_ps(dest, mask, source)

//Trig
#ifdef  __INTEL_COMPILER
#define AVX512_ARC_COSINE(X, order)                 _mm512_acos_ps(X)
#else
    //Other compilers do not have the vectorized arc-cosine
#define AVX512_ARC_COSINE(X, order)                  inv_cosine_avx512(X, order)
#endif

    //Max
#define AVX512_MAX_FLOATS(X,Y)               _mm512_max_ps(X,Y)


  //Absolute value
#define AVX512_ABS_FLOAT(X)                  _mm512_max_ps(_mm512_sub_ps(_mm512_setzero_ps(), X), X)
  
    //Casting (does not actual convert between types)
#define AVX512_CAST_FLOAT_TO_INT(X)          _mm512_castps_si256(X)
#define AVX512_CAST_INT_TO_FLOAT(X)          _mm512_castsi256_ps(X)

    //Streaming store
#define AVX512_STREAMING_STORE_FLOATS(X,Y)   _mm512_stream_ps(X,Y)
#define AVX512_STREAMING_STORE_INTS(X,Y)     _mm512_stream_si256(X,Y)

#else //DOUBLE PRECISION CALCULATIONS
  
#define DOUBLE                           double
#define AVX512_NVEC                         8    
#define AVX512_MASK                         __mmask8
#define AVX512_FLOATS                       __m512d
#define AVX512_INTS                         __m256i
#define AVX512_SET_INT(X)                   _mm256_set1_epi32(X)

#define AVX512_LOAD_FLOATS_UNALIGNED(X)     _mm512_loadu_pd(X)
#define AVX512_LOAD_FLOATS_ALIGNED(X)       _mm512_load_pd(X)
#define AVX512_MASK_LOAD_FLOATS_UNALIGNED(FALSEVALS, MASK, X)     _mm512_mask_loadu_pd(FALSEVALS, MASK, X)
#define AVX512_MASK_LOAD_FLOATS_ALIGNED(FALSEVALS, MASK, X)       _mm512_mask_load_pd(FALSEVALS, MASK, X)

#define AVX512_MULTIPLY_FLOATS(X,Y)         _mm512_mul_pd(X,Y)
#define AVX512_DIVIDE_FLOATS(X,Y)           _mm512_div_pd(X,Y)
#define AVX512_SUBTRACT_FLOATS(X,Y)         _mm512_sub_pd(X,Y)
#define AVX512_ADD_FLOATS(X,Y)              _mm512_add_pd(X,Y)
#define AVX512_FMA_ADD_FLOATS(X,Y,Z)            _mm512_fmadd_pd(X,Y,Z)
#define AVX512_MASK_FMA_ADD_FLOATS(X, MASK, Y, Z) _mm512_mask_fmadd_pd(X, MASK, Y, Z)

#ifdef  __INTEL_COMPILER
#define AVX512_HORIZONTAL_SUM_FLOATS(X)     _mm512_reduce_add_ps(X)
#define AVX512_MASK_HORIZONTAL_SUM_FLOATS(MASK, X)     _mm512_mask_reduce_add_ps(MASK, X)
#endif

#define AVX512_SQRT_FLOAT(X)                _mm512_sqrt_pd(X)
#define AVX512_MASK_SQRT_FLOAT(FALSEVALS, MASK, X)  _mm512_mask_sqrt_pd(FALSEVALS, MASK, X)
#define AVX512_SVML_SQRT_FLOAT(X)           _mm512_svml_sqrt_pd(X)
#define AVX512_TRUNCATE_FLOAT_TO_INT(X)     _mm512_cvttpd_epi32(X)
#define AVX512_STORE_FLOATS_TO_MEMORY(X,Y)  _mm512_storeu_pd(X,Y)
#define AVX512_SQUARE_FLOAT(X)              _mm512_mul_pd(X,X)
#define AVX512_LOG_FLOAT(X)                 _mm512_log_pd(X)
#define AVX512_LOG2_FLOAT(X)                _mm512_log2_pd(X)
#define AVX512_LOG10_FLOAT(X)                _mm512_log10_pd(X)
#define AVX512_RECIPROCAL_FLOATS(X)         _mm512_rcp_pd(X)

    // X OP Y
#define AVX512_COMPARE_FLOATS(X, Y, OP)      _mm512_cmp_pd_mask(X, Y, OP)
#define AVX512_MASK_COMPARE_FLOATS(M, X, Y, OP) _mm512_mask_cmp_pd_mask(M, X, Y, OP)

#define AVX512_BROADCAST_FLOAT(X)            _mm512_broadcast_sd(X)
#define AVX512_SET_FLOAT(X)                  _mm512_set1_pd(X)
//MoveMask
#define AVX512_TEST_COMPARISON(X)            _mm512_movemask_pd(X)

#define AVX512_BLEND_FLOATS_WITH_MASK(MASK, FALSEVALUE,TRUEVALUE) _mm512_mask_blend_pd(MASK, FALSEVALUE,TRUEVALUE)
#define AVX512_BLEND_INTS_WITH_MASK(MASK, FALSEVALUE,TRUEVALUE) _mm256_mask_blend_epi32(MASK, FALSEVALUE,TRUEVALUE)
#define AVX512_MASKSTORE_FLOATS(dest, mask, source)   _mm512_maskstore_pd(dest, mask, source)

//Trig
#ifdef  __INTEL_COMPILER
#define AVX512_ARC_COSINE(X, order)                 _mm512_acos_pd(X)
#else
#define AVX512_ARC_COSINE(X, order)                  inv_cosine_avx(X, order)
#endif

    //Max
#define AVX512_MAX_FLOATS(X,Y)               _mm512_max_pd(X,Y)

  //Absolute value
#define AVX512_ABS_FLOAT(X)                  _mm512_max_pd(_mm512_sub_pd(_mm512_setzero_pd(), X), X)
  
    //Casting (does not actual convert between types)
#define AVX512_CAST_FLOAT_TO_INT(X)          _mm512_castpd_si256(X)
#define AVX512_CAST_INT_TO_FLOAT(X)          _mm512_castsi256_pd(X)

    //Streaming store
#define AVX512_STREAMING_STORE_FLOATS(X,Y)   _mm512_stream_pd(X,Y)
#define AVX512_STREAMING_STORE_INTS(X,Y)     _mm512_stream_si512(X,Y)

#endif //DOUBLE_PREC

#ifndef  __INTEL_COMPILER
#include "fast_acos.h"
    
static inline AVX512_FLOATS inv_cosine_avx512(const AVX512_FLOATS X, const int order)
{
    union cos{
        AVX512_FLOATS m;
        DOUBLE x[AVX512_NVEC];
    };
    union cos union_costheta;
    union cos union_returnvalue;
    union_costheta.m = X;
    const DOUBLE minus_one = (DOUBLE) -1.0;
    const DOUBLE one = (DOUBLE) 1.0;

    //Force everything to be in range [0,1]
    for(int ii=0;ii<AVX512_NVEC;ii++) {
        const DOUBLE costheta = union_costheta.x[ii];
        union_costheta.x[ii] = costheta <= minus_one ? minus_one:costheta;
        union_costheta.x[ii] = costheta >= one ? one:costheta;
    }
    
    if(order == 0) {
        for(int ii=0;ii<AVX512_NVEC;ii++) {
            const DOUBLE costheta = union_costheta.x[ii];
            union_returnvalue.x[ii] = ACOS(costheta);
        }
    } else {
        //fast acos
        /*Taken from associated C++ code in http://www.geometrictools.com/GTEngine/Include/Mathematics/GteACosEstimate.h*/
        for(int ii=0;ii<AVX512_NVEC;ii++) {
            union_returnvalue.x[ii] = FAST_ACOS(union_costheta.x[ii]);
        }
    }
    return union_returnvalue.m;
  }

#endif

  extern const int64_t bits_set_in_avx512_mask_float[];
  extern const uint16_t masks_per_misalignment_value_float[];

  extern const int64_t bits_set_in_avx512_mask_double[];
  extern const uint8_t masks_per_misalignment_value_double[];

#ifdef __cplusplus
}
#endif
