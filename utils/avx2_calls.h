/* File: avx2_calls.h */
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
#define AVX2_BIT_COUNT_INT(X)                __builtin_popcount(X)
#else
#define AVX2_BIT_COUNT_INT(X)                _popcnt32(X)
#endif


#define AVX2_BIT_COUNT_LONG(X)               _popcnt64(X)
#define AVX2_BIT_COUNT_UNSIGNED_INT(X)       _mm_popcnt_u32(X)
#define AVX2_BIT_COUNT_UNSIGNED_LONG(X)      _mm_popcnt_u64(X)

#ifndef DOUBLE_PREC

#define DOUBLE                           float  
#define AVX2_NVEC                         8    
#define AVX2_INTS                         __m256i
#define AVX2_FLOATS                       __m256

#define AVX2_SETZERO_FLOAT()              _mm256_setzero_ps()
#define AVX2_LOAD_FLOATS_UNALIGNED(X)     _mm256_loadu_ps(X)
#define AVX2_LOAD_FLOATS_ALIGNED(X)       _mm256_load_ps(X)
#define AVX2_MULTIPLY_FLOATS(X,Y)         _mm256_mul_ps(X,Y)
#define AVX2_DIVIDE_FLOATS(X,Y)           _mm256_div_ps(X,Y)
#define AVX2_SUBTRACT_FLOATS(X,Y)         _mm256_sub_ps(X,Y)
#if defined(__FMA__) 
#define AVX2_FMA_ADD_FLOATS(X,Y,Z)        _mm256_fmadd_ps(X,Y,Z)
#elif defined(__FMA4__) /* on AMD cpus with FMA4 */
#define AVX2_FMA_ADD_FLOATS(X,Y,Z)        _mm256_macc_ps(X,Y,Z)
#else
#error Can not detect applicable FMA instruction 
#endif

#define AVX2_ADD_FLOATS(X,Y)              _mm256_add_ps(X,Y)

#define AVX2_SQRT_FLOAT(X)                _mm256_sqrt_ps(X)
#define AVX2_SVML_SQRT_FLOAT(X)           _mm256_svml_sqrt_ps(X)
#define AVX2_TRUNCATE_FLOAT_TO_INT(X)     _mm256_cvttps_epi32(X)
#define AVX2_STORE_FLOATS_TO_MEMORY(X,Y)  _mm256_storeu_ps(X,Y)
#define AVX2_SQUARE_FLOAT(X)              _mm256_mul_ps(X,X)
#define AVX2_LOG_FLOAT(X)                 _mm256_log_ps(X)
#define AVX2_LOG10_FLOAT(X)                 _mm256_log10_ps(X)
#define AVX2_LOG2_FLOAT(X)                _mm256_log2_ps(X)
#define AVX2_RECIPROCAL_FLOATS(X)         _mm256_rcp_ps(X)

#define AVX2_BROADCAST_FLOAT(X)           _mm256_broadcast_ss(X)
#define AVX2_SET_FLOAT(X)                 _mm256_set1_ps(X)
#define AVX2_SET_INT(X)                   _mm256_set1_epi32(X)

    // X OP Y
#define AVX2_COMPARE_FLOATS(X,Y,OP)        _mm256_cmp_ps(X,Y,OP)
#define AVX2_BITWISE_AND(X,Y)              _mm256_and_ps(X,Y)
#define AVX2_BITWISE_OR(X,Y)               _mm256_or_ps(X,Y)
#define AVX2_XOR_FLOATS(X,Y)               _mm256_xor_ps(X,Y)
#define AVX2_AND_NOT(X,Y)                  _mm256_andnot_ps(X,Y)  //~X & Y

//MoveMask
#define AVX2_TEST_COMPARISON(X)            _mm256_movemask_ps(X)

#define AVX2_BLEND_FLOATS_WITH_MASK(FALSEVALUE,TRUEVALUE,MASK) _mm256_blendv_ps(FALSEVALUE,TRUEVALUE,MASK)
#define AVX2_BLEND_INTS_WITH_MASK(FALSEVALUE, TRUEVALUE, MASK) _mm256_blend_epi32(FALSEVALUE,TRUEVALUE,MASK)

#define AVX2_MASKSTORE_FLOATS(dest, mask, source)   _mm256_maskstore_ps(dest, mask, source)

//Trig
#ifdef  __INTEL_COMPILER
#define AVX2_ARC_COSINE(X, order)                 _mm256_acos_ps(X)
#else
    //Other compilers do not have the vectorized arc-cosine
#define AVX2_ARC_COSINE(X, order)                  inv_cosine_avx2(X, order)
#endif

    //Max
#define AVX2_MAX_FLOATS(X,Y)               _mm256_max_ps(X,Y)


  //Absolute value
#define AVX2_ABS_FLOAT(X)                  _mm256_max_ps(_mm256_sub_ps(_mm256_setzero_ps(), X), X)
  
    //Casting (does not actual convert between types)
#define AVX2_CAST_FLOAT_TO_INT(X)          _mm256_castps_si256(X)
#define AVX2_CAST_INT_TO_FLOAT(X)          _mm256_castsi256_ps(X)

    //Streaming store
#define AVX2_STREAMING_STORE_FLOATS(X,Y)   _mm256_stream_ps(X,Y)
#define AVX2_STREAMING_STORE_INTS(X,Y)     _mm256_stream_si256(X,Y)

#else //DOUBLE PRECISION CALCULATIONS
  
#define DOUBLE                           double
#define AVX2_NVEC                         4    
#define AVX2_INTS                         __m128i
#define AVX2_FLOATS                       __m256d

#define AVX2_SETZERO_FLOAT()              _mm256_setzero_pd()    
#define AVX2_LOAD_FLOATS_UNALIGNED(X)     _mm256_loadu_pd(X)
#define AVX2_LOAD_FLOATS_ALIGNED(X)       _mm256_load_pd(X)
#define AVX2_MULTIPLY_FLOATS(X,Y)         _mm256_mul_pd(X,Y)
#define AVX2_DIVIDE_FLOATS(X,Y)           _mm256_div_pd(X,Y)
#define AVX2_SUBTRACT_FLOATS(X,Y)         _mm256_sub_pd(X,Y)
#define AVX2_ADD_FLOATS(X,Y)              _mm256_add_pd(X,Y)

#if defined(__FMA__) 
#define AVX2_FMA_ADD_FLOATS(X,Y,Z)        _mm256_fmadd_pd(X,Y,Z)
#elif defined(__FMA4__)
#define AVX2_FMA_ADD_FLOATS(X,Y,Z)        _mm256_macc_pd(X,Y,Z)
#else
#error Can not detect applicable FMA instruction 
#endif

#define AVX2_SQRT_FLOAT(X)                _mm256_sqrt_pd(X)
#define AVX2_SVML_SQRT_FLOAT(X)           _mm256_svml_sqrt_pd(X)
#define AVX2_TRUNCATE_FLOAT_TO_INT(X)     _mm256_cvttpd_epi32(X)
#define AVX2_STORE_FLOATS_TO_MEMORY(X,Y)  _mm256_storeu_pd(X,Y)
#define AVX2_SQUARE_FLOAT(X)              _mm256_mul_pd(X,X)
#define AVX2_LOG_FLOAT(X)                 _mm256_log_pd(X)
#define AVX2_LOG2_FLOAT(X)                _mm256_log2_pd(X)
#define AVX2_LOG10_FLOAT(X)                _mm256_log10_pd(X)
#define AVX2_RECIPROCAL_FLOATS(X)         _mm256_rcp_pd(X)

    // X OP Y
#define AVX2_COMPARE_FLOATS(X,Y,OP)        _mm256_cmp_pd(X,Y,OP)
#define AVX2_BITWISE_AND(X,Y)              _mm256_and_pd(X,Y)
#define AVX2_BITWISE_OR(X,Y)               _mm256_or_pd(X,Y)
#define AVX2_XOR_FLOATS(X,Y)               _mm256_xor_pd(X,Y)
#define AVX2_AND_NOT(X,Y)                  _mm256_andnot_pd((X),(Y))  //~X & Y

#define AVX2_BROADCAST_FLOAT(X)            _mm256_broadcast_sd(X)
#define AVX2_SET_FLOAT(X)                  _mm256_set1_pd(X)
#define AVX2_SET_INT(X)                    _mm_set1_epi32(X)

//MoveMask
#define AVX2_TEST_COMPARISON(X)            _mm256_movemask_pd(X)

#define AVX2_BLEND_FLOATS_WITH_MASK(FALSEVALUE,TRUEVALUE,MASK) _mm256_blendv_pd(FALSEVALUE,TRUEVALUE,MASK)
#define AVX2_BLEND_INTS_WITH_MASK(FALSEVALUE, TRUEVALUE, MASK) _mm_blend_epi32(FALSEVALUE,TRUEVALUE,MASK)

#define AVX2_MASKSTORE_FLOATS(dest, mask, source)   _mm256_maskstore_pd(dest, mask, source)

//Trig
#ifdef  __INTEL_COMPILER
#define AVX2_ARC_COSINE(X, order)                 _mm256_acos_pd(X)
#else
#define AVX2_ARC_COSINE(X, order)                  inv_cosine_avx(X, order) //avx and avx2 are the same wrt floats
#endif

//Max
#define AVX2_MAX_FLOATS(X,Y)               _mm256_max_pd(X,Y)

//Absolute value
#define AVX2_ABS_FLOAT(X)                  _mm256_max_pd(_mm256_sub_pd(_mm256_setzero_pd(), X), X)
  
    //Casting (does not actual convert between types)
#define AVX2_CAST_FLOAT_TO_INT(X)          _mm256_castpd_si256(X)
#define AVX2_CAST_INT_TO_FLOAT(X)          _mm256_castsi256_pd(X)

    //Streaming store
#define AVX2_STREAMING_STORE_FLOATS(X,Y)   _mm256_stream_pd(X,Y)
#define AVX2_STREAMING_STORE_INTS(X,Y)     _mm_stream_si128(X,Y)

#endif //DOUBLE_PREC

  //include all the avx matters including the declarations of union_int8 etc
#include "avx_calls.h"

#ifdef __cplusplus
}
#endif
