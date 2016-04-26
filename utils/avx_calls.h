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

#define AVX_NVEC                         8    
#define AVX_INTS                         __m256i
#define AVX_FLOATS                       __m256

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

#define AVX_BROADCAST_FLOAT(X)           _mm256_broadcast_ss(X);
#define AVX_SET_FLOAT(X)                 _mm256_set1_ps(X);


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
#define AVX_ARC_COSINE(X)                 _mm256_acos_ps(X)
#else
    //Other compilers do not have the vectorized arc-cosine
#define AVX_ARC_COSINE(X)                  inv_cosine(X)
#endif

    //Max
#define AVX_MAX_FLOATS(X,Y)               _mm256_max_ps(X,Y)

    //Casting (does not actual convert between types)
#define AVX_CAST_FLOAT_TO_INT(X)          _mm256_castps_si256(X)
#define AVX_CAST_INT_TO_FLOAT(X)          _mm256_castsi256_ps(X)

    //Streaming store
#define AVX_STREAMING_STORE_FLOATS(X,Y)   _mm256_stream_ps(X,Y)
#define AVX_STREAMING_STORE_INTS(X,Y)     _mm256_stream_si256(X,Y)

#else //DOUBLE PRECISION CALCULATIONS
#define AVX_NVEC                         4    
#define AVX_INTS                         __m128i
#define AVX_FLOATS                       __m256d
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

#define AVX_BROADCAST_FLOAT(X)            _mm256_broadcast_sd(X);
#define AVX_SET_FLOAT(X)                  _mm256_set1_pd(X);
//MoveMask
#define AVX_TEST_COMPARISON(X)            _mm256_movemask_pd(X)

#define AVX_BLEND_FLOATS_WITH_MASK(FALSEVALUE,TRUEVALUE,MASK) _mm256_blendv_pd(FALSEVALUE,TRUEVALUE,MASK)
#define AVX_MASKSTORE_FLOATS(dest, mask, source)   _mm256_maskstore_pd(dest, mask, source)

    //Trig
#ifdef  __INTEL_COMPILER
#define AVX_ARC_COSINE(X)                 _mm256_acos_pd(X)
#else
#define AVX_ARC_COSINE(X)                  inv_cosine(X)
#endif

    //Max
#define AVX_MAX_FLOATS(X,Y)               _mm256_max_pd(X,Y)

    //Casting (does not actual convert between types)
#define AVX_CAST_FLOAT_TO_INT(X)          _mm256_castpd_si256(X)
#define AVX_CAST_INT_TO_FLOAT(X)          _mm256_castsi256_pd(X)

    //Streaming store
#define AVX_STREAMING_STORE_FLOATS(X,Y)   _mm256_stream_pd(X,Y)
#define AVX_STREAMING_STORE_INTS(X,Y)     _mm_stream_si128(X,Y)

#endif


#ifndef  __INTEL_COMPILER
    static inline AVX_FLOATS inv_cosine(const AVX_FLOATS X)
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
        const DOUBLE zero = (DOUBLE) 0.0;

        for(int ii=0;ii<AVX_NVEC;ii++) {
            const DOUBLE costheta = union_costheta.x[ii];
            if(costheta < minus_one) {
                union_returnvalue.x[ii] = M_PI;
            } else if (costheta > one) {
                union_returnvalue.x[ii] = zero;
            } else {
                union_returnvalue.x[ii] = ACOS(costheta);
            }
        }
        return union_returnvalue.m;
    }
#endif

#ifdef __cplusplus
}
#endif
