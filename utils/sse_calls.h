/* File: sse_calls.h */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <nmmintrin.h>

#ifdef __cplusplus
extern "C" {
#endif

#include "function_precision.h"

#if defined(__GNUC__) || defined(__GNUG__)
#define SSE_BIT_COUNT_INT(X)                __builtin_popcount(X)
#else
#define SSE_BIT_COUNT_INT(X)                _popcnt32(X)
#endif


#ifndef DOUBLE_PREC

#define SSE_NVEC                         4
#define SSE_INTS                         __m128i
#define SSE_FLOATS                       __m128

#define SSE_SETZERO_FLOAT()              _mm_setzero_ps()    
#define SSE_LOAD_FLOATS_UNALIGNED(X)     _mm_loadu_ps(X)
#define SSE_MULTIPLY_FLOATS(X,Y)         _mm_mul_ps(X,Y)
#define SSE_SUBTRACT_FLOATS(X,Y)         _mm_sub_ps(X,Y)
#define SSE_ADD_FLOATS(X,Y)              _mm_add_ps(X,Y)
#define SSE_DIVIDE_FLOATS(X,Y)           _mm_div_ps(X,Y)
#define SSE_SQRT_FLOAT(X)                _mm_sqrt_ps(X)
#define SSE_TRUNCATE_FLOAT_TO_INT(X)     _mm_cvttps_epi32(X)
#define SSE_SQUARE_FLOAT(X)              _mm_mul_ps(X,X)
#define SSE_SET_FLOAT(X)                 _mm_set1_ps(X)

#define SSE_COMPARE_FLOATS_GE(X,Y)       _mm_cmpge_ps(X,Y)
#define SSE_COMPARE_FLOATS_LT(X,Y)       _mm_cmplt_ps(X,Y)
#define SSE_COMPARE_FLOATS_LE(X,Y)       _mm_cmple_ps(X,Y)    
#define SSE_COMPARE_FLOATS_GT(X,Y)       _mm_cmpgt_ps(X,Y)    
// X OP Y
//#define SSE_COMPARE_FLOATS(X,Y,OP)        _mm_cmp_ps(X,Y,OP)
#define SSE_BITWISE_AND(X,Y)              _mm_and_ps(X,Y)

//MoveMask
#define SSE_TEST_COMPARISON(X)            _mm_movemask_ps(X)

#define SSE_BLEND_FLOATS_WITH_MASK(FALSEVALUE,TRUEVALUE,MASK) _mm_blendv_ps(FALSEVALUE,TRUEVALUE,MASK)

//Max
#define SSE_MAX_FLOATS(X,Y)               _mm_max_ps(X,Y)

#define SSE_ABS_FLOAT(X)                  _mm_max_ps(_mm_sub_ps(_mm_setzero_ps(), X), X)
  
  
#ifdef  __INTEL_COMPILER
#define SSE_ARC_COSINE(X, order)                 _mm_acos_ps(X)
#else
#define SSE_ARC_COSINE(X, order)                  inv_cosine_sse(X, order)
#endif


#else //DOUBLE PRECISION CALCULATIONS

#define SSE_NVEC                         2
#define SSE_INTS                         __m128i
#define SSE_FLOATS                       __m128d

#define SSE_SETZERO_FLOAT()              _mm_setzero_pd()        
#define SSE_SET_FLOAT(X)                 _mm_set1_pd(X)
#define SSE_LOAD_FLOATS_UNALIGNED(X)     _mm_loadu_pd(X)
#define SSE_LOAD_FLOATS_ALIGNED(X)       _mm_load_pd(X)

//Math ops
#define SSE_MULTIPLY_FLOATS(X,Y)         _mm_mul_pd(X,Y)
#define SSE_SUBTRACT_FLOATS(X,Y)         _mm_sub_pd(X,Y)
#define SSE_ADD_FLOATS(X,Y)              _mm_add_pd(X,Y)
#define SSE_DIVIDE_FLOATS(X,Y)           _mm_div_pd(X,Y)
#define SSE_SQRT_FLOAT(X)                _mm_sqrt_pd(X)
#define SSE_SQUARE_FLOAT(X)              _mm_mul_pd(X,X)

//Memory stores
#define SSE_TRUNCATE_FLOAT_TO_INT(X)     _mm_cvttpd_epi32(X)
#define SSE_STORE_FLOATS_TO_MEMORY(X,Y)  _mm_storeu_pd(X,Y)

//The comparisons
#define SSE_COMPARE_FLOATS_GE(X,Y)       _mm_cmpge_pd(X,Y)
#define SSE_COMPARE_FLOATS_LT(X,Y)       _mm_cmplt_pd(X,Y)
#define SSE_COMPARE_FLOATS_LE(X,Y)       _mm_cmple_pd(X,Y)    
#define SSE_COMPARE_FLOATS_GT(X,Y)       _mm_cmpgt_pd(X,Y)    


//Bitwise and
#define SSE_BITWISE_AND(X,Y)              _mm_and_pd(X,Y)

//MoveMask
#define SSE_TEST_COMPARISON(X)            _mm_movemask_pd(X)

#define SSE_BLEND_FLOATS_WITH_MASK(FALSEVALUE,TRUEVALUE,MASK) _mm_blendv_pd(FALSEVALUE,TRUEVALUE,MASK)

#ifdef  __INTEL_COMPILER
#define SSE_ARC_COSINE(X, order)                 _mm_acos_pd(X)
#else
#define SSE_ARC_COSINE(X, order)                  inv_cosine_sse(X, order)
#endif

#define SSE_MAX_FLOATS(X,Y)               _mm_max_pd(X,Y)
#define SSE_ABS_FLOAT(X)                  _mm_max_pd(_mm_sub_pd(_mm_setzero_pd(), X), X)

#endif

#ifndef  __INTEL_COMPILER
#include "fast_acos.h"
    static inline SSE_FLOATS inv_cosine_sse(const SSE_FLOATS X, const int order)
{
    union cos{
        SSE_FLOATS m;
        DOUBLE x[SSE_NVEC];
    };
    union cos union_costheta;
    union cos union_returnvalue;
    union_costheta.m = X;
    const DOUBLE minus_one = (DOUBLE) -1.0;
    const DOUBLE one = (DOUBLE) 1.0;

    //Force everything to be in range [0,1]
    for(int ii=0;ii<SSE_NVEC;ii++) {
        const DOUBLE costheta = union_costheta.x[ii];
        union_costheta.x[ii] = costheta <= minus_one ? minus_one:costheta;
        union_costheta.x[ii] = costheta >= one ? one:costheta;
    }
    
    if(order==0) {
        for(int ii=0;ii<SSE_NVEC;ii++) {
            const DOUBLE costheta = union_costheta.x[ii];
            union_returnvalue.x[ii] = ACOS(costheta);
        }
    } else {
        //fast acos
        /*Taken from associated C++ code in http://www.geometrictools.com/GTEngine/Include/Mathematics/GteACosEstimate.h*/
        for(int ii=0;ii<SSE_NVEC;ii++) {
            union_returnvalue.x[ii] = FAST_ACOS(union_costheta.x[ii]);
        }
    }
    return union_returnvalue.m;
  }
#endif

  union int4 {
    SSE_INTS m_ibin;
    int ibin[SSE_NVEC];
  };
  
  union float4{
    SSE_FLOATS m_Dperp;
    DOUBLE Dperp[SSE_NVEC];
  };
  
  union float4_weights{
    SSE_FLOATS m_weights;
    DOUBLE weights[SSE_NVEC];
  };

#ifdef __cplusplus
}
#endif
    
