/* File: avx512_calls.h */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#pragma once

#if defined(__AVX512F__)

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
#define AVX512_MASK_BITWISE_NOT(X,Y)              _mm512_knot(X) //~X

    

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

#define AVX512_INTS                                      __m512i
#define AVX512_SET_INT(X)                                _mm512_set1_epi32(X)
#define AVX512_SETZERO_INT()                             _mm512_setzero_epi32()

        
#if 0
/* commenting out the integer math operations since they are either cumbersome or produce results of different SIMD widths*/    
#define AVX512_ADD_INTS(X, Y)                            _mm512_add_epi32(X, Y) 
#define AVX512_MASK_ADD_INTS(FALSEVALS, MASK, X, Y)      _mm512_mask_add_epi32(FALSEVALS, MASK, X, Y)
#define AVX512_MASKZ_ADD_INTS(MASK, X, Y)                _mm512_maskz_add_epi32(MASK, X, Y)

#define AVX512_MULTIPLY_INTS(X, Y)                                 _mm512_mul_epi32(X, Y) 
#define AVX512_MASK_MULTIPLY_INTS(FALSEVALS, MASK, X, Y)           _mm512_mask_mul_epi32(FALSEVALS, MASK, X, Y)
#define AVX512_MASKZ_MULTIPLY_INTS(MASK, X, Y)                     _mm512_maskz_mul_epi32(MASK, X, Y)
#define AVX512_MASKZ_MULTIPLY_INTS_LOW32(MASK, X, Y)               _mm512_maskz_mullo_epi32(MASK, X, Y)
#endif  /*end of integer math*/

#define AVX512_SETZERO_FLOAT()                                    _mm512_setzero_ps()
    
#define AVX512_LOAD_FLOATS_UNALIGNED(X)                            _mm512_loadu_ps(X)
#define AVX512_MASK_LOAD_FLOATS_UNALIGNED(FALSEVALS, MASK, X)      _mm512_mask_loadu_ps(FALSEVALS, MASK, X)
#define AVX512_MASKZ_LOAD_FLOATS_UNALIGNED(MASK, X)                _mm512_maskz_loadu_ps(MASK, X)

#define AVX512_LOAD_FLOATS_ALIGNED(X)                             _mm512_load_ps(X)
#define AVX512_MASK_LOAD_FLOATS_ALIGNED(FALSEVALS, MASK, X)       _mm512_mask_load_ps(FALSEVALS, MASK, X)
#define AVX512_MASKZ_LOAD_FLOATS_ALIGNED(MASK, X)                 _mm512_maskz_load_ps(MASK, X)


#define AVX512_MULTIPLY_FLOATS(X,Y)                               _mm512_mul_ps(X,Y)
#define AVX512_MASK_MULTIPLY_FLOATS(FALSEVALS, MASK, X,Y)         _mm512_mask_mul_ps(FALSEVALS, MASK, X,Y)
#define AVX512_MASKZ_MULTIPLY_FLOATS(MASK, X,Y)                   _mm512_maskz_mul_ps(MASK, X,Y)


#define AVX512_DIVIDE_FLOATS(X,Y)                                 _mm512_div_ps(X,Y)
#define AVX512_MASK_DIVIDE_FLOATS(FALSEVALS, MASK, X,Y)           _mm512_mask_div_ps(FALSEVALS, MASK, X,Y)
#define AVX512_MASKZ_DIVIDE_FLOATS(MASK, X,Y)                     _mm512_maskz_div_ps(MASK, X,Y)

#define AVX512_ADD_FLOATS(X,Y)                                    _mm512_add_ps(X,Y)
#define AVX512_MASK_ADD_FLOATS(FALSEVALS, MASK, X,Y)              _mm512_mask_add_ps(FALSEVALS, MASK, X,Y)
#define AVX512_MASKZ_ADD_FLOATS(MASK, X,Y)                        _mm512_maskz_add_ps(MASK, X,Y)

#define AVX512_SUBTRACT_FLOATS(X,Y)                               _mm512_sub_ps(X,Y)
#define AVX512_MASK_SUBTRACT_FLOATS(FALSEVALS, MASK, X,Y)         _mm512_mask_sub_ps(FALSEVALS, MASK, X,Y)
#define AVX512_MASKZ_SUBTRACT_FLOATS(MASK, X,Y)                   _mm512_maskz_sub_ps(MASK, X,Y)

  /* returns Z + XY*/
#define AVX512_FMA_ADD_FLOATS(X,Y,Z)                              _mm512_fmadd_ps(X,Y,Z)
#define AVX512_MASK_FMA_ADD_FLOATS(X, MASK, Y, Z)                 _mm512_mask_fmadd_ps(X, MASK, Y, Z)
#define AVX512_MASKZ_FMA_ADD_FLOATS(X, MASK, Y, Z)                _mm512_maskz_fmadd_ps(MASK, X, Y, Z)

  /* returns Z - XY*/
#define AVX512_FNMA_ADD_FLOATS(X, Y, Z)                           _mm512_fnmadd_ps(X, Y, Z) 
#define AVX512_MASK_FNMA_ADD_FLOATS(X, MASK, Y, Z)                _mm512_mask_fnmadd_ps(X, MASK, Y, Z) 
#define AVX512_MASKZ_FNMA_ADD_FLOATS(X, MASK, Y, Z)               _mm512_maskz_fnmadd_ps(MASK, X, Y, Z) 

  /* returns XY - Z */
#define AVX512_FMA_SUBTRACT_FLOATS(X,Y,Z)                         _mm512_fmsub_ps(X,Y,Z)
#define AVX512_MASK_FMA_SUBTRACT_FLOATS(X, MASK, Y, Z)            _mm512_mask_fmsub_ps(X, MASK, Y, Z)
#define AVX512_MASKZ_FMA_SUBTRACT_FLOATS(X, MASK, Y, Z)           _mm512_maskz_fmsub_ps(MASK, X, Y, Z)


#ifdef  __INTEL_COMPILER
#define AVX512_HORIZONTAL_SUM_FLOATS(X)                           _mm512_reduce_add_ps(X)
#define AVX512_MASK_HORIZONTAL_SUM_FLOATS(MASK, X)                _mm512_mask_reduce_add_ps(MASK, X)
#else
#define AVX512_HORIZONTAL_SUM_FLOATS(X)                           _horizontal_sum_floats(X)
#define AVX512_MASK_HORIZONTAL_SUM_FLOATS(MASK, X)                _horizontal_mask_sum_floats(X)
#endif

#define AVX512_SQRT_FLOAT(X)                                      _mm512_sqrt_ps(X)
#define AVX512_MASK_SQRT_FLOAT(FALSEVALS, MASK, X)                _mm512_mask_sqrt_ps(FALSEVALS, MASK, X)
#define AVX512_MASKZ_SQRT_FLOAT(MASK, X)                          _mm512_maskz_sqrt_ps(MASK, X)


#define AVX512_SVML_SQRT_FLOAT(X)                                 _mm512_svml_sqrt_ps(X)
#define AVX512_TRUNCATE_FLOAT_TO_INT(X)                           _mm512_cvttps_epi32(X)
#define AVX512_STORE_FLOATS_TO_MEMORY(X,Y)                        _mm512_storeu_ps(X,Y)
#define AVX512_SQUARE_FLOAT(X)                                    _mm512_mul_ps(X,X)
#define AVX512_MASKZ_SQUARE_FLOAT(MASK, X)                        _mm512_maskz_mul_ps(MASK, X,X)
#define AVX512_LOG_FLOAT(X)                                       _mm512_log_ps(X)
#define AVX512_LOG10_FLOAT(X)                                     _mm512_log10_ps(X)
#define AVX512_LOG2_FLOAT(X)                                      _mm512_log2_ps(X)

#define AVX512_RECIPROCAL_FLOATS(X)                               _mm512_rcp14_ps(X)
#define AVX512_MASK_RECIPROCAL_FLOATS(FALSEVALS, MASK, X)         _mm512_mask_rcp14_ps(FALSEVALS, MASK, X)
#define AVX512_MASKZ_RECIPROCAL_FLOATS(MASK, X)                   _mm512_maskz_rcp14_ps(MASK, X)

#define AVX512_SET_FLOAT(X)                                       _mm512_set1_ps(X)

// X OP Y
#define AVX512_COMPARE_FLOATS(X, Y, OP)                           _mm512_cmp_ps_mask(X, Y, OP)

//Mask operations (new in AVX512)
#define AVX512_MASK_COMPARE_FLOATS(M, X, Y, OP)                   _mm512_mask_cmp_ps_mask(M, X, Y, OP)
#define AVX512_BLEND_FLOATS_WITH_MASK(MASK, FALSE,TRUE)           _mm512_mask_blend_ps(MASK, FALSE,TRUE)
#define AVX512_BLEND_INTS_WITH_MASK(MASK, FALSE, TRUE)            _mm512_mask_blend_epi32(MASK, FALSE, TRUE)

//Trig
#ifdef  __INTEL_COMPILER
/* Needs SVML */
#define AVX512_ARC_COSINE(X, order)                 _mm512_acos_ps(X)  
#else
//Other compilers do not have the vectorized arc-cosine
#define AVX512_ARC_COSINE(X, order)                 inv_cosine_avx512(X, order)
#endif

//Max
#define AVX512_MAX_FLOATS(X,Y)               _mm512_max_ps(X,Y)


//Absolute value
#define AVX512_ABS_FLOAT(X)                  _mm512_abs_ps(X)
  
 //Casting (does not actual convert between types)
#define AVX512_CAST_FLOAT_TO_INT(X)          _mm512_castps_si512(X)
#define AVX512_CAST_INT_TO_FLOAT(X)          _mm512_castsi512_ps(X)


#else //DOUBLE PRECISION CALCULATIONS
  
#define DOUBLE                              double
#define AVX512_NVEC                         8    
#define AVX512_MASK                         __mmask8
#define AVX512_FLOATS                       __m512d

//This is AVX2 and not AVX512F 
#define AVX512_INTS                         __m256i

#define AVX512_SET_INT(X)                   _mm256_set1_epi32(X)
#define AVX512_SETZERO_INT()                _mm256_setzero_si256()
    
#if 0
#define AVX512_ADD_INTS(X, Y)               _mm256_add_epi32(X, Y) 
#define AVX512_MULTIPLY_INTS(X, Y)          _mm256_mul_epi32(X, Y) 

#if defined(__AVX512VL__)
#define AVX512_MASK_ADD_INTS(FALSEVALS, MASK, X, Y)               _mm256_mask_add_epi32(FALSEVALS, MASK, X, Y)
#define AVX512_MASKZ_ADD_INTS(MASK, X, Y)                         _mm256_maskz_add_epi32(MASK, X, Y)

#define AVX512_MASK_MULTIPLY_INTS(FALSEVALS, MASK, X, Y)          _mm256_mask_mul_epi32(FALSEVALS, MASK, X, Y)
#define AVX512_MASKZ_MULTIPLY_INTS(MASK, X, Y)                    _mm256_maskz_mul_epi32(MASK, X, Y)
#define AVX512_MASKZ_MULTIPLY_INTS_LOW32(MASK, X, Y)              _mm256_maskz_mullo_epi32(MASK, X, Y)
        
#elif defined(__AVX512F__)
#define AVX512_MASK_ADD_INTS(FALSEVALS, MASK, X, Y)               _mm512_castsi512_si256(_mm512_mask_add_epi32(_mm512_castsi256_si512(FALSEVALS), MASK, _mm512_castsi256_si512(X), _mm512_castsi256_si512(Y)))
#define AVX512_MASKZ_ADD_INTS(MASK, X, Y)                         _mm512_castsi512_si256(_mm512_maskz_add_epi32(MASK, _mm512_castsi256_si512(X), _mm512_castsi256_si512(Y)))

#define AVX512_MASK_MULTIPLY_INTS(FALSEVALS, MASK, X, Y)          _mm512_castsi512_si256(_mm512_mask_mul_epi32(_mm512_castsi256_si512(FALSEVALS), MASK, _mm512_castsi256_si512(X), _mm512_castsi256_si512(Y)))
#define AVX512_MASKZ_MULTIPLY_INTS(MASK, X, Y)                    _mm512_castsi512_si256(_mm512_maskz_mul_epi32(MASK, _mm512_castsi256_si512(X), _mm512_castsi256_si512(Y)))

#define AVX512_MASKZ_MULTIPLY_INTS_LOW32(MASK, X, Y)              _mm512_castsi512_si256(_mm512_maskz_mullo_epi32(MASK, _mm512_castsi256_si512(X), _mm512_castsi256_si512(Y)))

#endif
#endif /* commenting out the int math operations since they are either cumbersome or produce results of different SIMD widths*/

#define AVX512_SETZERO_FLOAT()                                    _mm512_setzero_pd()
    
#define AVX512_LOAD_FLOATS_UNALIGNED(X)                           _mm512_loadu_pd(X)
#define AVX512_MASK_LOAD_FLOATS_UNALIGNED(FALSEVALS, MASK, X)     _mm512_mask_loadu_pd(FALSEVALS, MASK, X)
#define AVX512_MASKZ_LOAD_FLOATS_UNALIGNED(MASK, X)               _mm512_maskz_loadu_pd(MASK, X)

#define AVX512_LOAD_FLOATS_ALIGNED(X)                             _mm512_load_pd(X)
#define AVX512_MASK_LOAD_FLOATS_ALIGNED(FALSEVALS, MASK, X)       _mm512_mask_load_pd(FALSEVALS, MASK, X)
#define AVX512_MASKZ_LOAD_FLOATS_ALIGNED(MASK, X)                 _mm512_maskz_load_pd(MASK, X)

#define AVX512_MULTIPLY_FLOATS(X,Y)                               _mm512_mul_pd(X,Y)
#define AVX512_MASK_MULTIPLY_FLOATS(FALSEVALS, MASK, X,Y)         _mm512_mask_mul_pd(FALSEVALS, MASK, X,Y)
#define AVX512_MASKZ_MULTIPLY_FLOATS(MASK, X,Y)                   _mm512_maskz_mul_pd(MASK, X,Y)

#define AVX512_DIVIDE_FLOATS(X,Y)                                 _mm512_div_pd(X,Y)
#define AVX512_MASK_DIVIDE_FLOATS(FALSEVALS, MASK, X,Y)           _mm512_mask_div_pd(FALSEVALS, MASK, X,Y)
#define AVX512_MASKZ_DIVIDE_FLOATS(MASK, X,Y)                     _mm512_maskz_div_pd(MASK, X,Y)

#define AVX512_ADD_FLOATS(X,Y)                                    _mm512_add_pd(X,Y)
#define AVX512_MASK_ADD_FLOATS(FALSEVALS, MASK, X,Y)              _mm512_mask_add_pd(FALSEVALS, MASK, X,Y)
#define AVX512_MASKZ_ADD_FLOATS(MASK, X,Y)                        _mm512_maskz_add_pd(MASK, X,Y)

#define AVX512_SUBTRACT_FLOATS(X,Y)                               _mm512_sub_pd(X,Y)
#define AVX512_MASK_SUBTRACT_FLOATS(FALSEVALS, MASK, X,Y)         _mm512_mask_sub_pd(FALSEVALS, MASK, X,Y)
#define AVX512_MASKZ_SUBTRACT_FLOATS(MASK, X,Y)                   _mm512_maskz_sub_pd(MASK, X,Y)



/* returns Z + XY*/
#define AVX512_FMA_ADD_FLOATS(X,Y,Z)                              _mm512_fmadd_pd(X,Y,Z)
#define AVX512_MASK_FMA_ADD_FLOATS(X, MASK, Y, Z)                 _mm512_mask_fmadd_pd(X, MASK, Y, Z)
#define AVX512_MASKZ_FMA_ADD_FLOATS(X, MASK, Y, Z)                _mm512_maskz_fmadd_pd(MASK, X, Y, Z)

/* returns Z - XY*/
#define AVX512_FNMA_ADD_FLOATS(X, Y, Z)                           _mm512_fnmadd_pd(X, Y, Z) 
#define AVX512_MASK_FNMA_ADD_FLOATS(X, MASK, Y, Z)                _mm512_mask_fnmadd_pd(X, MASK, Y, Z) 
#define AVX512_MASKZ_FNMA_ADD_FLOATS(X, MASK, Y, Z)               _mm512_maskz_fnmadd_pd(MASK, X, Y, Z) 

/* returns XY - Z */
#define AVX512_FMA_SUBTRACT_FLOATS(X,Y,Z)                         _mm512_fmsub_pd(X,Y,Z)
#define AVX512_MASK_FMA_SUBTRACT_FLOATS(X, MASK, Y, Z)            _mm512_mask_fmsub_pd(X, MASK, Y, Z)
#define AVX512_MASKZ_FMA_SUBTRACT_FLOATS(X, MASK, Y, Z)           _mm512_maskz_fmsub_pd(MASK, X, Y, Z)


#ifdef  __INTEL_COMPILER
#define AVX512_HORIZONTAL_SUM_FLOATS(X)                           _mm512_reduce_add_pd(X)
#define AVX512_MASK_HORIZONTAL_SUM_FLOATS(MASK, X)                _mm512_mask_reduce_add_pd(MASK, X)
#endif

#define AVX512_SQRT_FLOAT(X)                                      _mm512_sqrt_pd(X)
#define AVX512_MASK_SQRT_FLOAT(FALSEVALS, MASK, X)                _mm512_mask_sqrt_pd(FALSEVALS, MASK, X)
#define AVX512_MASKZ_SQRT_FLOAT(MASK, X)                          _mm512_maskz_sqrt_pd(MASK, X)


#define AVX512_SVML_SQRT_FLOAT(X)                                 _mm512_svml_sqrt_pd(X)
#define AVX512_TRUNCATE_FLOAT_TO_INT(X)                           _mm512_cvttpd_epi32(X) 
#define AVX512_STORE_FLOATS_TO_MEMORY(X,Y)                        _mm512_storeu_pd(X,Y)
#define AVX512_SQUARE_FLOAT(X)                                    _mm512_mul_pd(X,X)
#define AVX512_MASKZ_SQUARE_FLOAT(MASK, X)                        _mm512_maskz_mul_pd(MASK, X,X)
#define AVX512_LOG_FLOAT(X)                                       _mm512_log_pd(X)
#define AVX512_LOG2_FLOAT(X)                                      _mm512_log2_pd(X)
#define AVX512_LOG10_FLOAT(X)                                     _mm512_log10_pd(X)

#define AVX512_RECIPROCAL_FLOATS(X)                               _mm512_rcp14_pd(X)
#define AVX512_MASK_RECIPROCAL_FLOATS(FALSEVALS, MASK, X)         _mm512_mask_rcp14_pd(FALSEVALS, MASK, X)
#define AVX512_MASKZ_RECIPROCAL_FLOATS(MASK, X)                   _mm512_maskz_rcp14_pd(MASK, X)


    // X OP Y
#define AVX512_COMPARE_FLOATS(X, Y, OP)                           _mm512_cmp_pd_mask(X, Y, OP)
#define AVX512_MASK_COMPARE_FLOATS(M, X, Y, OP)                   _mm512_mask_cmp_pd_mask(M, X, Y, OP)

#define AVX512_SET_FLOAT(X)                                       _mm512_set1_pd(X)

#define AVX512_BLEND_FLOATS_WITH_MASK(MASK, FALSEVALUE,TRUEVALUE) _mm512_mask_blend_pd(MASK, FALSEVALUE,TRUEVALUE)

#if defined(__AVX512VL__)
#define AVX512_BLEND_INTS_WITH_MASK(MASK, FALSE,TRUE)             _mm256_mask_blend_epi32(MASK, FALSE,TRUE)
#elif defined(__AVX2__)
#define AVX512_BLEND_INTS_WITH_MASK(MASK, FALSE,TRUE)             _mm256_blend_epi32(FALSE, TRUE, MASK)//AVX2
#else
#define AVX512_BLEND_INTS_WITH_MASK(MASK, FALSE, TRUE)            _blend_epi32_with_mask(MASK, FALSE, TRUE)
#endif


//Trig
#ifdef  __INTEL_COMPILER
#define AVX512_ARC_COSINE(X, order)                 _mm512_acos_pd(X)
#else
#define AVX512_ARC_COSINE(X, order)                  inv_cosine_avx512(X, order)
#endif

    //Max
#define AVX512_MAX_FLOATS(X,Y)               _mm512_max_pd(X,Y)

//Absolute value
#if __GNUC__  <=  8 
    //there was a bug for the function proto-type
    //for _mm512_abs_pd -- see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=87467
    //Need to protect the user
#define AVX512_ABS_FLOAT(X)                  _mm512_max_pd(X,_mm512_sub_pd(_mm512_setzero_pd(), X))
#else
#define AVX512_ABS_FLOAT(X)                  _mm512_abs_pd(X)
#endif
  
 //Casting (does not actual convert between types)
#define AVX512_CAST_FLOAT_TO_INT(X)          _mm512_castpd_si512(X)
#define AVX512_CAST_INT_TO_FLOAT(X)          _mm512_castsi512_pd(_mm512_castsi256_si512(X))

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

  union int16 {
    AVX512_INTS m_ibin;
    int ibin[AVX512_NVEC];
  };
  union float16{
    AVX512_FLOATS m_Dperp;
    DOUBLE Dperp[AVX512_NVEC];
  };
  
  union float16_weights{
    AVX512_FLOATS m_weights;
    DOUBLE weights[AVX512_NVEC];
  };
  

    
#define CHECK_AND_FAST_DIVIDE_AVX512(result, numerator, denominator, mask, fast_divide_and_NR_steps) { \
      /* For double precision floats */                                 \
      if (fast_divide_and_NR_steps == 0) {                              \
          result = AVX512_MASKZ_DIVIDE_FLOATS(mask, numerator, denominator); \
      } else {                                                          \
          /* following blocks do an approximate reciprocal followed by two iterations of Newton-Raphson */ \
          const AVX512_FLOATS rc = AVX512_MASKZ_RECIPROCAL_FLOATS(mask, denominator); \
          /* We have the double->float->approx. reciprocal->double process done. */ \
          /* Now improve the accuracy of the divide with newton-raphson. */ \
          /* Ist iteration of NewtonRaphson */                          \
          const AVX512_FLOATS two = AVX512_SET_FLOAT((DOUBLE) 2.0);     \
          AVX512_FLOATS rc_iter = rc;                                   \
          for(unsigned int _ii=0;_ii<fast_divide_and_NR_steps;_ii++) {  \
              rc_iter = AVX512_MULTIPLY_FLOATS(rc_iter,                 \
                                               AVX512_FNMA_ADD_FLOATS(denominator, rc_iter, two));/*2.0 - l^2*rc */ \
          }                                                             \
          result = AVX512_MASKZ_MULTIPLY_FLOATS(mask, numerator, rc_iter); \
      } /* end of FAST_DIVIDE */                                        \
  }
    

#ifdef __cplusplus
}
#endif


#endif /* if defined(AVX512F) */
