/* File: ARM64_calls.h */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#pragma once

// Taken from https://github.com/noloader/SHA-Intrinsics/blob/master/sha256-arm.c
// MS 11th Oct, 2022
#if defined(__arm__) || defined(__aarch32__) || defined(__arm64__) || defined(__aarch64__) || defined(_M_ARM)
# if defined(__GNUC__)
#  include <stdint.h>
# endif
# if defined(__ARM_NEON) || defined(_MSC_VER) || defined(__GNUC__)
#  include <arm_neon.h>
# endif
/* GCC and LLVM Clang, but not Apple Clang */
# if defined(__GNUC__) && !defined(__apple_build_version__)
#  if defined(__ARM_ACLE) || defined(__ARM_FEATURE_CRYPTO)
#   include <arm_acle.h>
#  endif
# endif
#endif  /* ARM Headers */

#ifdef __cplusplus
extern "C" {
#endif

#include "function_precision.h" 


#if defined(__GNUC__) || defined(__GNUG__)
#define ARM64_BIT_COUNT_INT(X)                __builtin_popcount(X)
#else
#define ARM64_BIT_COUNT_INT(X)                _popcnt32(X)
//#define ARM64_BIT_COUNT_INT(X)                   vaddv_u32(vcnt_u8(vreinterpret_u8_s32(X))
#endif


// #define ARM64_BIT_COUNT_LONG(X)               _popcnt64(X)
// #define ARM64_BIT_COUNT_UNSIGNED_INT(X)       _mm_popcnt_u32(X)
// #define ARM64_BIT_COUNT_UNSIGNED_LONG(X)      _mm_popcnt_u64(X)
// #define ARM64_SET_INT(X)                      _mm256_set1_epi32(X)


#ifndef DOUBLE_PREC

#define ARM64_NVEC                         4    
#define ARM64_INTS                         int32x4_t
#define ARM64_UINTS                        uint32x4_t
#define ARM64_FLOATS                       float32x4_t
#define ARM64_INT_TYPE                     int32_t
#define ARM64_UINT_TYPE                    uint32_t
#define ARM64_INT_PRINT_FORMAT             "d"
#define ARM64_LAST_LANE_IMM                3


// or should this be vmovq_n_f32?
#define ARM64_SETZERO_FLOAT()                    vdupq_n_f32(0)
#define ARM64_SET_FLOAT(X)                       vdupq_n_f32(X)
#define ARM64_SET_INT(X)                         vdupq_n_s32(X)
#define ARM64_SET_UINT(X)                        vdupq_n_u32(X)

    
#define ARM64_LOAD_FLOATS(X)                     vld1q_f32(X)
#define ARM64_MULTIPLY_FLOATS(X,Y)               vmulq_f32(X,Y)
#define ARM64_DIVIDE_FLOATS(X,Y)                 vdivq_f32(X,Y)
#define ARM64_SUBTRACT_FLOATS(X,Y)               vsubq_f32(X,Y)
#define ARM64_ADD_FLOATS(X,Y)                    vaddq_f32(X,Y)
#define ARM64_SQRT_FLOAT(X)                      vsqrtq_f32(X)
#define ARM64_TRUNCATE_FLOAT_TO_INT(X)           vcvtq_s32_f32(X)
#define ARM64_SQUARE_FLOAT(X)                    vmulq_f32(X,X)
#define ARM64_FMA_FLOATS(X,Y,Z)                  vfmaq_f32(X,Y,Z)
#define ARM64_RECIPROCAL_FLOATS(X)               vrecpeq_f32(X)

//Casting (does not actual convert between types)
#define ARM64_CAST_FLOAT_TO_INT(X)          vreinterpretq_s32_f32(X)
#define ARM64_CAST_INT_TO_FLOAT(X)          vreinterpretq_f32_s32(X)
#define ARM64_CAST_UINT_TO_FLOAT(X)         vreinterpretq_f32_u32(X)

#define ARM64_GET_LANE_FROM_INTS(X, lane)    vgetq_lane_s32(X, lane)
#define ARM64_GET_LANE_FROM_UINTS(X, lane)   vgetq_lane_u32(X, lane)

// Store to memory
#define ARM64_STORE_FLOATS_TO_MEMORY(ptr, X)     vst1q_f32(ptr,X)
#define ARM64_STORE_INTS_TO_MEMORY(ptr, X)       vst1q_s32(ptr,X)
#define ARM64_STORE_UINTS_TO_MEMORY(ptr, X)      vst1q_u32(ptr,X)


// Compare ops
#define ARM64_COMPARE_FLOATS_GE(X,Y)       vcgeq_f32(X,Y)
#define ARM64_COMPARE_FLOATS_LT(X,Y)       vcltq_f32(X,Y)
#define ARM64_COMPARE_FLOATS_LE(X,Y)       vcleq_f32(X,Y)
#define ARM64_COMPARE_FLOATS_GT(X,Y)       vcgtq_f32(X,Y)

// logical ops
#define ARM64_BITWISE_AND(X,Y)              vandq_u32(X,Y)
#define ARM64_BITWISE_OR(X,Y)               vorrq_u32(X,Y)
#define ARM64_XOR_FLOATS(X,Y)               veorq_u32(X,Y)

static inline int ARM64_COUNT_MATCHES(ARM64_UINTS input)
{
    //There are only 16 possibilities -> 0x0000, 0x000F, 0x000F0, 0x
    const uint32x4_t mask = { 1, 2, 4, 8 };
    // const uint8x16_t val = vandq_u32(input, mask);//extract bits and cast to u8
    return vaddvq_u32(vandq_u32(input, mask));

//   const uint32x4_t mask = { 1, 2, 4, 8 };
//   return vaddvq_u32(vandq_u32((uint32x4_t)val, mask));    
    // uint32x4_t sum = vpaddq_u32(val, val);
    // sum = vpaddq_u32(sum, sum);
    // tmp = vpaddq_u8(tmp,tmp);
    // return vgetq_lane_u32(sum, 0);
    // return vaddvq_u32(vandq_u32(input, mask));

    // const uint32_t one   = vgetq_lane_u32(input, 0) & 1u;
    // const uint32_t two   = vgetq_lane_u32(input, 1) & 1u;
    // const uint32_t three = vgetq_lane_u32(input, 2) & 1u;
    // const uint32_t four  = vgetq_lane_u32(input, 3) & 1u;
    // return one + two + three + four;
}


//MoveMask
static inline ARM64_UINT_TYPE ARM64_TEST_COMPARISON(ARM64_UINTS input)
{
    // Taken from sse2neon.h https://github.com/DLTcollab/sse2neon/blob/master/sse2neon.h
    //Can't do a set to quad register from 4 floats with one NEON op -> hence this has to be
    //a function (rather than a macro)
    //One solution would be to pass the temporary variable (for shift) into the macro itself
    static const int32x4_t shift = {0, 1, 2, 3};
    uint32x4_t tmp = vshrq_n_u32(input, 31);
    return vaddvq_u32(vshlq_u32(tmp, shift));

    //https://github.com/WebAssembly/simd/issues/131
    // const uint32x4_t mask = { 1, 2, 4, 8 };
    // return vaddvq_u32(vandq_u32(input, mask));   

    // const uint64_t magic = 0x103071;//correct magic number for 32 bit uint
    // uint64_t lo = ((vgetq_lane_u64(input, 0) * magic) >> 56);
    // uint64_t hi = ((vgetq_lane_u64(input, 1) * magic) >> 48) & 0xFF00;
    // return (hi + lo);

    //Taken from https://github.com/simdjson/simdjson/discussions/1658
    // const uint32x4_t magic = { 0xf0f1f3f8, 0x0f1f3f80, 0xf0f1f3f8, 0x0f1f3f80 };
    // const uint8x8_t idx = { 3, 7, 11, 15, 19, 23, 27, 31 };

    // uint8x16x2_t tbl = {
    //     vpaddq_u32(vmulq_u32(vgetq_lane_u32(input, 0), magic), vmulq_u32(vgetq_lane_u32(input, 1), magic)),
    //     vpaddq_u32(vmulq_u32(vgetq_lane_u32(input, 2), magic), vmulq_u32(vgetq_lane_u32(input, 3), magic)),
    //  };
    // return vget_lane_u64(vreinterpret_u64_u8(vqtbl2_u8(tbl, idx)), 0);

}
#define ARM64_BLEND_FLOATS_WITH_MASK(FALSEVALUE,TRUEVALUE,MASK)      vbslq_f32(MASK, TRUEVALUE, FALSEVALUE)
#define ARM64_BLEND_INTS_WITH_MASK(FALSEVALUE, TRUEVALUE, MASK)      vbslq_s32(MASK, TRUEVALUE, FALSEVALUE)
#define ARM64_BLEND_UINTS_WITH_MASK(FALSEVALUE, TRUEVALUE, MASK)     vbslq_u32(MASK, TRUEVALUE, FALSEVALUE)
//#define ARM64_BLEND_INTS_WITH_MASK(FALSEVALUE, TRUEVALUE, MASK)      vbslq_s32(vshrq_n_u32(MASK, 31), TRUEVALUE, FALSEVALUE)


//#define ARM64_MASKSTORE_FLOATS(dest, mask, source)                   _mm256_maskstore_ps(dest, mask, source)

//Trig
#define ARM64_ARC_COSINE(X, order)                  inv_cosine_arm64(X, order)

//Max
#define ARM64_MAX_FLOATS(X,Y)               vmaxq_f32(X,Y)


 //Absolute value
//#define ARM64_ABS_FLOAT(X)                  _mm256_max_ps(_mm256_sub_ps(_mm256_setzero_ps(), X), X)
#define ARM64_ABS_FLOAT(X)                  vabsq_f32(X)

#else //DOUBLE PRECISION CALCULATIONS
  
#define DOUBLE                            double  
#define ARM64_NVEC                         2    
#define ARM64_INTS                         int64x2_t
#define ARM64_UINTS                        uint64x2_t
#define ARM64_FLOATS                       float64x2_t
#define ARM64_INT_TYPE                     int64_t
#define ARM64_UINT_TYPE                    uint64_t
#define ARM64_INT_PRINT_FORMAT             "lld"
#define ARM64_LAST_LANE_IMM                1

// or should this be vmovq_n_f32?
#define ARM64_SETZERO_FLOAT()                    vdupq_n_f64(0)
#define ARM64_SET_FLOAT(X)                       vdupq_n_f64(X)
#define ARM64_SET_INT(X)                         vdupq_n_s64(X)
#define ARM64_SET_UINT(X)                        vdupq_n_u64(X)

    
#define ARM64_LOAD_FLOATS(X)                     vld1q_f64(X)
#define ARM64_MULTIPLY_FLOATS(X,Y)               vmulq_f64(X,Y)
#define ARM64_DIVIDE_FLOATS(X,Y)                 vdivq_f64(X,Y)
#define ARM64_SUBTRACT_FLOATS(X,Y)               vsubq_f64(X,Y)
#define ARM64_ADD_FLOATS(X,Y)                    vaddq_f64(X,Y)
#define ARM64_SQRT_FLOAT(X)                      vsqrtq_f64(X)
#define ARM64_TRUNCATE_FLOAT_TO_INT(X)           vcvtq_s64_f64(X)
#define ARM64_SQUARE_FLOAT(X)                    vmulq_f64(X,X)
#define ARM64_FMA_FLOATS(X,Y,Z)                  vfmaq_f64(X,Y,Z) // x + y*z
#define ARM64_RECIPROCAL_FLOATS(X)               vrecpeq_f64(X)

//Casting (does not actual convert between types)
#define ARM64_CAST_FLOAT_TO_INT(X)          vreinterpretq_s64_f64(X)
#define ARM64_CAST_FLOAT_TO_UINT(X)         vreinterpretq_u64_f64(X)
#define ARM64_CAST_INT_TO_FLOAT(X)          vreinterpretq_f64_s64(X)
#define ARM64_CAST_UINT_TO_FLOAT(X)         vreinterpretq_f64_u64(X)

//Get individual elements out from a packed integer
#define ARM64_GET_LANE_FROM_INTS(X, lane)    vgetq_lane_s64(X, lane)
#define ARM64_GET_LANE_FROM_UINTS(X, lane)   vgetq_lane_u64(X, lane)


#define ARM64_STORE_FLOATS_TO_MEMORY(ptr, X)     vst1q_f64(ptr,X)
#define ARM64_STORE_INTS_TO_MEMORY(ptr, X)       vst1q_s64(ptr,X)
#define ARM64_STORE_UINTS_TO_MEMORY(ptr, X)      vst1q_u64(ptr,X)

// Compare operations
#define ARM64_COMPARE_FLOATS_GE(X,Y)       vcgeq_f64(X,Y)
#define ARM64_COMPARE_FLOATS_LT(X,Y)       vcltq_f64(X,Y)
#define ARM64_COMPARE_FLOATS_LE(X,Y)       vcleq_f64(X,Y)    
#define ARM64_COMPARE_FLOATS_GT(X,Y)       vcgtq_f64(X,Y)    

// Logical ops
#define ARM64_BITWISE_AND(X,Y)              vandq_u64(X,Y)
#define ARM64_BITWISE_OR(X,Y)               vorrq_u64(X,Y)
//#define ARM64_XOR_FLOATS(X,Y)               veorq_u64(X,Y)

//my version -> just do the bit-count
static inline int ARM64_COUNT_MATCHES(ARM64_UINTS input)
{
    //There are only four possibilities -> 0x00, 0x0F, 0xF0, 0xFF
    // const uint64_t nmatches = {0, 1, 1, 2};
    // const uint64x2_t mask = {1, 1};
    // return vaddvq_u64(vandq_u64(input, mask));

    const uint32_t lo = vgetq_lane_u64(input, 0) & 1u;
    const uint32_t hi = vgetq_lane_u64(input, 1) & 1u;
    return hi + lo;
}

//MoveMask
// Taken from sse2neon.h https://github.com/DLTcollab/sse2neon/blob/master/sse2neon.h
static inline ARM64_INT_TYPE ARM64_TEST_COMPARISON(ARM64_UINTS input)
{
    uint64x2_t high_bits = vshrq_n_u64(input, 63);
    return vgetq_lane_u64(high_bits, 0) | (vgetq_lane_u64(high_bits, 1) << 1);

    // const uint8x8_t res   = vshrn_n_u16(vreinterpretq_u16_u64(input), 4);
    // const uint64_t matches = vget_lane_u64(vreinterpret_u64_u8(res), 0);
    // // return matches & 0x8000000080000000ull;
    // return matches;

    //https://github.com/WebAssembly/simd/issues/131
    // const uint64x2_t mask = { 1, 2 };
    // return vaddvq_u64(vandq_u64(input, mask));

    // how about a table lookup -> this can only be 0x00, 0x0F, 0xF0 or 0xFF?

    //https://twitter.com/0b0000000000000/status/1376568414634840065
    // const uint64_t magic = 0x103070F1F3F80ULL;
    // uint64_t lo = ((vgetq_lane_u64(input, 0) * magic) >> 56);
    // uint64_t hi = ((vgetq_lane_u64(input, 1) * magic) >> 48) & 0xFF00;
    // return (hi + lo);
    
    //Again from twitter: https://twitter.com/0b0000000000000/status/1459018074233847813?s=20&t=eZhtdP4hFxVd6xIutO9fPg
    //0x103070F1F3F80ULL (decimal representation 284803830071168)
    // const uint64_t magic = 0x103070F1F3F80ULL;
    // uint64_t lo = ((vgetq_lane_u64(input, 0) * magic) >> 56);
    // uint64_t hi = ((vgetq_lane_u64(input, 1) * magic) >> 48) & 0xFF00;
    // return (hi + lo);
}

#define ARM64_BLEND_FLOATS_WITH_MASK(FALSEVALUE, TRUEVALUE, MASK)      vbslq_f64(MASK, TRUEVALUE, FALSEVALUE)
#define ARM64_BLEND_INTS_WITH_MASK(FALSEVALUE, TRUEVALUE, MASK)        vbslq_s64(MASK, TRUEVALUE, FALSEVALUE)
#define ARM64_BLEND_UINTS_WITH_MASK(FALSEVALUE, TRUEVALUE, MASK)       vbslq_u64(MASK, TRUEVALUE, FALSEVALUE)


//Trig
#define ARM64_ARC_COSINE(X, order)                  inv_cosine_arm64(X, order)

//Max
#define ARM64_MAX_FLOATS(X,Y)               vmaxq_f64(X,Y)

//Absolute value
#define ARM64_ABS_FLOAT(X)                  vabsq_f64(X)
#endif //DOUBLE_PREC

#include "fast_acos.h"
    
static inline ARM64_FLOATS inv_cosine_arm64_DOUBLE(const ARM64_FLOATS X, const int order)
{
    union cos{
        ARM64_FLOATS m;
        DOUBLE x[ARM64_NVEC];
    };
    union cos union_costheta;
    union cos union_returnvalue;
    union_costheta.m = X;
    const DOUBLE minus_one = (DOUBLE) -1.0;
    const DOUBLE one = (DOUBLE) 1.0;

    //Force everything to be in range [0,1]
    for(int ii=0;ii<ARM64_NVEC;ii++) {
        const DOUBLE costheta = union_costheta.x[ii];
        union_costheta.x[ii] = costheta <= minus_one ? minus_one:costheta;
        union_costheta.x[ii] = costheta >= one ? one:costheta;
    }
    
    if(order == 0) {
        for(int ii=0;ii<ARM64_NVEC;ii++) {
            const DOUBLE costheta = union_costheta.x[ii];
            union_returnvalue.x[ii] = ACOS(costheta);
        }
    } else {
        //fast acos
        /*Taken from associated C++ code in http://www.geometrictools.com/GTEngine/Include/Mathematics/GteACosEstimate.h*/
        for(int ii=0;ii<ARM64_NVEC;ii++) {
            union_returnvalue.x[ii] = FAST_ACOS(union_costheta.x[ii]);
        }
    }
    return union_returnvalue.m;
  }

  
//The three different unions used
//for computing rpavg and weightavg
union int4 {
  ARM64_INTS m_ibin;
  ARM64_INT_TYPE ibin[ARM64_NVEC];
};

union float4{
  ARM64_FLOATS m_Dperp;
  DOUBLE Dperp[ARM64_NVEC];
};

union float4_weights{
  ARM64_FLOATS m_weights;
  DOUBLE weights[ARM64_NVEC];
};


#ifdef __cplusplus
}
#endif
