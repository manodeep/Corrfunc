/* File: cpu_features.h */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/


  Adapted from Agner Fog's vectorclass: http://agner.org/
*/

#pragma once
#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
 extern "C" {
#endif

typedef enum {
  DEFAULT=-42,/* present simply to make the enum a signed int*/
  FALLBACK=0, /* No special options */
  SSE=1,  /* 64 bit vectors */
  SSE2=2, /* 128 bit vectors */
  SSE3=3, /* 128 bit vectors */
  SSSE3=4, /* 128 bit vectors */
  SSE4=5,/* 128bit vectors */
  SSE42=6, /* 128bit vectors with blend operations */
  AVX=7, /* 256bit vector width */
  AVX2=8,  /* AVX2 (integer operations)*/
  AVX512F=9,/* AVX 512 Foundation */
  ARM64=10, /* ARM64 on Apple M1/M2 */
  NUM_ISA  /*NUM_ISA will be the next integer after
            the last declared enum. AVX512F:=9 (so, NUM_ISA==10)*/
} isa;  //name for instruction sets -> corresponds to the return values for functions in cpu_features.c


static inline void cpuid (int output[4], int functionnumber) {	
#if defined (__ARM_NEON)
    (void) output;
    (void) functionnumber;
#else
#if defined(__GNUC__) || defined(__clang__)              // use inline assembly, Gnu/AT&T syntax

   int a, b, c, d;
   __asm("cpuid" : "=a"(a),"=b"(b),"=c"(c),"=d"(d) : "a"(functionnumber),"c"(0) );
   output[0] = a;
   output[1] = b;
   output[2] = c;
   output[3] = d;

#else                                                      // unknown platform. try inline assembly with masm/intel syntax

    __asm {
        mov eax, functionnumber
        xor ecx, ecx
        cpuid;
        mov esi, output
        mov [esi],    eax
        mov [esi+4],  ebx
        mov [esi+8],  ecx
        mov [esi+12], edx
    }

#endif
#endif
}

// Define interface to xgetbv instruction
static inline int64_t xgetbv (int ctr) {	
#if defined(__ARM_NEON)
    (void) ctr;
    return (int64_t) FALLBACK; /* use FALLBACK kernels until the ARM64 kernels are added to the pair-counters */
#else
#if (defined (__INTEL_COMPILER) && __INTEL_COMPILER >= 1200) //Intel compiler supporting _xgetbv intrinsic
    return _xgetbv(ctr);                                   // intrinsic function for XGETBV
#elif defined(__GNUC__)                                    // use inline assembly, Gnu/AT&T syntax
   uint32_t a, d;
   __asm("xgetbv" : "=a"(a),"=d"(d) : "c"(ctr) : );
   return a | (((uint64_t) d) << 32);
#else  
   uint32_t a, d;
    __asm {
        mov ecx, ctr
        _emit 0x0f
        _emit 0x01
        _emit 0xd0 ; // xgetbv
        mov a, eax
        mov d, edx
    }
    return a | (((uint64_t) d) << 32);

#endif
#endif
}

extern int runtime_instrset_detect(void);
extern int get_max_usable_isa(void);

#ifdef __cplusplus
}
#endif
