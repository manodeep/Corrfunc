/* File: cpu_features.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/

  Adapted from Agner Fog's vectorclass: http://agner.org/
*/

#include "cpu_features.h"
#include <stdio.h>

// Use CPUID to detect what instruction sets the CPU supports
// The compiler may not support all these features though!
// Use get_max_usable_isa() to find the max ISA supported
// by both the compiler and CPU
int runtime_instrset_detect(void)
{
    static int iset = -1;                                  // remember value for next call
    if (iset >= 0) {
        return iset;                                       // called before
    }
    iset = 0;                                              // default value
    int abcd[4] = {0,0,0,0};                               // cpuid results
    cpuid(abcd, 0);                                        // call cpuid function 0
    if (abcd[0] == 0) return iset;                         // no further cpuid function supported
    cpuid(abcd, 1);                                        // call cpuid function 1 for feature flags
    if ((abcd[3] & (1 <<  0)) == 0) return iset;           // no floating point
    if ((abcd[3] & (1 << 23)) == 0) return iset;           // no MMX
    if ((abcd[3] & (1 << 15)) == 0) return iset;           // no conditional move
    if ((abcd[3] & (1 << 24)) == 0) return iset;           // no FXSAVE
    if ((abcd[3] & (1 << 25)) == 0) return iset;           // no SSE
    iset = 1;                                              // 1: SSE supported

    if ((abcd[3] & (1 << 26)) == 0) return iset;           // no SSE2
    iset = 2;                                              // 2: SSE2 supported

    if ((abcd[2] & (1 <<  0)) == 0) return iset;           // no SSE3
    iset = 3;                                              // 3: SSE3 supported

    if ((abcd[2] & (1 <<  9)) == 0) return iset;           // no SSSE3
    iset = 4;                                              // 4: SSSE3 supported

    if ((abcd[2] & (1 << 19)) == 0) return iset;           // no SSE4.1
    iset = 5;                                              // 5: SSE4.1 supported

    if ((abcd[2] & (1 << 23)) == 0) return iset;           // no POPCNT
    if ((abcd[2] & (1 << 20)) == 0) return iset;           // no SSE4.2
    iset = 6;                                              // 6: SSE4.2 supported

    if ((abcd[2] & (1 << 27)) == 0) return iset;           // no OSXSAVE
    if ((xgetbv(0) & 6) != 6)       return iset;           // AVX not enabled in O.S.
    if ((abcd[2] & (1 << 28)) == 0) return iset;           // no AVX
    iset = 7;                                              // 7: AVX supported

    cpuid(abcd, 7);                                        // call cpuid leaf 7 for feature flags
    if ((abcd[1] & (1 <<  5)) == 0) return iset;           // no AVX2
    iset = 8;                                              // 8: AVX2 supported

    cpuid(abcd, 0xD);                                      // call cpuid leaf 0xD for feature flags
    if ((abcd[0] & 0x60) != 0x60)   return iset;           // no AVX512
    iset = 9;                                              // 9: AVX512F supported
    return iset;
}

// Report the max ISA supported by both the CPU and compiler
int get_max_usable_isa(void)
{
    static int iset = -1;                                  // remember value for next call
    if (iset >= 0) {
        return iset;                                       // called before
    }
    iset = runtime_instrset_detect();

    switch(iset){
        case 9:
#ifdef __AVX512F__
            iset = 9;
            break;
#elif defined(GAS_BUG_DISABLE_AVX512)
            fprintf(stderr, "[Warning] AVX512 is disabled due to a GNU Assembler bug.  Upgrade to binutils >= 2.32 to fix this.\n");
#else
            fprintf(stderr, "[Warning] The CPU supports AVX512 but the compiler does not.  Can you try another compiler?\n");
#endif
        case 8:
#ifdef __AVX2__
            iset = 8;
            break;
#else
            fprintf(stderr, "[Warning] The CPU supports AVX2 but the compiler does not.  Can you try another compiler?\n");
#endif
        case 7:
#ifdef __AVX__
            iset = 7;
            break;
#else
            fprintf(stderr, "[Warning] The CPU supports AVX but the compiler does not.  Can you try another compiler?\n");
#endif
        case 6:
#ifdef __SSE4_2__
            iset = 6;
            break;
#else
            fprintf(stderr, "[Warning] The CPU supports SSE4.2 but the compiler does not.  Can you try another compiler?\n");
#endif
        case 5:
#ifdef __SSE4_1__
            iset = 5;
            break;
#else
            fprintf(stderr, "[Warning] The CPU supports SSE4.1 but the compiler does not.  Can you try another compiler?\n");
#endif
        case 4:
#ifdef __SSSE3__
            iset = 4;
            break;
#else
            fprintf(stderr, "[Warning] The CPU supports SSSE3 but the compiler does not.  Can you try another compiler?\n");
#endif
        case 3:
#ifdef __SSE3__
            iset = 3;
            break;
#else
            fprintf(stderr, "[Warning] The CPU supports SSE3 but the compiler does not.  Can you try another compiler?\n");
#endif
        case 2:
#ifdef __SSE2__
            iset = 2;
            break;
#else
            fprintf(stderr, "[Warning] The CPU supports SSE2 but the compiler does not.  Can you try another compiler?\n");
#endif
        case 1:
#ifdef __SSE__
            iset = 1;
            break;
#else
            fprintf(stderr, "[Warning] The CPU supports SSE but the compiler does not.  Can you try another compiler?\n");
#endif
        case 0:
        default:
            iset = 0;
            break;
    }

    return iset;
}
