/* File: cpu_features.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/

  Adapted from Agner Fog's vectorclass: http://agner.org/
*/

#include <stdio.h>

#include "cpu_features.h"

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
    iset = FALLBACK;                                       // default value

#ifdef __ARM_NEON
    return iset; /* return FALLBACK as the ISA on ARM64 */
#endif

    int abcd[4] = {0,0,0,0};                               // cpuid results
    cpuid(abcd, 0);                                        // call cpuid function 0
    if (abcd[0] == 0) return iset;                         // no further cpuid function supported
    cpuid(abcd, 1);                                        // call cpuid function 1 for feature flags
    if ((abcd[3] & (1 <<  0)) == 0) return iset;           // no floating point
    if ((abcd[3] & (1 << 23)) == 0) return iset;           // no MMX
    if ((abcd[3] & (1 << 15)) == 0) return iset;           // no conditional move
    if ((abcd[3] & (1 << 24)) == 0) return iset;           // no FXSAVE
    if ((abcd[3] & (1 << 25)) == 0) return iset;           // no SSE
    iset = SSE;                                            // 1: SSE supported

    if ((abcd[3] & (1 << 26)) == 0) return iset;           // no SSE2
    iset = SSE2;                                           // 2: SSE2 supported

    if ((abcd[2] & (1 <<  0)) == 0) return iset;           // no SSE3
    iset = SSE3;                                           // 3: SSE3 supported

    if ((abcd[2] & (1 <<  9)) == 0) return iset;           // no SSSE3
    iset = SSSE3;                                          // 4: SSSE3 supported

    if ((abcd[2] & (1 << 19)) == 0) return iset;           // no SSE4.1
    iset = SSE4;                                           // 5: SSE4.1 supported

    if ((abcd[2] & (1 << 23)) == 0) return iset;           // no POPCNT
    if ((abcd[2] & (1 << 20)) == 0) return iset;           // no SSE4.2
    iset = SSE42;                                          // 6: SSE4.2 supported

    if ((abcd[2] & (1 << 27)) == 0) return iset;           // no OSXSAVE
    if ((xgetbv(0) & 6) != 6)       return iset;           // AVX not enabled in O.S.
    if ((abcd[2] & (1 << 28)) == 0) return iset;           // no AVX
    iset = AVX;                                            // 7: AVX supported

    cpuid(abcd, 7);                                        // call cpuid leaf 7 for feature flags
    if ((abcd[1] & (1 <<  5)) == 0) return iset;           // no AVX2
    iset = AVX2;                                           // 8: AVX2 supported

    cpuid(abcd, 0xD);                                      // call cpuid leaf 0xD for feature flags
    if ((abcd[0] & 0x60) != 0x60)   return iset;           // no AVX512
    iset = AVX512F;                                        // 9: AVX512F supported
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
        case AVX512F:
#ifdef __AVX512F__
            iset = AVX512F;
            break;
#elif defined(GAS_BUG_DISABLE_AVX512)
            fprintf(stderr, "[Warning] AVX512F is disabled due to a GNU Assembler bug.  Upgrade to binutils >= 2.32 to fix this.\n");
#else
            fprintf(stderr, "[Warning] The CPU supports AVX512F but the compiler does not.  Can you try another compiler?\n");
#endif
            // fall through
        case AVX2:
#ifdef __AVX2__
            iset = AVX2;
            break;
#else
            fprintf(stderr, "[Warning] The CPU supports AVX2 but the compiler does not.  Can you try another compiler?\n");
#endif
            // fall through
        case AVX:
#ifdef __AVX__
            iset = AVX;
            break;
#else
            fprintf(stderr, "[Warning] The CPU supports AVX but the compiler does not.  Can you try another compiler?\n");
#endif
            // fall through
        case SSE42:
#ifdef __SSE4_2__
            iset = SSE42;
            break;
#else
            fprintf(stderr, "[Warning] The CPU supports SSE4.2 but the compiler does not.  Can you try another compiler?\n");
#endif
            // fall through
        case SSE4:
#ifdef __SSE4_1__
            iset = SSE4;
            break;
#else
            fprintf(stderr, "[Warning] The CPU supports SSE4.1 but the compiler does not.  Can you try another compiler?\n");
#endif
            // fall through
        case SSSE3:
#ifdef __SSSE3__
            iset = SSSE3;
            break;
#else
            fprintf(stderr, "[Warning] The CPU supports SSSE3 but the compiler does not.  Can you try another compiler?\n");
#endif
            // fall through
        case SSE3:
#ifdef __SSE3__
            iset = SSE3;
            break;
#else
            fprintf(stderr, "[Warning] The CPU supports SSE3 but the compiler does not.  Can you try another compiler?\n");
#endif
            // fall through
        case SSE2:
#ifdef __SSE2__
            iset = SSE2;
            break;
#else
            fprintf(stderr, "[Warning] The CPU supports SSE2 but the compiler does not.  Can you try another compiler?\n");
#endif
            // fall through
        case SSE:
#ifdef __SSE__
            iset = SSE;
            break;
#else
            fprintf(stderr, "[Warning] The CPU supports SSE but the compiler does not.  Can you try another compiler?\n");
#endif
            // fall through
// This will need to be uncommented when the ARM64 kernels are added in
//         case ARM64:
// #ifdef __ARM_NEON
//             iset = ARM64;
//             break;
// #else
//             fprintf(stderr, "[Warning] The CPU supports NEON but the compiler does not.  Can you try another compiler?\n");
// #endif
        // fall through
        case FALLBACK:
        default:
            iset = FALLBACK;
            break;
    }

    return iset;
}
