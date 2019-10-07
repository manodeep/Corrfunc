/* File: cpu_features.c */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/

  Adapted from Agner Fog's vectorclass: http://agner.org/
*/

#include "cpu_features.h"

int instrset_detect(void)
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
#ifdef __SSE__
    iset = 1;                                              // 1: SSE supported
#endif

    if ((abcd[3] & (1 << 26)) == 0) return iset;           // no SSE2
#ifdef __SSE2__
    iset = 2;                                              // 2: SSE2 supported
#endif

    if ((abcd[2] & (1 <<  0)) == 0) return iset;           // no SSE3
#ifdef __SSE3__
    iset = 3;                                              // 3: SSE3 supported
#endif

    if ((abcd[2] & (1 <<  9)) == 0) return iset;           // no SSSE3
#ifdef __SSSE3__
    iset = 4;                                              // 4: SSSE3 supported
#endif

    if ((abcd[2] & (1 << 19)) == 0) return iset;           // no SSE4.1
#ifdef __SSE4_1__
    iset = 5;                                              // 5: SSE4.1 supported
#endif

    if ((abcd[2] & (1 << 23)) == 0) return iset;           // no POPCNT
    if ((abcd[2] & (1 << 20)) == 0) return iset;           // no SSE4.2
#ifdef __SSE4_2__
    iset = 6;                                              // 6: SSE4.2 supported
#endif

    if ((abcd[2] & (1 << 27)) == 0) return iset;           // no OSXSAVE
    if ((xgetbv(0) & 6) != 6)       return iset;           // AVX not enabled in O.S.
    if ((abcd[2] & (1 << 28)) == 0) return iset;           // no AVX
#ifdef __AVX__
    iset = 7;                                              // 7: AVX supported
#endif

    cpuid(abcd, 7);                                        // call cpuid leaf 7 for feature flags
    if ((abcd[1] & (1 <<  5)) == 0) return iset;           // no AVX2
#ifdef __AVX2__
    iset = 8;                                              // 8: AVX2 supported
#endif

    cpuid(abcd, 0xD);                                      // call cpuid leaf 0xD for feature flags
    if ((abcd[0] & 0x60) != 0x60)   return iset;           // no AVX512
#ifdef __AVX512F__
    iset = 9;                                              // 9: AVX512F supported
#endif
    return iset;
}

