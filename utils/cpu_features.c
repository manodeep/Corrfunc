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

