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


static inline void cpuid (int output[4], int functionnumber) {	
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
}

// Define interface to xgetbv instruction
static inline int64_t xgetbv (int ctr) {	
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
}

extern int instrset_detect(void);

#ifdef __cplusplus
}
#endif
