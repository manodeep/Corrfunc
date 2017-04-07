/* File: fast_acos.h */
/*
  This file is covered by the original Boost License
  for http://www.geometrictools.com/.

  Distributed as part of Corrfunc package
  at https://github.com/manodeep/Corrfunc/
*/

#pragma once

#include <stdio.h>
#include <math.h>
#include "function_precision.h"


#if 0
//~20% faster than math acos at the expense of accuracy. 

// Taken from http://stackoverflow.com/questions/11261170/c-and-maths-fast-approximation-of-a-trigonometric-function
/* minimax approximation to arcsin on [0, 0.5625] with rel. err. ~= 1.5e-11 */
static inline DOUBLE asin_core (DOUBLE x)
{
    const DOUBLE x2 = x * x;
    const DOUBLE x4 = x2 * x2;
    const DOUBLE x8 = x4 * x4;
    /* evaluate polynomial using a mix of Estrin's and Horner's scheme */
    return (((4.5334220547132049e-2 * x2 - 1.1226216762576600e-2) * x4 +
             (2.6334281471361822e-2 * x2 + 2.0596336163223834e-2)) * x8 +
            (3.0582043602875735e-2 * x2 + 4.4630538556294605e-2) * x4 +
            (7.5000364034134126e-2 * x2 + 1.6666666300567365e-1)) * x2 * x + x;
}

static inline DOUBLE FAST_ACOS(const DOUBLE x)
{
    const DOUBLE xa = FABS(x);
    const DOUBLE one = 1.0;
    const DOUBLE half = 0.5;
    const DOUBLE two = 2.0;
    /* arcsin(x) = pi/2 - 2 * arcsin (sqrt ((1-x) / 2)) 
     * arccos(x) = pi/2 - arcsin(x)
     * arccos(x) = 2 * arcsin (sqrt ((1-x) / 2))
     */
    DOUBLE t;
    if (xa > 0.5625) {
        t = two * asin_core (SQRT (half * (one - xa)));
    } else {
        t = M_PI_2 - asin_core (xa);
    }
    /* arccos (-x) = pi - arccos(x) */
    return (x < ZERO) ? (M_PI - t) : t;
}
   

#else 

static inline DOUBLE FAST_ACOS(const DOUBLE x)
{    
    /* Taken from http://www.geometrictools.com/GTEngine/Include/Mathematics/GteConstants.h */    
    /*Taken from associated C++ code in http://www.geometrictools.com/GTEngine/Include/Mathematics/GteACosEstimate.h*/
    /* 8th degree polynomial, valid on [0, 1], max absolute error 3.6e-9 */

    /* Faster but with a larger error bound than the alternate FAST_ACOS implementation above */
    
    DOUBLE poly;
#define GTE_C_ACOS_DEG8_C0 +1.5707963267948966
#define GTE_C_ACOS_DEG8_C1 -2.1460143648688035e-01
#define GTE_C_ACOS_DEG8_C2 +8.9034700107934128e-02
#define GTE_C_ACOS_DEG8_C3 -5.0625279962389413e-02
#define GTE_C_ACOS_DEG8_C4 +3.2683762943179318e-02
#define GTE_C_ACOS_DEG8_C5 -2.0949278766238422e-02
#define GTE_C_ACOS_DEG8_C6 +1.1272900916992512e-02
#define GTE_C_ACOS_DEG8_C7 -4.1160981058965262e-03
#define GTE_C_ACOS_DEG8_C8 +7.1796493341480527e-04
#define GTE_C_ACOS_DEG8_MAX_ERROR 3.6340015129032732e-9
    const DOUBLE xa = FABS(x);
    DOUBLE one = (DOUBLE) 1.0;
    poly = (DOUBLE)GTE_C_ACOS_DEG8_C8;
    poly = (DOUBLE)GTE_C_ACOS_DEG8_C7 + poly * xa;
    poly = (DOUBLE)GTE_C_ACOS_DEG8_C6 + poly * xa;
    poly = (DOUBLE)GTE_C_ACOS_DEG8_C5 + poly * xa;
    poly = (DOUBLE)GTE_C_ACOS_DEG8_C4 + poly * xa;
    poly = (DOUBLE)GTE_C_ACOS_DEG8_C3 + poly * xa;
    poly = (DOUBLE)GTE_C_ACOS_DEG8_C2 + poly * xa;
    poly = (DOUBLE)GTE_C_ACOS_DEG8_C1 + poly * xa;
    poly = (DOUBLE)GTE_C_ACOS_DEG8_C0 + poly * xa;
    poly = poly * SQRT(one - xa);

#undef GTE_C_ACOS_DEG8_C0 
#undef GTE_C_ACOS_DEG8_C1 
#undef GTE_C_ACOS_DEG8_C2 
#undef GTE_C_ACOS_DEG8_C3 
#undef GTE_C_ACOS_DEG8_C4 
#undef GTE_C_ACOS_DEG8_C5 
#undef GTE_C_ACOS_DEG8_C6 
#undef GTE_C_ACOS_DEG8_C7 
#undef GTE_C_ACOS_DEG8_C8 
#undef GTE_C_ACOS_DEG8_MAX_ERROR 

    return (x < ZERO) ? (M_PI - poly) : poly;
}

#endif
