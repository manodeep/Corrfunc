/* File: fast_acos.h */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#pragma once

#include <stdio.h>
#include "function_precision.h"

static inline DOUBLE FAST_ACOS(const DOUBLE x)
{    
    /* Taken from http://www.geometrictools.com/GTEngine/Include/Mathematics/GteConstants.h */    

    /*Taken from associated C++ code in http://www.geometrictools.com/GTEngine/Include/Mathematics/GteACosEstimate.h*/
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
    DOUBLE one = (DOUBLE) 1.0;
    poly = (DOUBLE)GTE_C_ACOS_DEG8_C8;
    poly = (DOUBLE)GTE_C_ACOS_DEG8_C7 + poly * x;
    poly = (DOUBLE)GTE_C_ACOS_DEG8_C6 + poly * x;
    poly = (DOUBLE)GTE_C_ACOS_DEG8_C5 + poly * x;
    poly = (DOUBLE)GTE_C_ACOS_DEG8_C4 + poly * x;
    poly = (DOUBLE)GTE_C_ACOS_DEG8_C3 + poly * x;
    poly = (DOUBLE)GTE_C_ACOS_DEG8_C2 + poly * x;
    poly = (DOUBLE)GTE_C_ACOS_DEG8_C1 + poly * x;
    poly = (DOUBLE)GTE_C_ACOS_DEG8_C0 + poly * x;
    poly = poly * SQRT(one - x);

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
    return poly;
}
