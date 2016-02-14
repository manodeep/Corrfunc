/* File: defs.h */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#pragma once


#define ADD_DIFF_TIME(t0,t1)     ((t1.tv_sec - t0.tv_sec) + 1e-6*(t1.tv_usec - t0.tv_usec))
#define ALIGNMENT                32
