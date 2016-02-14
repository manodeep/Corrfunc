/* File: progressbar.h */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif


    void init_my_progressbar(const int64_t N,int *interrupted);
    void my_progressbar(const int64_t curr_index,int *interrupted);
    void finish_myprogressbar(int *interrupted);

#ifdef __cplusplus
}
#endif
