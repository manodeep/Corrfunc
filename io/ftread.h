/* File: ftread.h */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#pragma once

#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

    int ftread(void *ptr,size_t size,size_t nitems,FILE *stream) __attribute__((warn_unused_result));
    int my_ftread(void *ptr,size_t size,size_t nitems,FILE *stream) __attribute__((warn_unused_result));

#ifdef __cplusplus
}
#endif
