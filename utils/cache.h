/* File: cache.h */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#pragma once


#ifdef __cplusplus
extern "C" {
#endif

#define DEFAULT_CACHE_LINE_SIZE  64

/*Taken from http://stackoverflow.com/questions/794632/programmatically-get-the-cache-line-size*/
#include <stddef.h>
size_t cache_line_size(void);

#ifdef __cplusplus
}
#endif
