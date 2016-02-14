/* File: io.h */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <stdarg.h>
#include <string.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#include "function_precision.h"
#include "ftread.h"
#include "utils.h"

    int64_t read_positions(const char *filename, const char *format, const size_t size, const int num_fields, ...) __attribute__((warn_unused_result));

#ifdef __cplusplus
}
#endif
