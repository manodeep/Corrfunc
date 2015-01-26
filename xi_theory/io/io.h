#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>


#ifdef __cplusplus
extern "C" {
#endif

int64_t read_positions(const char *filename, const char *format, void **xpos, void **ypos, void **zpos, const size_t size);

#ifdef __cplusplus
}
#endif
