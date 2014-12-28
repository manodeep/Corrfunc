#pragma once
#include <stdio.h>
#include <stdlib.h>

unsigned int read_positions(const char *filename, const char *format, void **xpos, void **ypos, void **zpos, const size_t size);
