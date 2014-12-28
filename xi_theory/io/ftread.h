#ifndef __FTREAD_H
#define __FTREAD_H
#include <stdio.h>
#include <stdlib.h>


int ftread(void *ptr,size_t size,size_t nitems,FILE *stream);
int my_ftread(void *ptr,size_t size,size_t nitems,FILE *stream);

#endif
