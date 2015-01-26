#ifndef __FTREAD_H
#define __FTREAD_H

#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

int ftread(void *ptr,size_t size,size_t nitems,FILE *stream);
int my_ftread(void *ptr,size_t size,size_t nitems,FILE *stream);

#ifdef __cplusplus
}
#endif
	
#endif
