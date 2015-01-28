/* File: ftread.h */
/*
		This file is a part of the Corrfunc package
		Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
		License: MIT LICENSE. See LICENSE file under the top-level
		directory at https://bitbucket.org/manodeep/corrfunc/
*/

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
