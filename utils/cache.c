/* File: cache.h */
/*
  This file is a part of the Corrfunc package
  Copyright (C) 2015-- Manodeep Sinha (manodeep@gmail.com)
  License: MIT LICENSE. See LICENSE file under the top-level
  directory at https://github.com/manodeep/Corrfunc/
*/

/*Taken from http://stackoverflow.com/questions/794632/programmatically-get-the-cache-line-size*/

// Author: Nick Strupat
// Date: October 29, 2010
// Returns the cache line size (in bytes) of the processor, or 0 on failure


/* BUG NOTE: Might be worthwhile to take the minimum cache linesize across all cpus
   just in case the cpus have different caches. This has already happened on ARM
   https://github.com/mono/mono/pull/3549
   http://lists.infradead.org/pipermail/linux-arm-kernel/2016-September/453859.html
*/


#if defined(__APPLE__)

#include <sys/sysctl.h>
size_t cache_line_size(void)
{
    size_t line_size = 0;
    size_t sizeof_line_size = sizeof(line_size);
    sysctlbyname("hw.cachelinesize", &line_size, &sizeof_line_size, 0, 0);
    return line_size;
}

#elif defined(__linux__)

#include <stdio.h>
size_t cache_line_size(void)
{
    FILE *fp = fopen("/sys/devices/system/cpu/cpu0/cache/index0/coherency_line_size", "r");
    size_t lineSize = DEFAULT_CACHE_LINE_SIZE;
    if (fp != NULL) {
        int nitems = fscanf(fp, "%zu", &lineSize);
        if(nitems !=1)  {
            linesize=0;
        }
        fclose(fp);
    }
    return lineSize;
}

#else

#warning Unrecognized platform for figuring out cache line size - returning a default value of DEFAULT_CACHE_LINE_SIZE
size_t cache_line_size(void)
{
    return DEFAULT_CACHE_LINE_SIZE;//set a default of 64 bytes for cache line size
}

#endif
