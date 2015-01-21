#pragma once


#define ADD_DIFF_TIME(t0,t1)     ((t1.tv_sec - t0.tv_sec) + 1e-6*(t1.tv_usec - t0.tv_usec))


