#### Science use-cases for Theory Correlation Functions
OPT = -DPERIODIC
#OPT += -DOUTPUT_RPAVG  ### Enabling this can cause up to a 2x performance hit


#### Extra options for Data Correlation Functions
DATA_OPT += -DLINK_IN_DEC
DATA_OPT += -DLINK_IN_RA

#### Code specs for both theory and data Correlation Functions
#OPT += -DDOUBLE_PREC
OPT += -DUSE_AVX
OPT += -DUSE_OMP

GSL_CFLAGS := $(shell gsl-config --cflags) 
GSL_LINK   := $(shell gsl-config --ldflags)
GSL_LIBDIR := $(shell gsl-config --prefix)/lib


### Set the compiler -- options are icc/gcc/clang
CC=gcc
#### Add any compiler specific flags you want
CFLAGS= 
#### Add any compiler specific link flags you want
CLINK=

#### Add any compiler specific flags you want
DATA_CFLAGS=
DATA_CLINK=

### You should NOT edit below this line
DISTNAME=corrfunc
MINOR=0
MAJOR=1

INCLUDE=-I../../io -I../../utils 
CFLAGS += -Wsign-compare -Wall -Wextra -Wshadow -Wunused -std=c99 -g -m64 -fPIC -O3  #-Werror


ifneq (USE_OMP,$(findstring USE_OMP,$(OPT)))
  ifneq (clang,$(findstring clang,$(CC)))
     $(warning Recommended compiler for a serial build is clang)
  endif
endif

ifeq (OUTPUT_RPAVG,$(findstring OUTPUT_RPAVG,$(OPT)))
  ifneq (DOUBLE_PREC,$(findstring DOUBLE_PREC,$(OPT)))
    $(error DOUBLE_PREC must be enabled with OUTPUT_RPAVG -- loss of precision will give you incorrect results for the outer bins (>=20-30 million pairs))
  endif
endif

ifneq (DOUBLE_PREC,$(findstring DOUBLE_PREC,$(OPT)))
	VECTOR_TYPE:=float
else
	VECTOR_TYPE:=double
endif


ifeq (icc,$(findstring icc,$(CC)))
  CFLAGS += -xhost -opt-prefetch -ipo  #-vec-report6  
  ifeq (USE_OMP,$(findstring USE_OMP,$(OPT)))
		CFLAGS += -openmp
		CLINK  += -openmp 
  endif
else

  ### compiler specific flags for gcc
  ifeq (gcc,$(findstring gcc,$(CC)))
		CFLAGS += -ftree-vectorize -funroll-loops -fprefetch-loop-arrays #-fprofile-use -fprofile-correction #-fprofile-generate
    ifeq (USE_OMP,$(findstring USE_OMP,$(OPT)))
			CFLAGS += -fopenmp
			CLINK  += -fopenmp
    endif
  endif

  ### compiler specific flags for clang
  ifeq (clang,$(findstring clang,$(CC)))
		CFLAGS += -funroll-loops
    ifeq (USE_OMP,$(findstring USE_OMP,$(OPT)))
      $(error clang does not support OpenMP - please use gcc/icc for compiling with openmp)
     endif
  endif

  ifeq (USE_AVX,$(findstring USE_AVX,$(OPT)))
    CFLAGS  +=  -mavx -mpopcnt
  endif

  #### common options for gcc and clang
  CFLAGS  += -march=native -Wformat=2  -Wpacked  -Wnested-externs -Wpointer-arith  -Wredundant-decls  -Wfloat-equal -Wcast-qual  
  CFLAGS  +=  -Wcast-align -Wmissing-declarations -Wmissing-prototypes  -Wnested-externs -Wstrict-prototypes  #-D_POSIX_C_SOURCE=2 -Wpadded -Wconversion
  CLINK += -lm
endif


