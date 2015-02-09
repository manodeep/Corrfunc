### Set the compiler -- options are icc/gcc/clang. 
CC=gcc
#### Add any compiler specific flags you want
CFLAGS=
#### Add any compiler specific link flags you want
CLINK=

### You should NOT edit below this line
DISTNAME=corrfunc
MINOR=0
MAJOR=1

INCLUDE=-I../../io -I../../utils 

### The POSIX_SOURCE flag is required to get the definition of strtok_r
CFLAGS += -Wsign-compare -Wall -Wextra -Wshadow -Wunused -std=c99 -g -m64 -fPIC -D_POSIX_SOURCE -D_DARWIN_C_SOURCE -O3 -Ofast
GSL_CFLAGS := $(shell gsl-config --cflags) 
GSL_LIBDIR := $(shell gsl-config --prefix)/lib
GSL_LINK   := $(shell gsl-config --libs) -Xlinker -rpath -Xlinker $(GSL_LIBDIR) 

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
  CFLAGS += -xhost -opt-prefetch  #-vec-report6  
  ifeq (USE_OMP,$(findstring USE_OMP,$(OPT)))
		CFLAGS += -openmp
		CLINK  += -openmp 
  endif
else

  ### compiler specific flags for gcc
  ifeq (gcc,$(findstring gcc,$(CC)))
		CFLAGS += -ftree-vectorize -funroll-loops -fprefetch-loop-arrays #-ftree-vectorizer-verbose=6 -fopt-info-vec-missed #-fprofile-use -fprofile-correction #-fprofile-generate
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
  # CFLAGS  += -march=native
	CFLAGS  += -Wformat=2  -Wpacked  -Wnested-externs -Wpointer-arith  -Wredundant-decls  -Wfloat-equal -Wcast-qual  
  CFLAGS  +=  -Wcast-align -Wmissing-declarations -Wmissing-prototypes  -Wnested-externs -Wstrict-prototypes  #-D_POSIX_C_SOURCE=2 -Wpadded -Wconversion
  CLINK += -lm
endif


ifeq (USE_MKL,$(findstring USE_MKL,$(OPT)))
	BLAS_INCLUDE:=-DMKL_ILP64 -m64 -I$(MKLROOT)/include 
  ##Use the Intel MKL library. Check the compiler + openmp
	ifneq (USE_OMP,$(findstring USE_OMP,$(OPT)))
    ##Link+include sequential libraries
		ifeq (icc,$(findstring icc,$(CC)))
      ##icc with Intel MKL
			BLAS_LINK:= -L$(MKLROOT)/lib/intel64 -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential -lpthread -lm
		else
	    ##gcc with Intel MKL
			BLAS_LINK:= -Wl,--no-as-needed -L$(MKLROOT)/lib -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential -lpthread -lm
		endif
	else
		ifeq (icc,$(findstring icc,$(CC)))
      ##icc with Intel MKL+OpenMP
			BLAS_LINK:= -L$(MKLROOT)/lib -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -lpthread -lm
		else
	    ##gcc with Intel MKL
			BLAS_LINK:= -Wl,--no-as-needed -L$(MKLROOT)/lib -lmkl_intel_ilp64 -lmkl_core -lmkl_gnu_thread -ldl -lpthread -lm
		endif
	endif
else
##Use some OpenMP parallel BLAS library (OpenBlas/ATLAS, for instance)
BLAS_INCLUDE:=
BLAS_LINK:=
endif

ifeq (OUTPUT_THETAAVG,$(findstring OUTPUT_THETAAVG,$(OPT)))
  ifneq (DOUBLE_PREC,$(findstring DOUBLE_PREC,$(OPT)))
    $(error DOUBLE_PREC must be enabled with OUTPUT_THETAAVG -- loss of precision will give you incorrect results for the outer bins (>=20-30 million pairs))
  endif
  ifeq (USE_AVX,$(findstring USE_AVX,$(OPT)))
     ifneq (icc,$(findstring icc,$(CC)))
        $(error OUTPUT_THETAAVG with AVX capabilties will only work with icc)
     endif
  endif
endif

ifeq (FAST_DIVIDE,$(findstring FAST_DIVIDE,$(OPT)))
  ifneq (USE_AVX,$(findstring USE_AVX,$(OPT)))
    $(warning Makefile option FAST_DIVIDE will not do anything unless USE_AVX is set)
  endif
endif
