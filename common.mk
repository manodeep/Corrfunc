### Set the default compiler -- options are icc/gcc/clang. 
CC:=gcc

#### Add any compiler specific flags you want
CFLAGS:=

#### Add any compiler specific link flags you want
CLINK:=

### You should NOT edit below this line
DISTNAME:=Corrfunc
MAJOR:=0
MINOR:=0
PATCHLEVEL:=1
VERSION:=$(MAJOR).$(MINOR).$(PATCHLEVEL)

INCLUDE:=-I../../io -I../../utils
### The POSIX_SOURCE flag is required to get the definition of strtok_r
CFLAGS += -Wsign-compare -Wall -Wextra -Wshadow -Wunused -std=c99 -g -m64 -fPIC -D_POSIX_SOURCE -D_DARWIN_C_SOURCE -O3 #-Ofast
GSL_CFLAGS := $(shell gsl-config --cflags) 
GSL_LIBDIR := $(shell gsl-config --prefix)/lib
GSL_LINK   := $(shell gsl-config --libs) -Xlinker -rpath -Xlinker $(GSL_LIBDIR) 


COMPILE_PYTHON_EXT := 1
PYTHON_VERSION_FULL := $(wordlist 2,4,$(subst ., ,$(shell python --version 2>&1)))
PYTHON_VERSION_MAJOR := $(word 1,${PYTHON_VERSION_FULL})
PYTHON_VERSION_MINOR := $(word 2,${PYTHON_VERSION_FULL})

## I only need this so that I can print out the full python version (correctly)
## in case of error
PYTHON_VERSION_PATCH := $(word 3,${PYTHON_VERSION_FULL})

## Check numpy version
NUMPY_VERSION_FULL :=  $(wordlist 1,3,$(subst ., ,$(shell python -c "from __future__ import print_function; import numpy; print(numpy.__version__)")))
NUMPY_VERSION_MAJOR := $(word 1,${NUMPY_VERSION_FULL})
NUMPY_VERSION_MINOR := $(word 2,${NUMPY_VERSION_FULL})

## Same reason as python patch level. 
NUMPY_VERSION_PATCH := $(word 3,${NUMPY_VERSION_FULL})

### Check for minimum python + numpy versions. In theory, I should also check
### that *any* python and numpy are available but that seems too much effort
MIN_PYTHON_MAJOR := 2
MIN_PYTHON_MINOR := 6

MIN_NUMPY_MAJOR  := 1
MIN_NUMPY_MINOR  := 7

PYTHON_AVAIL := $(shell [ $(PYTHON_VERSION_MAJOR) -gt $(MIN_PYTHON_MAJOR) -o \( $(PYTHON_VERSION_MAJOR) -eq $(MIN_PYTHON_MAJOR) -a $(PYTHON_VERSION_MINOR) -ge $(MIN_PYTHON_MINOR) \) ] && echo true)
NUMPY_AVAIL  := $(shell [ $(NUMPY_VERSION_MAJOR) -gt $(MIN_NUMPY_MAJOR) -o \( $(NUMPY_VERSION_MAJOR) -eq $(MIN_NUMPY_MAJOR) -a $(NUMPY_VERSION_MINOR) -ge $(MIN_NUMPY_MINOR) \) ] && echo true)

ifneq ($(PYTHON_AVAIL),true)
$(warning Found python version $(PYTHON_VERSION_MAJOR).$(PYTHON_VERSION_MINOR).$(PYTHON_VERSION_PATCH) but minimum required python is $(MIN_PYTHON_MAJOR).$(MIN_PYTHON_MINOR))
COMPILE_PYTHON_EXT := 0
endif

ifneq ($(NUMPY_AVAIL),true)
$(warning Found NUMPY version $(NUMPY_VERSION_MAJOR).$(NUMPY_VERSION_MINOR).$(NUMPY_VERSION_PATCH) but minimum required numpy is $(MIN_NUMPY_MAJOR).$(MIN_NUMPY_MINOR))
COMPILE_PYTHON_EXT := 0
endif

ifeq ($(PYTHON_VERSION_MAJOR), 2)
PYTHON_CONFIG_EXE:=python-config
else
PYTHON_CONFIG_EXE:=python3-config
endif
PYTHON_CFLAGS := $(shell $(PYTHON_CONFIG_EXE) --includes) $(shell python -c "from __future__ import print_function; import numpy; print('-isystem' + numpy.__path__[0] + '/core/include/numpy/')")
PYTHON_LIBDIR := $(shell $(PYTHON_CONFIG_EXE) --prefix)/lib
PYTHON_LIBS   := $(shell $(PYTHON_CONFIG_EXE) --libs)
PYTHON_LINK   := -L$(PYTHON_LIBDIR) $(PYTHON_LIBS) -Xlinker -rpath -Xlinker $(PYTHON_LIBDIR)
PYTHON_LIB_BASE := $(strip $(subst -l,lib, $(filter -lpython%,$(PYTHON_LIBS))))

### Check if conda is being used on OSX - then we need to fix python link libraries
UNAME := $(shell uname)
FIX_PYTHON_LINK := 0
ifeq ($(UNAME), Darwin)
## use the clang assembler instead of GNU assembler
## http://stackoverflow.com/questions/10327939/erroring-on-no-such-instruction-while-assembling-project-on-mac-os-x-lion
ifeq (gcc,$(findstring gcc,$(CC)))
  CFLAGS += -Wa,-q
endif
PATH_TO_PYTHON := $(shell which python)
ifeq (conda, $(findstring conda, $(PATH_TO_PYTHON)))
FIX_PYTHON_LINK := 1
PYTHON_LINK := $(filter-out -framework, $(PYTHON_LINK))
PYTHON_LINK := $(filter-out CoreFoundation, $(PYTHON_LINK))
PYTHON_LINK += -dynamiclib -Wl,-single_module -undefined dynamic_lookup -Wl,-compatibility_version,$(VERSION) -Wl,-current_version,$(VERSION) 
endif


### Another check for stack-size. travis ci chokes on this with gcc
# comma := ,
# PYTHON_LINK := $(filter-out -Wl$(comma)-stack_size$(comma)1000000$(comma), $(PYTHON_LINK))
# PYTHON_LINK := $(filter-out -Wl$(comma)-stack_size$(comma)1000000$(comma), $(PYTHON_LINK))
# PYTHON_LINK := $(filter-out -stack_size$(comma)1000000$(comma), $(PYTHON_LINK))
endif




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
  CFLAGS += -xhost -opt-prefetch -opt-prefetch-distance=16 #-vec-report6  
  ifeq (USE_OMP,$(findstring USE_OMP,$(OPT)))
		CFLAGS += -openmp
		CLINK  += -openmp 
  endif
else

  ### compiler specific flags for gcc
  ifeq (gcc,$(findstring gcc,$(CC)))
		CFLAGS += -ftree-vectorize -funroll-loops -fprefetch-loop-arrays --param simultaneous-prefetches=4 #-ftree-vectorizer-verbose=6 -fopt-info-vec-missed #-fprofile-use -fprofile-correction #-fprofile-generate
    ifeq (USE_OMP,$(findstring USE_OMP,$(OPT)))
			CFLAGS += -fopenmp
			CLINK  += -fopenmp
    endif
  endif

  ### compiler specific flags for clang
  ifeq (clang,$(findstring clang,$(CC)))
		CFLAGS += -funroll-loops
		ifeq (USE_OMP,$(findstring USE_OMP,$(OPT)))
      CLANG_VERSION:=$(shell $(CC) -dumpversion 2>&1)
      ifeq ($(CLANG_OMP_AVAIL),1)
			  CFLAGS += -fopenmp
			  CLINK  += -fopenmp=libomp
      else
        $(warning clang does not support OpenMP - please use gcc/icc for compiling with openmp. Removing USE_OMP from compile options)
        OPT:=$(filter-out -DUSE_OMP,$(OPT))
      endif
    endif
  endif

  ifeq (USE_AVX,$(findstring USE_AVX,$(OPT)))
    CFLAGS  +=  -mavx -mpopcnt
  endif

  #### common options for gcc and clang
  CFLAGS  += -march=native
	CFLAGS  += -Wformat=2  -Wpacked  -Wnested-externs -Wpointer-arith  -Wredundant-decls  -Wfloat-equal -Wcast-qual  
  CFLAGS  +=  -Wcast-align -Wmissing-declarations -Wmissing-prototypes  -Wnested-externs -Wstrict-prototypes  #-D_POSIX_C_SOURCE=2 -Wpadded -Wconversion
  CLINK += -lm
endif

ifeq (OUTPUT_THETAAVG,$(findstring OUTPUT_THETAAVG,$(OPT)))
  ifneq (DOUBLE_PREC,$(findstring DOUBLE_PREC,$(OPT)))
    $(error DOUBLE_PREC must be enabled with OUTPUT_THETAAVG -- loss of precision will give you incorrect results for the outer bins (>=20-30 million pairs))
  endif
  ifeq (USE_AVX,$(findstring USE_AVX,$(OPT)))
     ifneq (icc,$(findstring icc,$(CC)))
        $(warning WARNING: OUTPUT_THETAAVG with AVX capabilties is slow with gcc (disables AVX essentially) with gcc. Try to use icc if available)
     endif
  endif
endif

ifeq (FAST_DIVIDE,$(findstring FAST_DIVIDE,$(OPT)))
  ifneq (USE_AVX,$(findstring USE_AVX,$(OPT)))
    $(warning Makefile option FAST_DIVIDE will not do anything unless USE_AVX is set)
  endif
endif


### The following sections are currently not relevant for the Corrfunc package
### but I do not want to have to figure this out again! 
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

