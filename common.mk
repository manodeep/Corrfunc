### Set the default compiler -- options are icc/gcc/clang.
CC:=gcc

#### Add any compiler specific flags you want
CFLAGS:=

#### Add any compiler specific link flags you want
CLINK:=

### You should NOT edit below this line
DISTNAME:=Corrfunc
MAJOR:=1
MINOR:=1
PATCHLEVEL:=0
VERSION:=$(MAJOR).$(MINOR).$(PATCHLEVEL)

DO_CHECKS := 1
ifeq (clean,$(findstring clean,$(MAKECMDGOALS)))
  DO_CHECKS := 0
endif

ifeq (distclean,$(findstring distclean,$(MAKECMDGOALS)))
  DO_CHECKS := 0
endif

ifeq (realclean,$(findstring realclean,$(MAKECMDGOALS)))
  DO_CHECKS := 0
endif

## Only set everything if the command is not "make clean"
ifeq ($(DO_CHECKS), 1)
  ## First check make version. Versions of make older than 3.80 will crash
  ifneq (3.80,$(firstword $(sort $(MAKE_VERSION) 3.80)))
    ## Order-only attributes were added to make version 3.80
    $(warning $(ccmagenta)Please upgrade $(ccblue)make$(ccreset))
    ifeq ($(UNAME), Darwin)
      $(info $(ccmagenta)on Mac+homebrew, use $(ccgreen)"brew outdated xctool || brew upgrade xctool"$(ccreset))
      $(info $(ccmagenta)Otherwise, install $(ccblue)XCode command-line tools$(ccmagenta) directly: $(ccgreen)"xcode-select --install"$(ccreset))
      $(info $(ccmagenta)This link: $(ccgreen)"http://railsapps.github.io/xcode-command-line-tools.html"$(ccmagenta) has some more details$(ccreset))
    else
      $(info $(ccmagenta)On Linux: Try some variant of $(ccgreen)"sudo apt-get update && sudo apt-get upgrade"(ccreset))
    endif
    $(error $(ccmagenta)Project requires make >= 3.80 to compile.$(ccreset))
  endif
  #end of checks for make. 


  UNAME := $(shell uname)
  ## Colored text output
  ## Taken from: http://stackoverflow.com/questions/24144440/color-highlighting-of-makefile-warnings-and-errors
  ## Except, you have to use "echo -e" on linux and "echo" on Mac
  ECHO_COMMAND := echo -e
  ifeq ($(UNAME), Darwin)
    ECHO_COMMAND := echo
  endif
  ccreset :=$(shell $(ECHO_COMMAND) "\033[0;0m")
  ccred:=$(shell $(ECHO_COMMAND) "\033[0;31m")
  ccmagenta:=$(shell $(ECHO_COMMAND) "\033[0;35m")
  ccgreen:=$(shell $(ECHO_COMMAND) "\033[0;32m")
  ccblue:=$(shell $(ECHO_COMMAND) "\033[0;34m")
  ## end of colored text output


  ## Make clang the default compiler on Mac
  ## But first check for clang-omp, use that if available
  ifeq ($(UNAME), Darwin)
    CLANG_OMP_FOUND := $(shell clang-omp --version 2>/dev/null)
    ifndef CLANG_OMP_FOUND
      CC := clang
    else
      CC := clang-omp
    endif
  endif
  # Check if CPU supports AVX -> this trumps everything. For instance, compiler might
  # support AVX but the cpu might not. Then compilation will work fine but there will
  # be a runtime crash with "Illegal Instruction"
  ifeq ($(UNAME), Darwin)
    # On a MAC, best to use sysctl
    AVX_AVAIL := $(shell sysctl -n machdep.cpu.features 2>/dev/null | grep -o -i AVX | tr '[:lower:]' '[:upper:]')
  else
    # On Linux/Unix, just grep on /proc/cpuinfo
    # There might be multiple cores, so just take the first line
    # (Is it possible that someone has one core that has AVX and another that doesnt?)
    AVX_AVAIL := $(shell grep -o -i AVX /proc/cpuinfo 2>/dev/null | head -n 1 | tr '[:lower:]' '[:upper:]' )
  endif
  REMOVE_AVX :=0
  ifdef AVX_AVAIL
    ifneq ($(AVX_AVAIL) , AVX)
      REMOVE_AVX := 1
    endif
  else
    REMOVE_AVX :=1
  endif

  ifeq ($(REMOVE_AVX), 1)
    $(warning $(ccmagenta) CPU does not seem support AVX instructions. Removing USE_AVX from compile options. $(ccreset))
    OPT:=$(filter-out -DUSE_AVX,$(OPT))
  endif
  # end of checking if CPU supports AVX      


  # Now check if gcc is set to be the compiler but if clang is really under the hood.
  export CC_IS_CLANG ?= -1
  ifeq ($(CC_IS_CLANG), -1)
    CC_VERSION := $(shell $(CC) --version 2>/dev/null)
    ifndef CC_VERSION
      $(error $(ccred)Could not find $$CC = ${CC}$(ccreset))
    endif
    ifeq (clang,$(findstring clang,$(CC_VERSION)))
      export CC_IS_CLANG := 1
    else
      export CC_IS_CLANG := 0
    endif
    export CC_VERSION
	endif
  # Done with checking if clang is underneath gcc

  # CC is set at this point. In case the compiler on Mac is *not* clang under the hood
  # print an info message saying what to do in case of an error
  ifeq ($(UNAME), Darwin)
    ifneq ($(CC_IS_CLANG), 1)
      export CLANG_COMPILER_WARNING_PRINTED ?= 0
	    ifeq ($(CLANG_COMPILER_WARNING_PRINTED), 0)
        $(warning $(ccmagenta)Looks like clang (on Mac) is not set as the compiler. If you run into errors like $(ccred)"no such instruction: `vxorpd %xmm1, %xmm1,%xmm1'"$(ccmagenta), then please use clang as the compiler (directly invoke $(ccgreen)make$(ccmagenta), NOT $(ccred)make CC=gcc$(ccmagenta)) $(ccreset))
        export CLANG_COMPILER_WARNING_PRINTED := 1
      endif
    endif
  endif

  INCLUDE:=-I../../io -I../../utils
  ### The POSIX_SOURCE flag is required to get the definition of strtok_r
  CFLAGS 	+= -Wsign-compare -Wall -Wextra -Wshadow -Wunused -std=c99 -g -m64 -fPIC -D_POSIX_SOURCE -D_DARWIN_C_SOURCE -O3 #-Ofast
  GSL_FOUND := $(shell gsl-config --version)
  ifndef GSL_FOUND
    $(error $(ccred) GSL not found in path - please install GSL before installing $(DISTNAME).$(VERSION) $(ccreset))
  endif
  GSL_CFLAGS := $(shell gsl-config --cflags)
  GSL_LIBDIR := $(shell gsl-config --prefix)/lib
  GSL_LINK   := $(shell gsl-config --libs) -Xlinker -rpath -Xlinker $(GSL_LIBDIR)

  # Check if all progressbar output is to be suppressed
  OUTPUT_PGBAR := 1
	ifeq (SILENT, $(findstring SILENT, $(CFLAGS)))
    OUTPUT_PGBAR := 0
  endif

  ifeq (SILENT, $(findstring SILENT, $(OPT)))
    OUTPUT_PGBAR := 0
  endif
  #end of progressbar checks



  ## Check for conflicting options
  ifeq (OUTPUT_RPAVG,$(findstring OUTPUT_RPAVG,$(OPT)))
    ifneq (DOUBLE_PREC,$(findstring DOUBLE_PREC,$(OPT)))
      $(error $(ccred) DOUBLE_PREC must be enabled with OUTPUT_RPAVG -- loss of precision will give you incorrect results for the outer bins (>=20-30 million pairs) $(ccreset))
    endif
  endif

  ifeq (OUTPUT_THETAAVG,$(findstring OUTPUT_THETAAVG,$(OPT)))
    ifneq (DOUBLE_PREC,$(findstring DOUBLE_PREC,$(OPT)))
      $(error $(ccred) DOUBLE_PREC must be enabled with OUTPUT_THETAAVG -- loss of precision will give you incorrect results for the outer bins (>=20-30 million pairs) $(ccreset))
    endif
  endif

  ifeq (FAST_DIVIDE,$(findstring FAST_DIVIDE,$(OPT)))
    ifneq (USE_AVX,$(findstring USE_AVX,$(OPT)))
      $(warning $(ccmagenta) Makefile option FAST_DIVIDE will not do anything unless USE_AVX is set $(ccreset))
    endif
  endif
  ## done with check for conflicting options

  ## Set the vector type appropriately -> sed will use this to generate the header files correctly
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
    endif ##openmp with icc
  else ## not icc -> gcc or clang follow

    ## Warning that w(theta) with OUTPUT_THETAAVG is very slow without icc
    ## Someday I am going to fix that by linking with MKL 
    ifeq (USE_AVX,$(findstring USE_AVX,$(OPT)))
      ifeq (OUTPUT_THETAAVG,$(findstring OUTPUT_THETAAVG,$(OPT)))
        $(warning $(ccmagenta) WARNING: OUTPUT_THETAAVG with AVX capabilties is slow with gcc/clang (disables AVX essentially) with gcc/clang. Try to use icc if available $(ccreset))
      endif
	  endif

    ### GCC is slightly more complicated. CC might be called gcc but it might be clang underneath
    ### compiler specific flags for gcc
    ifneq ($(CC_IS_CLANG), 1)
      ## Real gcc here
      ifeq (gcc,$(findstring gcc,$(CC)))
        CFLAGS += -ftree-vectorize -funroll-loops -fprefetch-loop-arrays --param simultaneous-prefetches=4 #-ftree-vectorizer-verbose=6 -fopt-info-vec-missed #-fprofile-use -fprofile-correction #-fprofile-generate
        # Use the clang assembler on Mac.
        ifeq ($(UNAME), Darwin)
          CFLAGS += -Wa,-q
          export CLANG_ASM_WARNING_PRINTED ?= 0
	        ifeq ($(CLANG_ASM_WARNING_PRINTED), 0)
            $(warning $(ccmagenta) WARNING: gcc on Mac does not support intrinsics. Attempting to use the clang assembler $(ccreset))
            $(warning $(ccmagenta) If you see the error message $(ccred) "/opt/local/bin/as: assembler (/opt/local/bin/clang) not installed" $(ccmagenta) then try the following fix $(ccreset))
            $(warning $(ccmagenta) Either install clang ($(ccgreen)for Macports use, "sudo port install clang-3.8"$(ccmagenta)) or remove option $(ccgreen)"USE_AVX"$(ccmagenta) from both the $(ccgreen)"theory.options"$(ccmagenta) and $(ccgreen)"mocks.options"$(ccmagenta) files$(ccreset))
            export CLANG_ASM_WARNING_PRINTED := 1
          endif # warning printed
        endif

        ifeq (USE_OMP,$(findstring USE_OMP,$(OPT)))
          CFLAGS += -fopenmp
          CLINK  += -fopenmp
        endif #openmp with gcc
      endif #gcc findstring
    else ##CC is clang
      ### compiler specific flags for clang
      CLANG_OMP_AVAIL := false
      export APPLE_CLANG := 0
      ifeq (USE_OMP,$(findstring USE_OMP,$(OPT)))
        ifeq (clang-omp,$(findstring clang-omp,$(CC)))
          CLANG_OMP_AVAIL:=true
          CFLAGS += -fopenmp
          CLINK  += -liomp5
        else
          # Apple clang/gcc does not support OpenMP
          ifeq (Apple, $(findstring Apple, $(CC_VERSION)))
            CLANG_OMP_AVAIL:= false
            $(warning $(ccmagenta)Compiler is Apple clang and does not support OpenMP$(ccreset))
            $(info $(ccmagenta)If you want OpenMP support, please install clang with OpenMP support$(ccreset))
            $(info $(ccmagenta)For homebrew, use $(ccgreen)"brew update && (brew outdated xctool || brew upgrade xctool) && brew tap homebrew/versions && brew install clang-omp"$(ccreset))
            $(info $(ccmagenta)For Macports, use $(ccgreen)"sudo port install clang-3.8 +assertions +debug + openmp"$(ccreset))
            export APPLE_CLANG := 1
          else
            ## Need to do a version check clang >= 3.7 supports OpenMP. If it is Apple clang, then it doesn't support OpenMP.
            ## All of the version checks go here. If OpenMP is supported, update CLANG_OMP_AVAIL to 1.
            CLANG_VERSION_FULL := $(shell $(CC) --version |  egrep -o 'version (.*)' | awk "{print \$$2}")
            CLANG_VERSION_FULL :=  $(subst ., ,$(CLANG_VERSION_FULL))
            CLANG_VERSION_MAJOR := $(word 1,${CLANG_VERSION_FULL})
            CLANG_VERSION_MINOR := $(word 2,${CLANG_VERSION_FULL})
            CLANG_MAJOR_MIN_OPENMP := 3
            CLANG_MINOR_MIN_OPENMP := 7
            CLANG_OMP_AVAIL := $(shell [ $(CLANG_VERSION_MAJOR) -gt $(CLANG_MAJOR_MIN_OPENMP) -o \( $(CLANG_VERSION_MAJOR) -eq $(CLANG_MAJOR_MIN_OPENMP) -a $(CLANG_VERSION_MINOR) -ge $(CLANG_MINOR_MIN_OPENMP) \) ] && echo true)
            CLANG_IS_38 := $(shell [ $(CLANG_VERSION_MAJOR) -eq 3 -a $(CLANG_VERSION_MINOR) -eq 8  ] && echo true)
            CFLAGS += -fopenmp=libomp
            CLINK  += -fopenmp=libomp
          endif #Apple check
        endif  #clang-omp check

        ifeq ($(CLANG_OMP_AVAIL),true)
          ifeq ($(APPLE_CLANG),0)
            ifeq ($(UNAME), Darwin)
              $(info $(ccmagenta)Enabling OpenMP with clang.$(ccreset))
              CLANG_LD_ERROR := "dyld: Library not loaded: @rpath/libLLVM.dylib\nReferenced from: /opt/local/libexec/llvm-3.8/lib/libLTO.dylib\nReason: image not found\n"
              export CLANG_LD_WARNING_PRINTED ?= 0
              ifeq ($(CLANG_LD_WARNING_PRINTED), 0)
                ifeq ($(CLANG_IS_38), true)
                  $(warning With $(ccgreen)clang-3.8$(ccreset), You might see this $(ccred)$(CLANG_LD_ERROR)$(ccmagenta) error with the final linking step.$(ccreset))
                  $(info $(ccmagenta)Use the following to fix the issue $(ccgreen) "sudo install_name_tool -change @executable_path/../lib/libLTO.dylib @rpath/../lib/libLTO.dylib /opt/local/libexec/ld64/ld-latest"$(ccreset))
                  $(info $(ccmagenta)You can see the bug report here $(ccgreen)"https://trac.macports.org/ticket/50853"$(ccreset))
                  export CLANG_LD_WARNING_PRINTED := 1
                endif #clang-3.8
              endif #clang warning printed
            endif #Darwin
          endif #Apple clang. If at some point Apple clang supports OpenMP, then there will need to be an else above this endif. 
        else
          # I dislike being warned multiple times but the compiler warning will not
          # be visible if the entire codebase is being compiled. 
          # export WARNING_PRINTED ?= 0
	        # ifeq ($(WARNING_PRINTED), 0)
          $(warning $(ccmagenta) $$CC = ${CC} does not support OpenMP - please use gcc/icc for compiling with openmp. Removing USE_OMP from compile options. $(ccreset))
          infovar := "OPT:=$$(filter-out -DUSE_OMP,$$(OPT))"
          $(info $(ccmagenta)If you are sure your version of clang ($(ccblue) must be >= 3.7, NOT Apple clang $(ccmagenta)) does support OpenMP, then comment out the line $(ccred) $(infovar) $(ccmagenta) in the file $(ccgreen)"common.mk"$(ccreset))
          $(info $(ccmagenta)You might have to add in the include path (path to omp.h) to CFLAGS and the runtime library path to CLINK at the top of $(ccgreen)"common.mk"$(ccreset))

          # comment out the following line if your version of clang definitely supports OpenMP
          OPT:=$(filter-out -DUSE_OMP,$(OPT))
          # export WARNING_PRINTED := 1
         endif # CLANG_OMP_AVAIL is not 1
      endif # USE_OMP
    endif # CC is clang

    #### common options for gcc and clang
    ifeq (USE_AVX,$(findstring USE_AVX,$(OPT)))
      CFLAGS  +=  -mavx -mpopcnt
    endif

    CFLAGS += -funroll-loops
    CFLAGS += -march=native -fno-strict-aliasing
    CFLAGS += -Wformat=2  -Wpacked  -Wnested-externs -Wpointer-arith  -Wredundant-decls  -Wfloat-equal -Wcast-qual
    CFLAGS +=  -Wcast-align -Wmissing-declarations -Wmissing-prototypes  -Wnested-externs -Wstrict-prototypes  #-D_POSIX_C_SOURCE=2 -Wpadded -Wconversion
    CLINK += -lm 
  endif #not icc


  # All of the python/numpy checks follow
  export PYTHON_CHECKED ?= 0
  ifeq ($(PYTHON_CHECKED), 0)
    export COMPILE_PYTHON_EXT := 1
    export PYTHON_VERSION_FULL := $(wordlist 2,4,$(subst ., ,$(shell python --version 2>&1)))
    export PYTHON_VERSION_MAJOR := $(word 1,${PYTHON_VERSION_FULL})
    export PYTHON_VERSION_MINOR := $(word 2,${PYTHON_VERSION_FULL})

    ## I only need this so that I can print out the full python version (correctly)
    ## in case of error
    PYTHON_VERSION_PATCH := $(word 3,${PYTHON_VERSION_FULL})

    ## Check numpy version
    export NUMPY_VERSION_FULL :=  $(wordlist 1,3,$(subst ., ,$(shell python -c "from __future__ import print_function; import numpy; print(numpy.__version__)")))
    export NUMPY_VERSION_MAJOR := $(word 1,${NUMPY_VERSION_FULL})
    export NUMPY_VERSION_MINOR := $(word 2,${NUMPY_VERSION_FULL})

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
      $(warning $(ccmagenta) Found python version $(PYTHON_VERSION_MAJOR).$(PYTHON_VERSION_MINOR).$(PYTHON_VERSION_PATCH) but minimum required python is $(MIN_PYTHON_MAJOR).$(MIN_PYTHON_MINOR) $(ccreset))
	    COMPILE_PYTHON_EXT := 0
    endif

    ifneq ($(NUMPY_AVAIL),true)
      $(warning $(ccmagenta) Found NUMPY version $(NUMPY_VERSION_MAJOR).$(NUMPY_VERSION_MINOR).$(NUMPY_VERSION_PATCH) but minimum required numpy is $(MIN_NUMPY_MAJOR).$(MIN_NUMPY_MINOR) $(ccreset))
      COMPILE_PYTHON_EXT := 0
    endif

    ifeq ($(PYTHON_VERSION_MAJOR), 2)
      PYTHON_CONFIG_EXE:=python-config
    else
	    PYTHON_CONFIG_EXE:=python3-config
    endif
    export PYTHON_CFLAGS := $(shell $(PYTHON_CONFIG_EXE) --includes) $(shell python -c "from __future__ import print_function; import numpy; print('-isystem' + numpy.__path__[0] + '/core/include/numpy/')")
    export PYTHON_LIBDIR := $(shell $(PYTHON_CONFIG_EXE) --prefix)/lib
    export PYTHON_LIBS   := $(shell $(PYTHON_CONFIG_EXE) --libs)
    export PYTHON_LINK   := -L$(PYTHON_LIBDIR) $(PYTHON_LIBS) -Xlinker -rpath -Xlinker $(PYTHON_LIBDIR)
    export PYTHON_LIB_BASE := $(strip $(subst -l,lib, $(filter -lpython%,$(PYTHON_LIBS))))

    ### Check if conda is being used on OSX - then we need to fix python link libraries
    export FIX_PYTHON_LINK := 0
    ifeq ($(UNAME), Darwin)
      PATH_TO_PYTHON := $(shell which python)
      ifeq (conda, $(findstring conda, $(PATH_TO_PYTHON)))
	      FIX_PYTHON_LINK := 1
      endif
      PYTHON_LINK := $(filter-out -framework, $(PYTHON_LINK))
      PYTHON_LINK := $(filter-out -ldl, $(PYTHON_LINK))
      PYTHON_LINK := $(filter-out CoreFoundation, $(PYTHON_LINK))
      PYTHON_LINK += -dynamiclib -Wl,-compatibility_version,$(MAJOR).$(MINOR) -Wl,-current_version,$(VERSION)
      PYTHON_LINK += -headerpad_max_install_names

      ### Another check for stack-size. travis ci chokes on this with gcc
      # comma := ,
      # PYTHON_LINK := $(filter-out -Wl$(comma)-stack_size$(comma)1000000$(comma), $(PYTHON_LINK))
      # PYTHON_LINK := $(filter-out -Wl$(comma)-stack_size$(comma)1000000$(comma), $(PYTHON_LINK))
      # PYTHON_LINK := $(filter-out -stack_size$(comma)1000000$(comma), $(PYTHON_LINK))
    endif #Darwin checks
    export PYTHON_CHECKED:=1
  endif
  ### Done with python checks


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
endif ## make clean is not under effect (the if condition starts all the way at the top with variable DO_CHECKS)
