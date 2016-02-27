include ../../mocks.options
ifneq (clean,$(findstring clean,$(MAKECMDGOALS)))
  include ../../common.mk
endif

target=DDtheta_mocks wtheta

# ifeq (USE_MKL,$(findstring USE_MKL,$(OPT)))
# MKL_CHECK=mkl
# endif

INCL   = countpairs_theta_mocks.h ../../io/io.h ../../io/ftread.h ../../io/io.h ../../utils/utils.h ../../utils/gridlink_mocks.h \
        ../../utils/function_precision.h ../../utils/cellarray_mocks.h ../../utils/avx_calls.h \
        ../../utils/defs.h ../../utils/progressbar.h

SRC1= DDtheta_mocks.c  countpairs_theta_mocks.c ../../utils/gridlink_mocks.c \
			../../utils/utils.c ../../io/io.c ../../io/ftread.c ../../utils/progressbar.c
OBJS1=$(SRC1:.c=.o) 

all: $(target) $(SRC1) $(INCL) ../../mocks.options ../../common.mk Makefile

DDtheta_mocks: $(OBJS1) $(INCL) ../../mocks.options ../../common.mk Makefile
	$(CC) $(OBJS1) $(CLINK) -o $@ 

SRC2=wtheta.c ../../utils/utils.c
OBJS2=$(SRC2:.c=.o)
wtheta: $(OBJS2) $(INCL) ../../mocks.options ../../common.mk Makefile
	$(CC) $(OBJS2) $(CLINK) -o $@ 

%.o: %.c $(INCL) ../../mocks.options ../../common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

# countpairs_theta_brute_force.o: countpairs_theta_brute_force.c $(INCL) ../../mocks.options ../../common.mk Makefile $(MKL_CHECK)
# 	$(CC) $(OPT) $(CFLAGS) $(BLAS_INCLUDE) $(INCLUDE) -c $< -o $@

# mkl:
# ifndef MKLROOT
# 	$(error $(ccred) You have requested MKL but MKLROOT environment variable is not defined $(ccreset))
# endif

libs:lib
lib: libcountpairs_theta_mocks.a | ../../lib

OBJS3=countpairs_theta_mocks.o ../../utils/gridlink_mocks.o ../../utils/utils.o ../../utils/progressbar.o
libcountpairs_theta_mocks.a: $(OBJS3) $(INCL) ../../mocks.options ../../common.mk Makefile | ../../lib ../../include
	ar rcs $@ $(OBJS3)
	cp -p $@ ../../lib/
	sed -e "s/DOUBLE/$(VECTOR_TYPE)/g" countpairs_theta_mocks.h > ../../include/countpairs_theta_mocks.h

../../lib ../../bin ../../include:
	mkdir -p $@

install: $(target) | ../../bin
	cp -p $(target) ../../bin/

.PHONY: clean clena celan distclean realclean install lib libs tests

tests:
	$(MAKE) -C ../tests wtheta

clean:
	$(RM) $(target) $(OBJS1) $(OBJS2) libcountpairs_theta_mocks.a

clena: clean
celan: clean

distclean:clean | ../../lib ../../include ../../bin
	$(RM) ../../lib/libcountpairs_theta_mocks.a ../../include/countpairs_theta_mocks.h
	cd ../../bin && $(RM) $(target)

realclean: distclean
distclena: distclean
realclena: distclean

