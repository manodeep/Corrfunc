include ../../mocks.options ../../common.mk

# ifeq (USE_MKL,$(findstring USE_MKL,$(OPT)))
# MKL_CHECK=mkl
# endif

targets := DDtheta_mocks wtheta
SRC1   := DDtheta_mocks.c  countpairs_theta_mocks.c ../../utils/gridlink_mocks.c \
			../../utils/utils.c ../../io/io.c ../../io/ftread.c
INCL   := countpairs_theta_mocks.h ../../io/io.h ../../io/ftread.h ../../io/io.h ../../utils/utils.h ../../utils/gridlink_mocks.h \
        ../../utils/function_precision.h ../../utils/cellarray_mocks.h ../../utils/avx_calls.h \
        ../../utils/defs.h

ifeq ($(OUTPUT_PGBAR), 1)
SRC1 += ../../utils/progressbar.c
INCL += ../../utils/progressbar.h 
endif

OBJS1=$(SRC1:.c=.o) 

all: $(targets) $(SRC1) $(INCL) ../../mocks.options ../../common.mk Makefile

DDtheta_mocks: $(OBJS1) $(INCL) ../../mocks.options ../../common.mk Makefile
	$(CC) $(OBJS1) $(CLINK) -o $@ 

SRC2  :=wtheta.c ../../utils/utils.c
OBJS2 :=$(SRC2:.c=.o)
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

OBJS3 := countpairs_theta_mocks.o ../../utils/gridlink_mocks.o ../../utils/utils.o ../../utils/progressbar.o
libcountpairs_theta_mocks.a: $(OBJS3) $(INCL) ../../mocks.options ../../common.mk Makefile | ../../lib ../../include
	ar rcs $@ $(OBJS3)
	cp -p $@ ../../lib/
	sed -e "s/DOUBLE/$(VECTOR_TYPE)/g" countpairs_theta_mocks.h > ../../include/countpairs_theta_mocks.h

../../lib ../../bin ../../include:
	mkdir -p $@

install: ../../bin/DDtheta_mocks ../../bin/wtheta | ../../bin

../../bin/DDtheta_mocks:DDtheta_mocks
	cp -p $< $@

../../bin/wtheta: wtheta
	cp -p $< $@

.PHONY: clean clena celan distclean realclean install lib libs tests

tests:
	$(MAKE) -C ../tests wtheta

clean:
	$(RM) $(targets) $(OBJS1) $(OBJS2) libcountpairs_theta_mocks.a

clena: clean
celan: clean

distclean:clean | ../../lib ../../include ../../bin
	$(RM) ../../lib/libcountpairs_theta_mocks.a ../../include/countpairs_theta_mocks.h
	cd ../../bin && $(RM) $(targets)

realclean: distclean
distclena: distclean
realclena: distclean

