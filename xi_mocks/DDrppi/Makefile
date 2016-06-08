include ../../mocks.options ../../common.mk

targets := DDrppi_mocks wprp_mocks
SRC1   := DDrppi_mocks.c  countpairs_rp_pi_mocks.c ../../utils/gridlink_mocks.c \
          ../../utils/set_cosmo_dist.c ../../utils/cosmology_params.c \
			    ../../utils/utils.c ../../io/io.c ../../io/ftread.c
INCL   := countpairs_rp_pi_mocks.h ../../io/io.h ../../io/ftread.h ../../io/io.h  ../../utils/gridlink_mocks.h ../../utils/cellarray_mocks.h \
				../../utils/set_cosmo_dist.h ../../utils/cosmology_params.h  \
				../../utils/utils.h	../../utils/function_precision.h ../../utils/avx_calls.h ../../utils/defs.h

ifeq ($(OUTPUT_PGBAR), 1)
SRC1 += ../../utils/progressbar.c
INCL += ../../utils/progressbar.h 
endif


OBJS1=$(SRC1:.c=.o) 

all: $(targets) $(SRC1) $(INCL) ../../mocks.options ../../common.mk Makefile

DDrppi_mocks: $(OBJS1) $(INCL) ../../mocks.options ../../common.mk Makefile
	$(CC) $(OBJS1) $(CLINK) $(GSL_LINK) -o $@ 

countpairs_rp_pi_mocks.o: countpairs_rp_pi_mocks.c ../../mocks.options ../../common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) $(GSL_CFLAGS) $(INCLUDE) -c $< -o $@

SRC2=wprp_mocks.c ../../utils/utils.c
OBJS2=$(SRC2:.c=.o)
wprp_mocks: $(OBJS2) $(INCL) ../../mocks.options ../../common.mk Makefile
	$(CC) $(OBJS2) $(CLINK) -o $@ 

%.o: %.c $(INCL) ../../mocks.options ../../common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@


../../utils/set_cosmo_dist.o:../../utils/set_cosmo_dist.c ../../mocks.options ../../common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) $(GSL_CFLAGS) $(INCLUDE) -c $< -o $@

libs: lib
lib: libcountpairs_rp_pi_mocks.a | ../../lib

OBJS3 :=countpairs_rp_pi_mocks.o ../../utils/gridlink_mocks.o ../../utils/utils.o ../../utils/progressbar.o \
			../../utils/set_cosmo_dist.o ../../utils/cosmology_params.o 

libcountpairs_rp_pi_mocks.a: $(OBJS3) $(INCL) ../../mocks.options ../../common.mk Makefile | ../../lib ../../include
	ar rcs $@ $(OBJS3)
	cp -p $@ ../../lib/
	sed -e "s/DOUBLE/$(VECTOR_TYPE)/g" countpairs_rp_pi_mocks.h > ../../include/countpairs_rp_pi_mocks.h

../../lib ../../bin ../../include:
	mkdir -p $@

install: ../../bin/DDrppi_mocks ../../bin/wprp_mocks | ../../bin

../../bin/%_mocks: %_mocks
	cp -p $< $@

.PHONY: clean clena celan distclean realclean install tests

tests: 
	$(MAKE) -C ../tests DDrppi

clean:
	$(RM) $(targets) $(OBJS1) $(OBJS2)  libcountpairs_rp_pi_mocks.a

clena: clean
celan: clean

distclean:clean | ../../lib ../../include ../../bin
	$(RM) ../../lib/libcountpairs_rp_pi_mocks.a ../../include/countpairs_rp_pi_mocks.h
	cd ../../bin && $(RM) $(targets)

realclean: distclean
distclena: distclean
realclena: distclean

