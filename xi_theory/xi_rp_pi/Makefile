include ../../theory.options ../../common.mk

targets := DDrppi wprp

SRC1   := DDrppi.c countpairs_rp_pi.c countpairs_rp_pi_driver.c ../../utils/utils.c ../../utils/gridlink.c ../../io/io.c ../../io/ftread.c 
INCL   := countpairs_rp_pi.h ../../io/ftread.h ../../io/io.h ../../utils/utils.h ../../utils/gridlink.h \
		../../utils/function_precision.h ../../utils/cellarray.h ../../utils/avx_calls.h \
		../../utils/defs.h countpairs_rp_pi_kernels.c countpairs_rp_pi_driver.h

ifeq ($(OUTPUT_PGBAR), 1)
SRC1 += ../../utils/progressbar.c
INCL += ../../utils/progressbar.h 
endif
OBJS1 := $(SRC1:.c=.o)

all: $(targets) ../../theory.options ../../common.mk Makefile 

DDrppi: $(OBJS1) $(INCL) ../../theory.options ../../common.mk Makefile
	$(CC) $(OBJS1) $(CLINK) -o $@

%.o: %.c $(INCL) ../../theory.options ../../common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

OBJS2 := wprp.o ../../utils/utils.o 
wprp: $(OBJS2) ../../theory.options ../../common.mk Makefile
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $(OBJS2) $(CLINK)

libs: lib
lib: ../../lib/libcountpairs_rp_pi.a | ../../lib

OBJS3 := countpairs_rp_pi.o countpairs_rp_pi_driver.o ../../utils/gridlink.o ../../utils/utils.o
ifeq ($(OUTPUT_PGBAR), 1)
  OBJS3 += ../../utils/progressbar.o
endif

libcountpairs_rp_pi.a: $(OBJS3) $(INCL) ../../theory.options ../../common.mk Makefile | ../../lib ../../include
	ar rcs $@ $(OBJS3)

../../lib/libcountpairs_rp_pi.a: libcountpairs_rp_pi.a
	cp -p $< $@
	sed -e "s/DOUBLE/$(VECTOR_TYPE)/g" countpairs_rp_pi.h > ../../include/countpairs_rp_pi.h

../../lib ../../bin ../../include:
	mkdir -p $@

install: ../../bin/DDrppi ../../bin/wprp | ../../bin 

../../bin/DDrppi: DDrppi
	cp -p $< $@
../../bin/wprp: wprp
	cp -p $< $@

.PHONY: clean clena celan install lib tests

clean:
	$(RM) $(targets) $(OBJS1) $(OBJS2) $(OBJS3) libcountpairs_rp_pi.a

clena: clean
celan: celan

distclean:clean | ../../lib ../../include ../../bin
	$(RM) ../../lib/libcountpairs_rp_pi.a ../../include/countpairs_rp_pi.h
	cd ../../bin && $(RM) $(targets)

tests: $(SRC1)
	$(MAKE) -C ../tests xi_rp_pi



