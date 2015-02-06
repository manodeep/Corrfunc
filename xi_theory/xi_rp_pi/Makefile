include ../../theory.options ../../common.mk

target=DDrppi wprp
all: $(target) 


SRC1= DDrppi.c countpairs_rp_pi.c ../../utils/utils.c ../../utils/gridlink.c ../../utils/progressbar.c ../../io/io.c ../../io/ftread.c 
OBJS1  = $(SRC1:.c=.o)
INCL   = countpairs_rp_pi.h ../../io/ftread.h ../../io/io.h ../../utils/utils.h ../../utils/gridlink.h \
		../../utils/function_precision.h ../../utils/cellarray.h ../../utils/avx_calls.h \
		../../utils/defs.h ../../utils/progressbar.h 


DDrppi: $(OBJS1) $(INCL) ../../common.mk Makefile
	$(CC) $(CLINK) $(OBJS1) -o $@

%.o: %.c $(INCL) ../../common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

OBJS2= wprp.o ../../utils/utils.o
wprp: $(OBJS2) Makefile 
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $(OBJS2) $(CLINK)

lib: libcountpairs_rp_pi.a | ../../lib

OBJS3=countpairs_rp_pi.o ../../utils/gridlink.o ../../utils/utils.o ../../utils/progressbar.o 
libcountpairs_rp_pi.a: $(OBJS3) | ../../lib ../../include
	ar rcs $@ $(OBJS3)
	cp -p $@ ../../lib/
	sed -e "s/DOUBLE/$(VECTOR_TYPE)/g" countpairs_rp_pi.h > ../../include/countpairs_rp_pi.h

../../lib:
	mkdir -p ../../lib

../../bin:
	mkdir -p ../../bin

../../include:
	mkdir -p ../../include

install: $(target) | ../../bin 
	cp -p $(target) ../../bin/

.PHONY: clean clena celan install lib tests

clean:
	$(RM) $(target) $(OBJS1) $(OBJS2) $(OBJS3) libcountpairs_rp_pi.a

clena: clean

celan: celan

distclean:clean | ../../lib ../../include ../../bin
	$(RM) ../../lib/libcountpairs_rp_pi.a ../../include/countpairs_rp_pi.h
	cd ../../bin && $(RM) $(target)

tests: ../tests/tests_periodic ../tests/tests_nonperiodic 
	$(MAKE) -C ../tests xi_rp_pi

../tests/tests_periodic:
	$(MAKE) -C ../tests

../tests/tests_nonperiodic:
	$(MAKE) -C ../tests


