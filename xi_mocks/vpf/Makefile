include ../../mocks.options ../../common.mk

target=vpf_mocks 
SRC1   = vpf_mocks.c countspheres_mocks.c ../../utils/gridlink.c ../../utils/utils.c ../../utils/progressbar.c ../../io/ftread.c ../../io/io.c \
				../../utils/set_cosmo_dist.c ../../utils/cosmology_params.c 
OBJS1  = $(SRC1:.c=.o)
INCL   = countspheres_mocks.h ../../io/ftread.h ../../io/io.h ../../utils/utils.h ../../utils/gridlink.h \
				../../utils/function_precision.h ../../utils/cellarray.h ../../utils/avx_calls.h \
				../../utils/defs.h ../../utils/progressbar.h ../../utils/set_cosmo_dist.h ../../utils/cosmology_params.h

all: $(target) $(SRC1) $(INCL)

$(target): $(OBJS1) $(INCL) ../../common.mk Makefile 
	$(CC) $(OBJS1) $(CLINK) $(GSL_LINK) -o $@

%.o: %.c $(INCL) ../../common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

countspheres_mocks.o: countspheres_mocks.c $(INCL) ../../common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) $(GSL_CFLAGS) $(INCLUDE) -c $< -o $@



libs: lib
lib:  libcountspheres_mocks.a | ../../lib

SRC2   = countspheres_mocks.c ../../utils/gridlink.c ../../utils/utils.c ../../utils/progressbar.c \
				../../utils/set_cosmo_dist.c ../../utils/cosmology_params.c

OBJS2  = $(SRC2:.c=.o)
libcountspheres_mocks.a: $(OBJS2) | ../../lib ../../include
	ar rcs $@ $(OBJS2)
	cp -p $@ ../../lib/
	sed -e "s/DOUBLE/$(VECTOR_TYPE)/g" countspheres_mocks.h > ../../include/countspheres_mocks.h 

.PHONY: clean clena celan install lib tests 

../../lib:
	mkdir -p ../../lib

../../include:
	mkdir -p ../../include

../../bin:
	mkdir -p ../../bin

install: $(target) | ../../bin
	cp -p $(target) ../../bin/

distclean: clean
	$(RM) ../../bin/$(target)
	$(RM) ../../include/countspheres_mocks.h
	$(RM) ../../lib/libcountspheres_mocks.a

clean:
	$(RM) $(target) $(OBJS1) $(OBJS2) libcountspheres_mocks.a

clena: clean
celan: clean


tests: ../tests/tests_periodic ../tests/tests_nonperiodic 
	$(MAKE) -C ../tests	vpf

../tests/tests_periodic:
	$(MAKE) -C ../tests

../tests/tests_nonperiodic:
	$(MAKE) -C ../tests

