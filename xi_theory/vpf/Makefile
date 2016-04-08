include ../../theory.options ../../common.mk


target := vpf 
SRC1   := vpf.c countspheres.c ../../utils/gridlink.c ../../utils/utils.c ../../utils/progressbar.c ../../io/ftread.c ../../io/io.c 
OBJS1  := $(SRC1:.c=.o)
INCL   := countspheres.h ../../io/ftread.h ../../io/io.h ../../utils/utils.h ../../utils/gridlink.h \
		../../utils/function_precision.h ../../utils/cellarray.h ../../utils/avx_calls.h \
		../../utils/defs.h ../../utils/progressbar.h 

all: $(target) $(SRC1) $(INCL) ../../theory.options ../../common.mk Makefile 

$(target): $(OBJS1) $(INCL) ../../theory.options ../../common.mk Makefile 
	$(CC) $(OBJS1) $(CLINK) $(GSL_LINK) -o $@

countspheres.o: countspheres.c $(INCL) ../../theory.options ../../common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) $(GSL_CFLAGS) $(INCLUDE) -c $< -o $@

%.o: %.c $(INCL) ../../theory.options ../../common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

SRC2   = countspheres.c ../../utils/gridlink.c ../../utils/utils.c ../../utils/progressbar.c 
OBJS2  = $(SRC2:.c=.o)

libs: lib
lib:  libcountspheres.a | ../../lib ../../include 


libcountspheres.a: $(OBJS2) ../../theory.options ../../common.mk Makefile | ../../lib ../../include
	ar rcs $@ $(OBJS2)
	cp -p $@ ../../lib/
	sed -e "s/DOUBLE/$(VECTOR_TYPE)/g" countspheres.h > ../../include/countspheres.h 

../../bin ../../include ../../lib:
	mkdir -p $@

install: ../../bin/$(target) | ../../bin

../../bin/$(target): $(target)
	cp -p $(target) ../../bin/

distclean: clean | ../../lib ../../include ../../bin
	$(RM) ../../lib/libcountspheres.a
	$(RM) ../../include/countspheres.h
	$(RM) ../../bin/$(target)

clean:
	$(RM) $(target) $(OBJS1) $(OBJS2) libcountspheres.a

clena: clean
celan: clean

.PHONY: clean clena celan install lib tests distclean

tests: ../tests/tests_periodic ../tests/tests_nonperiodic 
	$(MAKE) -C ../tests	vpf

../tests/tests_periodic:
	$(MAKE) -C ../tests

../tests/tests_nonperiodic:
	$(MAKE) -C ../tests

