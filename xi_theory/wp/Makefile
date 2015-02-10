include ../../theory.options ../../common.mk

target=wp
SRC1= wp.c countpairs_wp.c ../../utils/utils.c ../../utils/gridlink.c ../../utils/progressbar.c ../../io/io.c ../../io/ftread.c 
OBJS1  = $(SRC1:.c=.o)
INCL   = countpairs_wp.h ../../io/ftread.h ../../io/io.h ../../utils/utils.h ../../utils/gridlink.h \
		../../utils/function_precision.h ../../utils/cellarray.h ../../utils/avx_calls.h \
		../../utils/defs.h ../../utils/sglib.h ../../utils/progressbar.h 

all: $(target) $(SRC1) $(INCL)

$(target): $(OBJS1) $(INCL) ../../theory.options ../../common.mk Makefile 
	$(CC) $(CLINK) $(OBJS1) -o $@

%.o: %.c $(INCL) ../../theory.options ../../common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

OBJS2=countpairs_wp.o ../../utils/gridlink.o ../../utils/utils.o ../../utils/progressbar.o 

libs: lib
lib: libcountpairs_wp.a | ../../lib ../../include

libcountpairs_wp.a: $(OBJS2) ../../theory.options ../../common.mk Makefile | ../../lib ../../include
	ar rcs $@ $(OBJS2)
	cp -p $@ ../../lib/
	sed -e "s/DOUBLE/$(VECTOR_TYPE)/g" countpairs_wp.h > ../../include/countpairs_wp.h

.PHONY: clean clena celan install lib tests distclean

../../lib:
	mkdir -p ../../lib

../../bin:
	mkdir -p ../../bin

../../include:
	mkdir -p ../../include

install: $(target) | ../../bin
	cp -p $(target) ../../bin/

distclean: clean | ../../lib ../../include ../../bin
	$(RM) ../../lib/libcountpairs_wp.a ../../include/countpairs_wp.h
	cd ../../bin && $(RM) $(target)

clean:
	$(RM) $(target) $(OBJS1) $(OBJS2) libcountpairs_wp.a

clena: clean

celan: celan

tests: tests_periodic 

tests_periodic: ../tests/tests_periodic
	$(MAKE) -C ../tests wp

../tests/tests_periodic:
	$(MAKE) -C ../tests


