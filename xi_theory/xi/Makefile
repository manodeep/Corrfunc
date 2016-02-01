include ../../theory.options
ifneq (clean,$(findstring clean,$(MAKECMDGOALS)))
  include ../../common.mk
endif

target=xi 
SRC1   = xi.c countpairs_xi.c ../../utils/gridlink.c ../../utils/utils.c ../../utils/progressbar.c ../../io/ftread.c ../../io/io.c 
OBJS1  = $(SRC1:.c=.o)
INCL   = countpairs_xi.h ../../io/ftread.h ../../io/io.h ../../utils/utils.h ../../utils/gridlink.h \
		../../utils/function_precision.h ../../utils/cellarray.h ../../utils/avx_calls.h \
		../../utils/defs.h ../../utils/progressbar.h 

all: $(target) $(SRC1) $(INCL) ../../theory.options ../../common.mk Makefile 

$(target): $(OBJS1) $(INCL) ../../common.mk ../../theory.options Makefile 
	$(CC) $(OBJS1) $(CLINK) -o $@

%.o: %.c $(INCL) ../../theory.options ../../common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

SRC2   = countpairs_xi.c ../../utils/gridlink.c ../../utils/utils.c ../../utils/progressbar.c 
OBJS2  = $(SRC2:.c=.o)

libs: lib
lib:  libcountpairs_xi.a | ../../lib ../../include 

libcountpairs_xi.a: $(OBJS2) ../../common.mk ../../theory.options Makefile | ../../lib ../../include
	ar rcs $@ $(OBJS2)
	cp -p $@ ../../lib/
	sed -e "s/DOUBLE/$(VECTOR_TYPE)/g" countpairs_xi.h > ../../include/countpairs_xi.h 

.PHONY: clean clena celan install lib tests distclean

../../lib:
	mkdir -p ../../lib

../../include:
	mkdir -p ../../include

../../bin:
	mkdir -p ../../bin

install: $(target) | ../../bin
	cp -p $(target) ../../bin/

distclean: clean | ../../lib ../../include ../../bin
	$(RM) ../../lib/libcountpairs_xi.a ../../include/countpairs_xi.h
	cd ../../bin && $(RM) $(target)

clean:
	$(RM) $(target) $(OBJS1) $(OBJS2) libcountpairs_xi.a

clena: clean
celan: clean

tests: tests_periodic 

tests_periodic: ../tests/tests_periodic
	$(MAKE) -C ../tests xi

../tests/tests_periodic:
	$(MAKE) -C ../tests

