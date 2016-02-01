include ../../theory.options
ifneq (clean,$(findstring clean,$(MAKECMDGOALS)))
  include ../../common.mk
endif


target=DD 
SRC1   = DD.c countpairs.c ../../utils/gridlink.c ../../utils/utils.c ../../utils/progressbar.c ../../io/ftread.c ../../io/io.c 
OBJS1  = $(SRC1:.c=.o)
INCL   = countpairs.h ../../io/ftread.h ../../io/io.h ../../utils/utils.h ../../utils/gridlink.h \
		../../utils/function_precision.h ../../utils/cellarray.h ../../utils/avx_calls.h \
		../../utils/defs.h ../../utils/progressbar.h 

all: $(target) $(SRC1) $(INCL) ../../theory.options ../../common.mk Makefile 

$(target): $(OBJS1) $(INCL) ../../theory.options ../../common.mk Makefile 
	$(CC) $(OBJS1) $(CLINK) -o $@

%.o: %.c $(INCL) ../../theory.options ../../common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

SRC2   = countpairs.c ../../utils/gridlink.c ../../utils/utils.c ../../utils/progressbar.c 
OBJS2  = $(SRC2:.c=.o)

libs: lib
lib:  libcountpairs.a | ../../lib

libcountpairs.a: $(OBJS2) ../../theory.options ../../common.mk Makefile | ../../lib ../../include
	ar rcs $@ $(OBJS2)
	cp -p $@ ../../lib/
	sed -e "s/DOUBLE/$(VECTOR_TYPE)/g" countpairs.h > ../../include/countpairs.h

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
	$(RM) ../../lib/libcountpairs.a ../../include/countpairs.h
	cd ../../bin && $(RM) $(target)

clean:
	$(RM) $(target) $(OBJS1) $(OBJS2) libcountpairs.a

clena: clean
celan: clean


tests: ../tests/tests_periodic ../tests/tests_nonperiodic 
	$(MAKE) -C ../tests	xi_of_r

../tests/tests_periodic:
	$(MAKE) -C ../tests

../tests/tests_nonperiodic:
	$(MAKE) -C ../tests

