include ../common.mk

target=DD 
SRC1   = DD.c countpairs.c ../../utils/gridlink.c ../../utils/utils.c ../../utils/progressbar.c ../../io/ftread.c ../../io/io.c 
OBJS1  = $(SRC1:.c=.o)
INCL   = countpairs.h ../../io/ftread.h ../../io/io.h ../../utils/utils.h ../../utils/gridlink.h \
		../../utils/function_precision.h ../../utils/cellarray.h ../../utils/avx_calls.h \
		../../utils/defs.h ../../utils/progressbar.h 

all: $(target) $(SRC1) $(INCL)

$(target): $(OBJS1) $(INCL) ../common.mk Makefile 
	$(CC) $(CLINK) $(OBJS1) -o $@

%.o: %.c $(INCL) ../common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

SRC2   = countpairs.c ../../utils/gridlink.c ../../utils/utils.c ../../utils/progressbar.c 
OBJS2  = $(SRC2:.c=.o)

lib:  libcountpairs.a | ../lib


libcountpairs.a: $(OBJS2) | ../lib
	ar rcs $@ $(OBJS2)
	cp -p $@ ../lib/
	sed -e "s/DOUBLE/$(VECTOR_TYPE)/g" countpairs.h > ../include/countpairs.h 

.PHONY: clean clena celan install lib tests 

../lib:
	mkdir -p ../lib

../bin:
	mkdir -p ../bin

install: $(target) | ../bin
	cp -p $(target) ../bin/

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

