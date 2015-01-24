include ../common.mk

target=DDrppi wprp
all: $(target) 


SRC1= DDrppi.c countpairs_rp_pi.c ../utils/utils.c ../utils/gridlink.c ../utils/progressbar.c ../io/io.c ../io/ftread.c 
OBJS1  = $(SRC1:.c=.o)
INCL   = countpairs_rp_pi.h ../io/ftread.h ../io/io.h ../utils/utils.h ../utils/gridlink.h \
					../include/function_precision.h ../include/cellarray.h ../include/avx_calls.h \
					../include/defs.h ../utils/progressbar.h 


DDrppi: $(OBJS1) $(INCL) ../common.mk Makefile
	$(CC) $(CLINK) $(OBJS1) -o $@

%.o: %.c $(INCL) ../common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

OBJS2= wprp.o ../utils/utils.o
wprp: $(OBJS2) Makefile 
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $(OBJS2) $(CLINK)

OBJS3=countpairs_rp_pi.o ../utils/gridlink.o ../utils/utils.o ../utils/progressbar.o 
lib: libcountpairs_rp_pi.a | ../lib

libcountpairs_rp_pi.a: $(OBJS3) | ../lib
	ar rcs $@ $(OBJS3)
	cp -p $@ ../lib/
	sed -e "s/DOUBLE/$(VECTOR_TYPE)/g" countpairs_rp_pi.h > ../include/countpairs_rp_pi.h

../lib:
	mkdir -p ../lib

../bin:
	mkdir -p ../bin

install: $(target) | ../bin 
	cp -p $(target) ../bin/

.PHONY: clean clena celan install lib tests

clean:
	$(RM) $(target) $(OBJS1) $(OBJS2) $(OBJS3) libcountpairs_rp_pi.a

clena: clean

celan: celan

tests: tests_periodic tests_nonperiodic 

tests_periodic: ../tests/tests_periodic
	cd ../tests && ./test_periodic 1 3 4 5 

tests_nonperiodic: ../tests/tests_nonperiodic
	cd ../tests && ./test_nonperiodic 1 2

../tests/tests_periodic:
	$(MAKE) -C ../tests

../tests/tests_nonperiodic:
	$(MAKE) -C ../tests


