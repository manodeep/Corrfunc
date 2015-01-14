include ../common.mk

target=DDrppi wprp
all: $(target) 


SRC1= DDrppi.c countpairs_rp_pi.c ../utils/utils.c ../utils/gridlink.c ../io/io.c ../io/ftread.c 
OBJS1  = $(SRC1:.c=.o)
INCL   = countpairs_rp_pi.h ../io/ftread.h ../io/io.h ../utils/utils.h ../utils/gridlink.h \
					../include/function_precision.h ../include/cellarray.h ../include/avx_calls.h \
					../include/defs.h 


DDrppi: $(OBJS1) $(INCL) ../common.mk Makefile
	$(CC) $(CLINK) $(OBJS1) -o $@

.c.o: $(INCL) ../common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

OBJS2= wprp.o ../utils/utils.o
wprp: $(OBJS2) Makefile 
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $(OBJS2) $(CLINK)

OBJS3=countpairs_rp_pi.o ../utils/gridlink.o ../utils/utils.o 
lib: libcountpairs_rp_pi.a

libcountpairs_rp_pi.a: ../lib $(OBJS3)
	ar rcs $@ $(OBJS3)
	cp -p $@ ../lib/
	cp -p countpairs_rp_pi.h ../include 

../lib:
	mkdir -p ../lib

.PHONY: clean clena celan 

clean:
	rm -f $(target) $(OBJS1) $(OBJS2) $(OBJS3) libcountpairs_rp_pi.a

clena: clean

celan: celan

install: ../bin
	cp -p $(target) ../bin/

../bin:
	mkdir -p ../bin
