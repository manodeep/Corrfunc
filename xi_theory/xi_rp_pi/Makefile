include ../common.mk

target=DDrppi wprp
all: $(target) 


SRC1= DDrppi.c countpairs_rp_pi.c ../utils/utils.c ../utils/gridlink.c ../io/io.c ../io/ftread.c 
OBJS1  = $(SRC1:.c=.o)
INCL   = countpairs_rp_pi.h ../io/ftread.h ../io/io.h ../utils/utils.h ../utils/gridlink.h \
					../include/function_precision.h ../include/cellarray.h ../include/avx_calls.h \
					../include/defs.h 


DDrppi: $(OBJS1) $(INCL) ../common.mk Makefile
	$(CC) $(CFLAGS) $(OPT) $(CLINK) $(OBJS1) -o $@

.c.o: $(INCL) ../common.mk Makefile
	$(CC) $(CFLAGS) $(OPT) $(INCLUDE) -c $< -o $@

OBJS2= wprp.o ../utils/utils.o
wprp: $(OBJS2) Makefile 
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $(OBJS2) $(CLINK)

OBJS3=countpairs_rp_pi.o ../utils/gridlink.o ../utils/utils.o 
lib: $(OBJS3)
	mkdir -p ../lib
# 	$(CC) $(LDFLAGS) $(OBJS3) -shared -Wl,-soname,libcountpairs_rp_pi.so.$(MAJOR) -o libcountpairs_rp_pi.so.$(MAJOR).0.$(MINOR) 
# 	ldconfig -n ./ 
# 	ln -s libcountpairs_rp_pi.so.$(MAJOR) libcountpairs_rp_pi.so
# 	cp -p libcountpairs_rp_pi* ../lib
# 	cp -p countpairs_rp_pi.h ../include 

	$(CC) $(CFLAGS) $(CLINK) $(OPT) $(LDFLAGS) $(OBJS3) -shared -o libcountpairs_rp_pi.so.$(MAJOR).0.$(MINOR) 
	cp -p libcountpairs_rp_pi.so.$(MAJOR).0.$(MINOR) ../lib
	rm -f ../lib/libcountpairs_rp_pi.so 
	ln -s ../lib/libcountpairs_rp_pi.so.$(MAJOR).0.$(MINOR) ../lib/libcountpairs_rp_pi.so
	cp -p countpairs_rp_pi.h ../include 

.PHONY: clean clena celan 

clean:
	rm -f $(target) $(OBJS1) $(OBJS2) $(OBJS3) libcountpairs_rp_pi.so 

clena:
	rm -f $(target) $(OBJS1) $(OBJS2) $(OBJS3) libcountpairs_rp_pi.so 

celan:
	rm -f $(target) $(OBJS1) $(OBJS2) $(OBJS3) libcountpairs_rp_pi.so

install:
	mkdir -p ../bin
	cp -p $(target) ../bin/


