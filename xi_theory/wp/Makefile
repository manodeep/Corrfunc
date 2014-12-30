include ../common.mk

target=wp
SRC1= wp.c countpairs_wp.c ../utils/utils.c ../utils/gridlink.c ../io/io.c ../io/ftread.c 
OBJS1  = $(SRC1:.c=.o)
INCL   = countpairs_wp.h ../io/ftread.h ../io/io.h ../utils/utils.h ../utils/gridlink.h \
					../include/function_precision.h ../include/cellarray.h ../include/avx_calls.h \
					../include/defs.h ../include/sglib.h 

all: $(target) $(SRC1) $(INCL)

$(target): $(OBJS1) $(INCL) ../common.mk Makefile 
	$(CC) $(CFLAGS) $(OPT) $(CLINK) $(OBJS1) -o $@

.c.o: $(INCL) ../common.mk Makefile
	$(CC) $(CFLAGS) $(OPT) $(INCLUDE) -c $< -o $@

OBJS2=countpairs_wp.o ../utils/gridlink.o ../utils/utils.o 
lib: $(OBJS2)
	mkdir -p ../lib
	$(CC) $(LDFLAGS) $(OBJS2) -shared -Wl,-soname,libcountpairs_wp.so.$(MAJOR) -o libcountpairs_wp.so.$(MAJOR).0.$(MINOR) 
	ldconfig -n ./ 
	ln -s libcountpairs_wp.so.$(MAJOR) libcountpairs_wp.so
	mv libcountpairs_wp* ../lib
	cp -p countpairs_wp.h ../include 

.PHONY: clean clena celan install lib

clean:
	rm -f $(target) $(OBJS1) $(OBJS2) libcountpairs_wp.so 

clena:
	rm -f $(target) $(OBJS1) $(OBJS2) libcountpairs_wp.so 

celan:
	rm -f $(target) $(OBJS1) $(OBJS2) libcountpairs_wp.so

install:
	mkdir -p ../bin
	cp -p $(target) ../bin/


