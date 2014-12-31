include ../common.mk

target=DD 
SRC1   = DD.c countpairs.c ../utils/gridlink.c ../utils/utils.c ../io/ftread.c ../io/io.c 
OBJS1  = $(SRC1:.c=.o)
INCL   = countpairs.h ../io/ftread.h ../io/io.h ../utils/utils.h ../utils/gridlink.h \
					../include/function_precision.h ../include/cellarray.h ../include/avx_calls.h \
					../include/defs.h 

all: $(target) $(SRC1) $(INCL)

$(target): $(OBJS1) $(INCL) ../common.mk Makefile 
	$(CC) $(CFLAGS) $(OPT) $(CLINK) $(OBJS1) -o $@

.c.o: $(INCL) ../common.mk Makefile
	$(CC) $(CFLAGS) $(OPT) $(INCLUDE) -c $< -o $@

SRC2   = countpairs.c ../utils/gridlink.c ../utils/utils.c
OBJS2  = $(SRC2:.c=.o)
lib: $(OBJS2)
	mkdir -p ../lib
# 	$(CC) $(CLINK) $(LDFLAGS) $(OBJS2) -shared -Wl,-soname,libcountpairs.so.$(MAJOR) -o libcountpairs.so.$(MAJOR).0.$(MINOR) 
# 	ldconfig -n ./ 
# 	ln -s libcountpairs.so.$(MAJOR) libcountpairs.so
# 	cp -p libcountpairs* ../lib
# 	cp -p countpairs.h ../include

	$(CC) $(CFLAGS) $(OPT) $(CLINK) $(LDFLAGS) $(OBJS2) -shared -o libcountpairs.so.$(MAJOR).0.$(MINOR) 
	rm -f libcountpairs.so 
	ln -s libcountpairs.so.$(MAJOR).0.$(MINOR) libcountpairs.so
	cp -p libcountpairs.so.$(MAJOR).0.$(MINOR) libcountpairs.so ../lib
	cp -p countpairs.h ../include


.PHONY: clean clena celan install lib tests 

clean:
	rm -f $(target) $(OBJS1) $(OBJS2) libcountpairs.so 

clena:
	rm -f $(target) $(OBJS1) $(OBJS2) libcountpairs.so 

celan:
	rm -f $(target) $(OBJS1) $(OBJS2) libcountpairs.so 

install:
	mkdir -p ../bin
	cp -p $(target) ../bin/


tests: tests_periodic tests_nonperiodic tests_avx tests_omp
	$(MAKE) clean
	$(MAKE)

tests_periodic: 
	$(info $$OPT is [$(OPT)])
	$(MAKE) clean
	$(MAKE) 

	./DD ../tests/data/gals_Mr19.ff f ../tests/data/gals_Mr19.ff f ../tests/bins 4 > new_DD
	if cmp -s ../tests/Mr19_DD_correct new_DD  ; then
		echo "Passed DD test"
		rm -f "new_DD"
	else
		echo "Failed DD test"
		exit
	endif

tests_nonperiodic:
	$(info $$OPT is [$(OPT)])
	$(MAKE) clean
	$(MAKE) 

tests_avx:
	$(info $$OPT is [$(OPT)])
	$(MAKE) clean
	$(MAKE) 

tests_omp:
	$(info $$OPT is [$(OPT)])
	$(MAKE) clean
	$(MAKE) 

