include ../common.mk

target=DD 
SRC1   = DD.c countpairs.c ../utils/gridlink.c ../utils/utils.c ../io/ftread.c ../io/io.c 
OBJS1  = $(SRC1:.c=.o)
INCL   = countpairs.h ../io/ftread.h ../io/io.h ../utils/utils.h ../utils/gridlink.h \
					../include/function_precision.h ../include/cellarray.h ../include/avx_calls.h \
					../include/defs.h 

all: $(target) $(SRC1) $(INCL)

$(target): $(OBJS1) $(INCL) ../common.mk Makefile 
	$(CC) $(CLINK) $(OBJS1) -o $@

%.o: %.c $(INCL) ../common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

SRC2   = countpairs.c ../utils/gridlink.c ../utils/utils.c
OBJS2  = $(SRC2:.c=.o)
lib:  libcountpairs.a

libcountpairs.a: $(OBJS2)
	ar rcs $@ $(OBJS2)
	cp -p $@ ../lib/
	cp -p countpairs.h ../include


.PHONY: clean clena celan install lib tests 

clean:
	$(RM) $(target) $(OBJS1) $(OBJS2) libcountpairs.a

clena: clean


celan: clean

install: 
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
		$(RM) "new_DD"
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

