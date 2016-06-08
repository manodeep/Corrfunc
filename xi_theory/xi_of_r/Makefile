include ../../theory.options ../../common.mk
#OPT+=-DSILENT

target := DD 
SRC1   := DD.c countpairs.c countpairs_driver.c ../../utils/gridlink.c ../../utils/utils.c ../../io/ftread.c ../../io/io.c 
INCL   := countpairs.h ../../io/ftread.h ../../io/io.h ../../utils/utils.h ../../utils/gridlink.h \
		../../utils/function_precision.h ../../utils/cellarray.h ../../utils/avx_calls.h \
		../../utils/defs.h countpairs_kernels.c

ifeq ($(OUTPUT_PGBAR), 1)
SRC1 += ../../utils/progressbar.c
INCL += ../../utils/progressbar.h 
endif
OBJS1  := $(SRC1:.c=.o)

all: $(target) $(SRC1) $(INCL) ../../theory.options ../../common.mk Makefile 

$(target): $(OBJS1) $(INCL) ../../theory.options ../../common.mk Makefile 
	$(CC) $(OBJS1) $(CLINK) -o $@

%.o: %.c $(INCL) ../../theory.options ../../common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

SRC2   = countpairs.c countpairs_driver.c ../../utils/gridlink.c ../../utils/utils.c
ifeq ($(OUTPUT_PGBAR), 1)
SRC2 += ../../utils/progressbar.c
endif

OBJS2  = $(SRC2:.c=.o)
libs: lib
lib:  ../../lib/libcountpairs.a | ../../lib

libcountpairs.a: $(OBJS2) $(INCL) ../../theory.options ../../common.mk Makefile | ../../lib ../../include
	ar rcs $@ $(OBJS2)

../../lib/libcountpairs.a: libcountpairs.a
	cp -p $< $@
	sed -e "s/DOUBLE/$(VECTOR_TYPE)/g" countpairs.h > ../../include/countpairs.h

.PHONY: clean clena celan install lib tests distclean

../../lib ../../bin ../../include:
	mkdir -p $@

install: ../../bin/$(target) | ../../bin

../../bin/$(target): $(target)
	cp -p $< $@

distclean: clean | ../../lib ../../include ../../bin
	$(RM) ../../lib/libcountpairs.a ../../include/countpairs.h
	cd ../../bin && $(RM) $(target)

clean:
	$(RM) $(target) $(OBJS1) $(OBJS2) libcountpairs.a

clena: clean
celan: clean


tests: 
	$(MAKE) -C ../tests xi_of_r


