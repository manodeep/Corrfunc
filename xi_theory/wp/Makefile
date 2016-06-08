include ../../theory.options ../../common.mk

target := wp
SRC1   := wp.c countpairs_wp.c wp_driver.c ../../utils/utils.c ../../utils/gridlink.c ../../io/io.c ../../io/ftread.c 
INCL   := countpairs_wp.h wp_driver.h ../../io/ftread.h ../../io/io.h ../../utils/utils.h ../../utils/gridlink.h \
		../../utils/function_precision.h ../../utils/cellarray.h ../../utils/avx_calls.h \
		../../utils/defs.h ../../utils/sglib.h wp_kernels.c 

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

OBJS2 := countpairs_wp.o wp_driver.o ../../utils/gridlink.o ../../utils/utils.o 
ifeq ($(OUTPUT_PGBAR), 1)
  OBJS2 += ../../utils/progressbar.o
endif

libs: lib
lib: ../../lib/libcountpairs_wp.a | ../../lib ../../include

libcountpairs_wp.a: $(OBJS2) $(INCL) ../../theory.options ../../common.mk Makefile | ../../lib ../../include
	ar rcs $@ $(OBJS2)

../../lib/libcountpairs_wp.a: libcountpairs_wp.a
	cp -p $< $@
	sed -e "s/DOUBLE/$(VECTOR_TYPE)/g" countpairs_wp.h > ../../include/countpairs_wp.h

.PHONY: clean clena celan install lib tests distclean

../../lib ../../bin ../../include:
	mkdir -p $@

install: ../../bin/$(target) | ../../bin

../../bin/$(target): $(target)
	cp -p $< $@

distclean: clean | ../../lib ../../include ../../bin
	$(RM) ../../lib/libcountpairs_wp.a ../../include/countpairs_wp.h
	cd ../../bin && $(RM) $(target)

clean:
	$(RM) $(target) $(OBJS1) $(OBJS2) libcountpairs_wp.a

clena: clean

celan: celan

tests: tests_periodic 

tests_periodic: 
	$(MAKE) -C ../tests wp

