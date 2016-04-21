include ../../theory.options
ifneq (clean,$(findstring clean,$(MAKECMDGOALS)))
  include ../../common.mk
endif

target := xi 
SRC1   := xi.c countpairs_xi.c xi_driver.c ../../utils/gridlink.c ../../utils/utils.c ../../io/ftread.c ../../io/io.c 
INCL   := countpairs_xi.h ../../io/ftread.h ../../io/io.h ../../utils/utils.h ../../utils/gridlink.h \
		../../utils/function_precision.h ../../utils/cellarray.h ../../utils/avx_calls.h \
		../../utils/defs.h xi_driver.h ../../utils/sglib.h xi_kernels.c 

ifeq ($(OUTPUT_PGBAR), 1)
SRC1 += ../../utils/progressbar.c
INCL += ../../utils/progressbar.h 
endif
OBJS1  := $(SRC1:.c=.o)

all: $(target) $(SRC1) $(INCL) ../../theory.options ../../common.mk Makefile 

$(target): $(OBJS1) $(INCL) ../../common.mk ../../theory.options Makefile 
	$(CC) $(OBJS1) $(CLINK) -o $@

%.o: %.c $(INCL) ../../theory.options ../../common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

SRC2  := countpairs_xi.c xi_driver.c ../../utils/gridlink.c ../../utils/utils.c
ifeq ($(OUTPUT_PGBAR), 1)
SRC2 += ../../utils/progressbar.c
endif
OBJS2 := $(SRC2:.c=.o)

libs: lib
lib:  ../../lib/libcountpairs_xi.a | ../../lib ../../include 

libcountpairs_xi.a: $(OBJS2) $(INCL) ../../common.mk ../../theory.options Makefile | ../../lib ../../include
	ar rcs $@ $(OBJS2)

../../lib/libcountpairs_xi.a: libcountpairs_xi.a
	cp -p $< $@
	sed -e "s/DOUBLE/$(VECTOR_TYPE)/g" countpairs_xi.h > ../../include/countpairs_xi.h 

.PHONY: clean clena celan install lib tests distclean

../../lib ../../bin ../../include:
	mkdir -p $@

install: ../../bin/$(target) | ../../bin

../../bin/$(target): $(target)
	cp -p $< $@

distclean: clean | ../../lib ../../include ../../bin
	$(RM) ../../lib/libcountpairs_xi.a ../../include/countpairs_xi.h
	cd ../../bin && $(RM) $(target)

clean:
	$(RM) $(target) $(OBJS1) $(OBJS2) libcountpairs_xi.a

clena: clean
celan: clean

tests: tests_periodic 

tests_periodic: 
	$(MAKE) -C ../tests xi

