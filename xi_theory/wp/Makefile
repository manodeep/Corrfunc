include ../../theory.options ../../common.mk

target := wp
SRC1   := wp.c countpairs_wp.c countpairs_wp_impl_double.c countpairs_wp_impl_float.c \
          ../../io/io.c ../../io/ftread.c \
          ../../utils/utils.c ../../utils/gridlink_impl_float.c ../../utils/gridlink_impl_double.c \
          ../../utils/progressbar.c ../../utils/cpu_features.c
INCL   := wp_kernels_double.c wp_kernels_float.c countpairs_wp.h countpairs_wp_impl_float.h countpairs_wp_impl_double.h \
          ../../io/ftread.h ../../io/io.h ../../utils/utils.h \
          ../../utils/gridlink_impl_float.h ../../utils/gridlink_impl_double.h \
          ../../utils/function_precision.h ../../utils/cellarray_double.h ../../utils/cellarray_float.h ../../utils/avx_calls.h \
          ../../utils/defs.h ../../utils/sglib.h ../../utils/progressbar.h ../../utils/cpu_features.h
OBJS1  := $(SRC1:.c=.o)

all: $(target) $(SRC1) $(INCL) ../../theory.options ../../common.mk Makefile 

$(target): $(OBJS1) $(INCL) ../../theory.options ../../common.mk Makefile 
	$(CC) $(OBJS1) $(CLINK) -o $@

%_float.c: %.c.src Makefile
	@$(ECHO_COMMAND) '/* This file is auto-generated from $< */' > $@
	sed -e "s/DOUBLE/float/g" $< >> $@
%_float.h: %.h.src Makefile
	@$(ECHO_COMMAND) '/* This file is auto-generated from $< */' > $@
	sed -e "s/DOUBLE/float/g"  $< >> $@

%_double.c: %.c.src Makefile
	@$(ECHO_COMMAND) "/* This file is auto-generated from $< */" > $@
	sed -e "s/DOUBLE/double/g"  $< >> $@
%_double.h: %.h.src Makefile
	@$(ECHO_COMMAND) '/* This file is auto-generated from $< */' > $@
	sed -e "s/DOUBLE/double/g"  $< >> $@

%_double.o: %_double.c $(INCL) ../../theory.options ../../common.mk Makefile
	$(CC) -DDOUBLE_PREC $(CFLAGS) $(INCLUDE) -c $< -o $@

%_float.o: %_float.c %_float.h $(INCL)
	$(CC) -DNDOUBLE_PREC $(CFLAGS) $(INCLUDE) -c $< -o $@

countpairs_wp.o: countpairs_wp.c countpairs_wp.h countpairs_wp_impl_double.h countpairs_wp_impl_float.h ../../theory.options ../../common.mk Makefile
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

$(target).o: $(target).c ../../theory.options ../../common.mk Makefile $(INCL)
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

%.o: %.c $(INCL) ../../theory.options ../../common.mk Makefile
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@


OBJS2 := countpairs_wp.o countpairs_wp_impl_double.o countpairs_wp_impl_float.o \
         ../../utils/gridlink_impl_double.o ../../utils/gridlink_impl_float.o \
         ../../utils/utils.o ../../utils/progressbar.o ../../utils/cpu_features.o
libs: lib
lib: ../../lib/libcountpairs_wp.a | ../../lib ../../include

libcountpairs_wp.a: $(OBJS2) $(INCL) ../../theory.options ../../common.mk Makefile | ../../lib ../../include
	ar rcs $@ $(OBJS2)

../../lib/libcountpairs_wp.a: libcountpairs_wp.a
	cp -p $< $@
	cp -p countpairs_wp.h > ../../include/countpairs_wp.h

.PHONY: clean clena celan install lib tests distclean all

../../lib ../../bin ../../include:
	mkdir -p $@

install: ../../bin/$(target) | ../../bin

../../bin/$(target): $(target)
	cp -p $< $@

distclean: clean | ../../lib ../../include ../../bin
	$(RM) ../../lib/libcountpairs_wp.a ../../include/countpairs_wp.h
	cd ../../bin && $(RM) $(target)

clean:
	$(RM) $(target) $(OBJS1) $(OBJS2) libcountpairs_wp.a *_float.[ch] *_double.[ch]

clena: clean

celan: celan

tests: tests_periodic 

tests_periodic: 
	$(MAKE) -C ../tests wp

