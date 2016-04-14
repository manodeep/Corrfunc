include ../../mocks.options
ifneq (clean,$(findstring clean,$(MAKECMDGOALS)))
  include ../../common.mk
endif

target := vpf_mocks 
SRC1   := vpf_mocks.c countspheres_mocks.c ../../utils/gridlink_mocks.c ../../utils/utils.c \
				../../io/ftread.c ../../io/io.c ../../utils/set_cosmo_dist.c ../../utils/cosmology_params.c 
INCL   := countspheres_mocks.h ../../io/ftread.h ../../io/io.h ../../utils/utils.h ../../utils/gridlink_mocks.h \
				../../utils/function_precision.h ../../utils/cellarray.h ../../utils/avx_calls.h \
				../../utils/defs.h ../../utils/set_cosmo_dist.h ../../utils/cosmology_params.h

ifeq ($(OUTPUT_PGBAR), 1)
  SRC1 += ../../utils/progressbar.c
  INCL += ../../utils/progressbar.h 
endif

OBJS1  := $(SRC1:.c=.o)

all: $(target) $(SRC1) $(INCL) Makefile ../../mocks.options ../../common.mk

$(target): $(OBJS1) $(INCL) ../../mocks.options ../../common.mk Makefile 
	$(CC) $(OBJS1) $(CLINK) $(GSL_LINK) -o $@

%.o: %.c $(INCL) ../../mocks.options ../../common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

countspheres_mocks.o: countspheres_mocks.c $(INCL) ../../common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) $(GSL_CFLAGS) $(INCLUDE) -c $< -o $@

../../utils/set_cosmo_dist.o:../../utils/set_cosmo_dist.c ../../mocks.options ../../common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) $(GSL_CFLAGS) $(INCLUDE) -c $< -o $@

libs: lib
lib:  libcountspheres_mocks.a | ../../lib ../../include

SRC2   = countspheres_mocks.c ../../utils/gridlink_mocks.c ../../utils/utils.c ../../utils/progressbar.c \
				../../utils/set_cosmo_dist.c ../../utils/cosmology_params.c

OBJS2  = $(SRC2:.c=.o)
libcountspheres_mocks.a: $(OBJS2) ../../mocks.options ../../common.mk Makefile | ../../lib ../../include
	ar rcs $@ $(OBJS2)
	cp -p $@ ../../lib/
	sed -e "s/DOUBLE/$(VECTOR_TYPE)/g" countspheres_mocks.h > ../../include/countspheres_mocks.h 

.PHONY: clean clena celan install lib tests distclean

../../lib ../../bin ../../include:
	mkdir -p $@

install: ../../bin/$(target) | ../../bin

../../bin/$(target): $(target)
	cp -p $< $@

distclean: clean | ../../include ../../bin ../../lib
	$(RM) ../../bin/$(target)
	$(RM) ../../include/countspheres_mocks.h
	$(RM) ../../lib/libcountspheres_mocks.a

clean:
	$(RM) $(target) $(OBJS1) $(OBJS2) libcountspheres_mocks.a

clena: clean
celan: clean


tests: 
	$(MAKE) -C ../tests	vpf



