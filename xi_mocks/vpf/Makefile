include ../../common.mk

target=DDrppi_data
all: $(target) 

dirs: 
	mkdir -p ../../bin ../../lib ../../include

INCL   = countpairs_data.h ../../io/ftread.h ../../io/io.h ../../utils/utils.h ../../utils/gridlink_data.h \
        ../../utils/function_precision.h ../../utils/cellarray.h ../../utils/avx_calls.h \
        ../../utils/defs.h ../../utils/progressbar.h

SRC1=DDrppi.c  countpairs_data.c ../../utils/gridlink_data.c \
	../../utils.c ../../io/ftread.c ../../io/io.c ../../utils/progressbar.c \
	../../utils/set_cosmo_dist.c ../../utils/cosmology_params.c
OBJS1=$(SRC1:.c=.o)

$(target): $(OBJS1) $(INCL) Makefile ../../common.mk
	$(CC) $(OBJS1) $(CLINK) $(GSL_LINK) -o $@ 

%.o: %.c $(INCL) ../../common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

set_cosmology.o: set_cosmology.c ../../common.mk Makefile
	$(CC) $(OPT) $(GSL_CFLAGS) $(CFLAGS) $(INCLUDE) -c $< -o $@

.PHONY: clean clena distclean install

install: ../../bin/$(target)  | dirs
	cp -p $(target) ../../bin/

distclean: clean | dirs
	cd ../../bin && $(RM) $(target)

clean:
	rm -f $(target) $(OBJS1) 

clena: clean

