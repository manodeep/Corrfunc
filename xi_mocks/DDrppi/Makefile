include ../../mocks.options ../../common.mk

target=DDrppi wprp
all: $(target) 

$(info $$INCLUDE is [${INCLUDE}])

INCL   = countpairs_mocks.h ../../io/io.h ../../io/ftread.h ../../io/io.h  ../../utils/gridlink_mocks.h ../../utils/cellarray_mocks.h \
				../../utils/utils.h	../../utils/function_precision.h ../../utils/avx_calls.h ../../utils/defs.h ../../utils/progressbar.h

SRC1= DDrppi.c  countpairs_mocks.c ../../utils/gridlink_mocks.c \
			../../utils/utils.c ../../io/io.c ../../io/ftread.c ../../utils/progressbar.c

OBJS1=$(SRC1:.c=.o) 
DDrppi: $(OBJS1) $(INCL) ../../mocks.options ../../common.mk Makefile
	$(CC) $(OBJS1) $(CLINK) -o $@ 

SRC2=wprp.c ../../utils/utils.c
OBJS2=$(SRC2:.c=.o)
wprp: $(OBJS2) $(INCL) ../../mocks.options ../../common.mk Makefile
	$(CC) $(OBJS2) $(CLINK) -o $@ 

%.o: %.c $(INCL) ../../mocks.options ../../common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

.PHONY: clean clena celan

clean:
	$(RM) $(target) $(OBJS1) $(OBJS2) 

clena: clean
celan: clean
