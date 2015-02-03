include ../../common.mk

target=DDtheta wtheta
all: $(target) 

INCL   = countpairs_theta.h ../../io/io.h ../../io/ftread.h ../../io/io.h ../../utils/utils.h ../../utils/gridlink_mocks.h \
        ../../utils/function_precision.h ../../utils/cellarray.h ../../utils/avx_calls.h \
        ../../utils/defs.h ../../utils/progressbar.h

SRC1= DDtheta.c  countpairs_theta.c ../../utils/gridlink_mocks.c \
	../../utils/utils.c ../../io/io.c ../../io/ftread.c ../../utils/progressbar.c  
OBJS1=$(SRC1:.c=.o) 
DDtheta: $(OBJS1) $(INCL) ../../common.mk Makefile
	$(CC) $(OBJS1) $(CLINK) -o $@ 

SRC2=wtheta.c ../../utils/utils.c
OBJS2=$(SRC2:.c=.o)
wtheta: $(OBJS2) $(INCL) ../../common.mk Makefile
	$(CC) $(OBJS2) $(CLINK) -o $@ 

%.o: %.c $(INCL) ../../common.mk Makefile
	$(CC) $(DATA_OPT) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

.PHONY: clean clena celan

clean:
	rm -f $(target) $(OBJS1) $(OBJS2) 

clena: clean
celan: clean
