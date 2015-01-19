include ../common.mk

target = run_correlations

SRC1   = run_correlations.c ../io/ftread.c ../io/io.c  ../utils/utils.c
OBJS1  = $(SRC1:.c=.o)
INCL   = ../io/ftread.h ../io/io.h ../utils/utils.h \
				../include/defs.h ../include/countpairs.h ../include/countpairs_rp_pi.h ../include/countpairs_wp.h 
C_LIBRARIES = ../lib/libcountpairs.a ../lib/libcountpairs_rp_pi.a ../lib/libcountpairs_wp.a 

all: $(target) $(SRC1) $(INCL)

$(target): $(OBJS1) $(INCL) $(C_LIBRARIES) ../common.mk Makefile 
	$(CC) $(OBJS1) $(C_LIBRARIES) $(CLINK) -o $@

%.o: %.c $(INCL) ../common.mk Makefile
	$(CC) $(CFLAGS) $(OPT) $(INCLUDE) -c $< -o $@

.PHONY: clean clena celan 

clean:
	$(RM) $(target) $(OBJS1) 

clena: clean

celan: clean



