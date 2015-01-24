include ../common.mk

target = run_correlations

SRC1   := run_correlations.c ../io/ftread.c ../io/io.c  ../utils/utils.c ../utils/progressbar.c 
OBJS1  := $(SRC1:.c=.o)
C_LIBRARIES = ../lib/libcountpairs.a ../lib/libcountpairs_rp_pi.a ../lib/libcountpairs_wp.a 

all: $(target) 

$(target): Makefile ../common.mk $(C_LIBRARIES) $(OBJS1) 
	$(CC) $(OBJS1) $(C_LIBRARIES) $(CLINK) -o $@

../lib/libcountpairs.a: 
	$(MAKE) -C ../xi_of_r/ lib

../lib/libcountpairs_rp_pi.a: 
	$(MAKE) -C ../xi_rp_pi/ lib

../lib/libcountpairs_wp.a: 
	$(MAKE) -C ../wp/ lib

%.o: %.c ../common.mk Makefile
	$(CC) $(CFLAGS) $(OPT) $(INCLUDE) -c $< -o $@

.PHONY: clean clena celan 

clean:
	$(RM) $(target) $(OBJS1) 

clena: clean

celan: clean



