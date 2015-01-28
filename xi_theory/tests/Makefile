include ../common.mk

target = test_periodic test_nonperiodic
SRC1   = test_periodic.c ../utils/utils.c ../utils/progressbar.c 
OBJS1  = $(SRC1:.c=.o)

SRC2   = test_nonperiodic.c ../utils/utils.c ../utils/progressbar.c 
OBJS2  = $(SRC2:.c=.o)

all: $(target) $(SRC1) $(SRC2)

test_periodic: $(OBJS1) ../common.mk Makefile \
	../xi_of_r/countpairs.c ../xi_rp_pi/countpairs_rp_pi.c ../wp/countpairs_wp.c
	$(CC) $(CLINK) $(OBJS1) -o $@

test_nonperiodic: $(OBJS2) ../common.mk Makefile \
	../xi_of_r/countpairs.c ../xi_rp_pi/countpairs_rp_pi.c 
	$(CC) $(CLINK) $(OBJS2) -o $@

%.o: %.c $(INCL) ../common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@


.PHONY: clean clena celan install lib tests 

clean:
	$(RM) $(target) $(OBJS1) 

clena: clean

celan: clean


