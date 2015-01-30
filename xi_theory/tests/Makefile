include ../../common.mk

target = test_periodic test_nonperiodic
SRC1   = test_periodic.c ../../utils/utils.c ../../utils/progressbar.c 
OBJS1  = $(SRC1:.c=.o)
INCL  = ../../utils/gridlink.c ../../io/io.c ../../io/ftread.c \
				../xi_of_r/countpairs.c ../xi_rp_pi/countpairs_rp_pi.c ../wp/countpairs_wp.c 


SRC2   = test_nonperiodic.c ../../utils/utils.c ../../utils/progressbar.c 
OBJS2  = $(SRC2:.c=.o)

# $(info $$INCLUDE is [${INCLUDE}])
# $(info $$INCL is [${INCL}])


all: $(target) $(INCL)

test_periodic: $(OBJS1) ../../common.mk Makefile \
	../xi_of_r/countpairs.c ../xi_rp_pi/countpairs_rp_pi.c ../wp/countpairs_wp.c
	$(CC) $(CLINK) $(OBJS1) -o $@

test_nonperiodic: $(OBJS2) ../../common.mk Makefile \
	../xi_of_r/countpairs.c ../xi_rp_pi/countpairs_rp_pi.c 
	$(CC) $(CLINK) $(OBJS2) -o $@

%.o: %.c $(INCL) ../../common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

tests: $(target)
	./test_periodic
	./test_nonperiodic

wp: test_periodic
	./test_periodic 2

xi_rp_pi: test_periodic test_nonperiodic
	./test_periodic 1 3 4 5
	./test_nonperiodic 1 2 

xi_of_r: test_periodic test_nonperiodic
	./test_periodic 0
	./test_nonperiodic 0

%.o: %.c ../common.mk Makefile $(INCL)
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

.PHONY: clean clena celan install lib tests 

clean:
	$(RM) $(target) $(OBJS1) $(OBJS2)

clena: clean

celan: clean


