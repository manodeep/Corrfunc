include ../../theory.options ../../common.mk

target = test_periodic test_nonperiodic
SRC1   = test_periodic.c ../../utils/utils.c ../../utils/progressbar.c 
OBJS1  = $(SRC1:.c=.o)
INCL  = ../../utils/gridlink.c ../../io/io.c ../../io/ftread.c \
				../xi_of_r/countpairs.c ../xi_rp_pi/countpairs_rp_pi.c ../wp/countpairs_wp.c \
				../vpf/countspheres.c ../xi/countpairs_xi.c


SRC2   = test_nonperiodic.c ../../utils/utils.c ../../utils/progressbar.c 
OBJS2  = $(SRC2:.c=.o)

all: $(target) $(INCL)

test_periodic: $(OBJS1) ../../theory.options ../../common.mk Makefile \
	../xi_of_r/countpairs.c ../xi_rp_pi/countpairs_rp_pi.c ../wp/countpairs_wp.c
	$(CC) $(OBJS1) $(CLINK) $(GSL_LINK) -lgsl -lgslcblas -o $@

test_nonperiodic: $(OBJS2) ../../theory.options ../../common.mk Makefile \
	../xi_of_r/countpairs.c ../xi_rp_pi/countpairs_rp_pi.c 
	$(CC) $(CLINK) $(OBJS2) -o $@

%.o: %.c $(INCL) ../../theory.options ../../common.mk Makefile
	$(CC) $(OPT) $(GSL_CFLAGS) $(CFLAGS) $(INCLUDE) -c $< -o $@

tests: $(target)
	./test_periodic
	./test_nonperiodic

wp: test_periodic
	./test_periodic 2

xi_rp_pi: test_periodic test_nonperiodic
	./test_periodic 1 5 6 7
	./test_nonperiodic 1 2 

xi_of_r: test_periodic test_nonperiodic
	./test_periodic 0
	./test_nonperiodic 0

vpf: test_periodic
	./test_periodic 3

xi: test_periodic
	./test_periodic 4

.PHONY: clean clena celan install lib tests 

clean:
	$(RM) $(target) $(OBJS1) $(OBJS2)

clena: clean

celan: clean


