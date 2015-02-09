include ../../mocks.options ../../common.mk

target = tests_nonperiodic

INCL  = ../../utils/gridlink.c ../../io/io.c ../../io/ftread.c \
				../xi_of_r/countpairs.c ../xi_rp_pi/countpairs_rp_pi.c ../wp/countpairs_wp.c \
				../vpf/countspheres.c 

SRC1   = tests_nonperiodic.c ../../utils/utils.c ../../utils/progressbar.c 
OBJS1  = $(SRC1:.c=.o)

all: $(target) $(INCL) $(SRC1) 

tests_nonperiodic: $(OBJS1) $(INCL) ../mocks.options ../../common.mk Makefile 
	$(CC) $(OBJS1) $(CLINK) $(GSL_LINK) -lgsl -lgslcblas -o $@

%.o: %.c $(INCL) ../mocks.options ../../common.mk Makefile
	$(CC) $(OPT) $(GSL_CFLAGS) $(CFLAGS) $(INCLUDE) -c $< -o $@

tests: tests_nonperiodic
	./test_nonperiodic

DDrppi: tests_nonperiodic
	./test_nonperiodic 1 2  

wtheta: tests_nonperiodic
	./test_nonperiodic 3 4 

vpf: tests_nonperiodic
	./test_nonperiodic 5

%.o: %.c ../../mocks.options ../common.mk Makefile $(INCL)
	$(CC) $(OPT) $(CFLAGS) $(GSL_CFLAGS) $(INCLUDE) -c $< -o $@

.PHONY: clean clena celan install tests 

clean:
	$(RM) $(target) $(OBJS1) $(OBJS2)

clena: clean

celan: clean


