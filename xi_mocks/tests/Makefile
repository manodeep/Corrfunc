include ../../mocks.options ../../common.mk

target = tests_mocks

SRC1   = tests_mocks.c ../../utils/utils.c ../../utils/progressbar.c ../../utils/set_cosmo_dist.c ../../utils/cosmology_params.c 
OBJS1  = $(SRC1:.c=.o)
INCL  = ../../utils/gridlink_mocks.c ../../io/io.c ../../io/ftread.c \
				../DDrppi/countpairs_rp_pi_mocks.c ../wtheta/countpairs_theta_mocks.c \
				../vpf/countspheres_mocks.c ../../utils/set_cosmo_dist.h ../../utils/cosmology_params.h 

all: $(target) $(INCL) 

tests_mocks: $(OBJS1) $(INCL) ../../mocks.options ../../common.mk Makefile 
	$(CC) $(OBJS1) $(CLINK) $(GSL_LINK) -lgsl -lgslcblas -o $@

%.o: %.c $(INCL) ../../mocks.options ../../common.mk Makefile
	$(CC) $(OPT) $(GSL_CFLAGS) $(CFLAGS) $(INCLUDE) -c $< -o $@

tests: tests_mocks
	./tests_mocks

DDrppi: tests_mocks
	./tests_mocks 0 3  

wtheta: tests_mocks
	./tests_mocks 1 4 

vpf: tests_mocks
	./tests_mocks 2 5

.PHONY: clean clena celan install tests 

clean:
	$(RM) $(target) $(OBJS1)

clena: clean

celan: clean


