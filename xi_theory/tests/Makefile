include ../../theory.options
ifneq (clean,$(findstring clean,$(MAKECMDGOALS)))
  include ../../common.mk
endif

target := test_periodic test_nonperiodic
SRC1   := test_periodic.c ../../utils/utils.c ../../utils/progressbar.c 
OBJS1  := $(SRC1:.c=.o)
INCL   := ../../utils/gridlink.c ../../io/io.c ../../io/ftread.c \
				../xi_of_r/countpairs.c ../xi_rp_pi/countpairs_rp_pi.c ../wp/countpairs_wp.c \
				../vpf/countspheres.c ../xi/countpairs_xi.c

SRC2   := test_nonperiodic.c ../../utils/utils.c ../../utils/progressbar.c 
OBJS2  := $(SRC2:.c=.o)

all: $(target) $(INCL) uncompress

test_periodic: $(OBJS1) $(INCL) ../../theory.options ../../common.mk Makefile 
	$(CC) $(OBJS1) $(GSL_LINK) -lgsl -lgslcblas $(CLINK) -o $@

test_nonperiodic: $(OBJS2) $(INCL) ../../theory.options ../../common.mk Makefile 
	$(CC) $(OBJS2) $(CLINK) -o $@

%.o: %.c $(INCL) ../../theory.options ../../common.mk Makefile
	$(CC) $(OPT) $(GSL_CFLAGS) $(CFLAGS) $(INCLUDE) -c $< -o $@

tests: $(target) 
	./test_nonperiodic
	./test_periodic


uncompress: | data
	@{\
		cd data ; \
		for f in $$(ls -f *.gz 2>/dev/null) ; do \
			gunzip $$f ; \
		done ;\
		cd ..; \
	}
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

.PHONY: clean clena celan install tests uncompress

clean:
	$(RM) $(target) $(OBJS1) $(OBJS2)

clena: clean

celan: clean


