include ../../theory.options
ifneq (clean,$(findstring clean,$(MAKECMDGOALS)))
  include ../../common.mk
endif

targets := test_periodic test_nonperiodic
ifneq ($(COMPILE_PYTHON_EXT), 0)
targets += python_lib
else
$(warning $(ccmagenta) Skipping python tests since python or numpy is unavailable $(ccreset))
endif

SRC1   := test_periodic.c ../../utils/utils.c ../../utils/progressbar.c 
OBJS1  := $(SRC1:.c=.o)
INCL   := ../../utils/gridlink.c ../../io/io.c ../../io/ftread.c \
	../xi_of_r/countpairs.c ../xi_of_r/countpairs_driver.c ../xi_of_r/countpairs_kernels.c \
        ../xi_rp_pi/countpairs_rp_pi.c ../xi_rp_pi/countpairs_rp_pi_driver.c ../xi_rp_pi/countpairs_rp_pi_kernels.c \
        ../wp/countpairs_wp.c ../wp/wp_driver.c ../wp/wp_kernels.c \
	../vpf/countspheres.c \
        ../xi/countpairs_xi.c ../xi/xi_driver.c ../xi/xi_kernels.c

SRC2   := test_nonperiodic.c ../../utils/utils.c ../../utils/progressbar.c 
OBJS2  := $(SRC2:.c=.o)

all: tests $(targets) $(INCL) uncompress ../../theory.options ../../common.mk Makefile

test_periodic: $(OBJS1) $(INCL) ../../theory.options ../../common.mk Makefile 
	$(CC) $(OBJS1) $(GSL_LINK) -lgsl -lgslcblas $(CLINK) -o $@

test_nonperiodic: $(OBJS2) $(INCL) ../../theory.options ../../common.mk Makefile 
	$(CC) $(OBJS2) $(CLINK) -o $@

%.o: %.c $(INCL) ../../theory.options ../../common.mk Makefile
	$(CC) $(OPT) $(GSL_CFLAGS) $(CFLAGS) $(INCLUDE) -c $< -o $@

python_lib: tests $(OBJS1) $(INCL) ../../theory.options ../../common.mk Makefile | ../../lib
	@echo 
	@echo "All THEORY tests are done. Now checking that the C extensions work."
	@echo 
	$(MAKE) -C ../python_bindings tests

tests: test_periodic test_nonperiodic
	./test_nonperiodic
	./test_periodic

../../lib:
	mkdir -p $@

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
	$(RM) $(targets) $(OBJS1) $(OBJS2)

clena: clean

celan: clean


