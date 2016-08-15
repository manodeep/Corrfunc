ROOT_DIR := ../..
INSTALL_HEADERS_DIR := $(ROOT_DIR)/include
INSTALL_LIB_DIR := $(ROOT_DIR)/lib
INSTALL_BIN_DIR := $(ROOT_DIR)/bin
UTILS_DIR := $(ROOT_DIR)/utils
IO_DIR := $(ROOT_DIR)/io

THEORY_DIR := $(ROOT_DIR)/xi_theory
DD_DIR := $(THEORY_DIR)/xi_of_r
DDrppi_DIR := $(THEORY_DIR)/xi_rp_pi
WP_DIR := $(THEORY_DIR)/wp
XI_DIR := $(THEORY_DIR)/xi
VPF_DIR := $(THEORY_DIR)/vpf

DD_LIB := countpairs
DDrppi_LIB := countpairs_rp_pi
WP_LIB := countpairs_wp
XI_LIB := countpairs_xi
VPF_LIB := countspheres

include $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk

TARGETS := test_periodic test_nonperiodic
ifneq ($(COMPILE_PYTHON_EXT), 0)
TARGETS += python_lib
else
$(warning $(ccmagenta) Skipping python tests since python or numpy is unavailable $(ccreset))
endif

SRC1   := test_periodic.c $(IO_DIR)/io.c $(IO_DIR)/ftread.c $(UTILS_DIR)/utils.c
OBJS1  := $(SRC1:.c=.o)
C_LIBRARY := $(DD_DIR)/lib$(DD_LIB).a $(DDrppi_DIR)/lib$(DDrppi_LIB).a $(WP_DIR)/lib$(WP_LIB).a \
             $(XI_DIR)/lib$(XI_LIB).a $(VPF_DIR)/lib$(VPF_LIB).a
INCL   := $(IO_DIR)/io.h $(IO_DIR)/ftread.h $(UTILS_DIR)/utils.h \
          $(DD_DIR)/$(DD_LIB).h $(DDrppi_DIR)/$(DDrppi_LIB).h $(WP_DIR)/$(WP_LIB).h \
          $(XI_DIR)/$(XI_LIB).h $(VPF_DIR)/$(VPF_LIB).h

SRC2   := test_nonperiodic.c $(UTILS_DIR)/utils.c $(IO_DIR)/io.c $(IO_DIR)/ftread.c
OBJS2  := $(SRC2:.c=.o)

all: tests $(TARGETS) $(INCL) uncompress $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile

test_periodic: $(OBJS1) $(C_LIBRARY) $(INCL) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile 
	$(CC) $(OBJS1) $(C_LIBRARY) $(GSL_LINK) $(CLINK) -o $@

test_nonperiodic: $(OBJS2) $(C_LIBRARY) $(INCL) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile 
	$(CC) $(OBJS2) $(C_LIBRARY) $(CLINK) -o $@

%.o: %.c $(INCL) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile
	$(CC) $(GSL_CFLAGS) $(CFLAGS) $(INCLUDE) -c $< -o $@

python_lib: tests $(OBJS1) $(INCL) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile | $(ROOT_DIR)/lib
	@echo 
	@echo "All THEORY tests are done. Now checking that the C extensions work."
	@echo 
	$(MAKE) -C ../python_bindings tests

tests: test_periodic test_nonperiodic
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

clean:
	$(RM) $(targets) $(OBJS1) $(OBJS2)


