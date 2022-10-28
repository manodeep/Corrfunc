ROOT_DIR := ../..
INSTALL_HEADERS_DIR := $(ROOT_DIR)/include
INSTALL_LIB_DIR := $(ROOT_DIR)/lib
INSTALL_BIN_DIR := $(ROOT_DIR)/bin
UTILS_DIR := $(ROOT_DIR)/utils
IO_DIR := $(ROOT_DIR)/io

THEORY_DIR := $(ROOT_DIR)/theory
DD_DIR := $(THEORY_DIR)/DD
DDrppi_DIR := $(THEORY_DIR)/DDrppi
DDsmu_DIR := $(THEORY_DIR)/DDsmu
WP_DIR := $(THEORY_DIR)/wp
XI_DIR := $(THEORY_DIR)/xi
VPF_DIR := $(THEORY_DIR)/vpf

DD_LIB := countpairs
DDrppi_LIB := countpairs_rp_pi
DDsmu_LIB := countpairs_s_mu
WP_LIB := countpairs_wp
XI_LIB := countpairs_xi
VPF_LIB := countspheres

GSL_REQUIRED := true
include $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk

TARGETS := test_periodic test_nonperiodic

ifneq ($(COMPILE_PYTHON_EXT), 0)
  TARGETS += python_lib
else
  $(warning $(ccmagenta) Skipping python tests since python or numpy is unavailable $(ccreset))
endif

SRC1   := test_periodic.c $(IO_DIR)/io.c $(IO_DIR)/ftread.c $(UTILS_DIR)/utils.c
OBJS1  := $(SRC1:.c=.o)
C_LIBRARIES := $(DD_DIR)/lib$(DD_LIB).a $(DDrppi_DIR)/lib$(DDrppi_LIB).a $(DDsmu_DIR)/lib$(DDsmu_LIB).a $(WP_DIR)/lib$(WP_LIB).a \
             $(XI_DIR)/lib$(XI_LIB).a $(VPF_DIR)/lib$(VPF_LIB).a
INCL   := $(IO_DIR)/io.h $(IO_DIR)/ftread.h $(UTILS_DIR)/utils.h $(UTILS_DIR)/tests_common.h \
          $(DD_DIR)/$(DD_LIB).h $(DDrppi_DIR)/$(DDrppi_LIB).h $(DDsmu_DIR)/$(DDsmu_LIB).h  $(WP_DIR)/$(WP_LIB).h \
          $(XI_DIR)/$(XI_LIB).h $(VPF_DIR)/$(VPF_LIB).h

SRC2   := test_nonperiodic.c $(UTILS_DIR)/utils.c $(IO_DIR)/io.c $(IO_DIR)/ftread.c
OBJS2  := $(SRC2:.c=.o)

all: tests $(TARGETS) $(INCL) uncompress $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile

test_periodic: $(OBJS1) $(C_LIBRARIES) $(INCL) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile
	$(CC) $(OBJS1) $(C_LIBRARIES) $(GSL_LINK) $(CLINK) -o $@

test_nonperiodic: $(OBJS2) $(C_LIBRARIES) $(INCL) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile
	$(CC) $(OBJS2) $(C_LIBRARIES) $(CLINK) -o $@

%.o: %.c $(INCL) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile
	$(CC) $(GSL_CFLAGS) $(CFLAGS) $(INCLUDE) -c $< -o $@

$(DD_DIR)/lib$(DD_LIB).a: $(DD_DIR)/*.c $(DD_DIR)/*.c.src $(DD_DIR)/*.h.src $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk
	$(MAKE) -C $(DD_DIR) libs

$(DDrppi_DIR)/lib$(DDrppi_LIB).a: $(DDrppi_DIR)/*.c $(DDrppi_DIR)/*.c.src $(DDrppi_DIR)/*.h.src $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk
	$(MAKE) -C $(DDrppi_DIR) libs

$(DDsmu_DIR)/lib$(DDsmu_LIB).a: $(DDsmu_DIR)/*.c $(DDsmu_DIR)/*.c.src $(DDsmu_DIR)/*.h.src $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk
	$(MAKE) -C $(DDsmu_DIR) libs

$(WP_DIR)/lib$(WP_LIB).a: $(WP_DIR)/*.c $(WP_DIR)/*.c.src $(WP_DIR)/*.h.src $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk
	$(MAKE) -C $(WP_DIR) libs

$(XI_DIR)/lib$(XI_LIB).a: $(XI_DIR)/*.c $(XI_DIR)/*.c.src $(XI_DIR)/*.h.src $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk
	$(MAKE) -C $(XI_DIR) libs

$(VPF_DIR)/lib$(VPF_LIB).a: $(VPF_DIR)/*.c $(VPF_DIR)/*.c.src $(VPF_DIR)/*.h.src $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk
	$(MAKE) -C $(VPF_DIR) libs

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

DDrppi: test_periodic test_nonperiodic
	./test_periodic 0 6 7 8
	./test_nonperiodic 1 3

DDsmu: test_periodic test_nonperiodic
	./test_periodic 5
	./test_nonperiodic 2

DD: test_periodic test_nonperiodic
	./test_periodic 1
	./test_nonperiodic 0

vpf: test_periodic
	./test_periodic 3

xi: test_periodic
	./test_periodic 4

.PHONY: celna clena celan clean

celna clena celan: clean
clean:
	$(RM) $(TARGETS) $(OBJS1) $(OBJS2)
	$(RM) -R *.dSYM
