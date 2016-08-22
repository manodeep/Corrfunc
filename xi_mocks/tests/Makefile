ROOT_DIR := ../..
INSTALL_HEADERS_DIR := $(ROOT_DIR)/include
INSTALL_LIB_DIR := $(ROOT_DIR)/lib
INSTALL_BIN_DIR := $(ROOT_DIR)/bin
UTILS_DIR := $(ROOT_DIR)/utils
IO_DIR := $(ROOT_DIR)/io

MOCKS_DIR := $(ROOT_DIR)/xi_mocks
DDrppi_DIR := $(MOCKS_DIR)/DDrppi
WTHETA_DIR := $(MOCKS_DIR)/wtheta
VPF_DIR := $(MOCKS_DIR)/vpf

DDrppi_LIB := countpairs_rp_pi_mocks
WTHETA_LIB := countpairs_theta_mocks
VPF_LIB := countspheres_mocks

include $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk

TARGETS := test_mocks
ifneq ($(COMPILE_PYTHON_EXT), 0)
TARGETS += python_lib
else
$(warning $(ccmagenta) Skipping python tests since python or numpy is unavailable $(ccreset))
endif

SRC1   := tests_mocks.c $(IO_DIR)/io.c $(IO_DIR)/ftread.c $(UTILS_DIR)/utils.c $(UTILS_DIR)/cosmology_params.c
OBJS1  := $(SRC1:.c=.o)
C_LIBRARIES := $(DDrppi_DIR)/lib$(DDrppi_LIB).a $(WTHETA_DIR)/lib$(WTHETA_LIB).a \
             $(VPF_DIR)/lib$(VPF_LIB).a
INCL   := $(IO_DIR)/io.h $(IO_DIR)/ftread.h $(UTILS_DIR)/utils.h \
          $(DDrppi_DIR)/$(DDrppi_LIB).h $(WTHETA_DIR)/$(WTHETA_LIB).h $(VPF_DIR)/$(VPF_LIB).h

all: tests $(TARGETS) $(INCL) uncompress $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk Makefile

test_mocks: $(OBJS1) $(INCL) Makefile 
	$(CC) $(OBJS1) $(C_LIBRARIES) $(CLINK) $(GSL_LINK) -lgsl -lgslcblas -o $@

%.o: %.c $(INCL) $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk Makefile
	$(CC) -DDOUBLE_PREC $(GSL_CFLAGS) $(CFLAGS) $(INCLUDE) -c $< -o $@

$(DDrppi_DIR)/$(DDrppi_LIB).h $(DDrppi_DIR)/lib$(DDrppi_LIB).a: $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk Makefile
	$(MAKE) -C $(DDrppi_DIR) libs

$(WTHETA_DIR)/$(WTHETA_LIB).h $(WTHETA_DIR)/lib$(WTHETA_LIB).a: $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk Makefile
	$(MAKE) -C $(WTHETA_DIR) libs

$(VPF_DIR)/$(VPF_LIB).h $(VPF_DIR)/lib$(VPF_LIB).a: $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk Makefile
	$(MAKE) -C $(VPF_DIR) libs


python_lib: tests $(OBJS1) $(INCL) $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk Makefile | $(ROOT_DIR)/lib
	@echo 
	@echo "All MOCKS tests are done. Now checking that the C extensions work."
	@echo 
	$(MAKE) -C ../python_bindings tests

tests: test_mocks
	./test_mocks

uncompress: | data
	@{\
		cd data ; \
		for f in $$(ls -f *.gz 2>/dev/null) ; do \
			gunzip $$f ; \
		done ;\
		cd ..; \
	}
DDrppi: tests_mocks
	./tests_mocks 0 3  

wtheta: tests_mocks
	./tests_mocks 1 4 

vpf: tests_mocks
	./tests_mocks 2 5

clean:
	$(RM) $(targets) $(OBJS1) 


