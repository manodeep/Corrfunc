ROOT_DIR := ../..
INSTALL_HEADERS_DIR := $(ROOT_DIR)/include
INSTALL_LIB_DIR := $(ROOT_DIR)/lib
INSTALL_BIN_DIR := $(ROOT_DIR)/bin
UTILS_DIR := $(ROOT_DIR)/utils
IO_DIR := $(ROOT_DIR)/io

MOCKS_DIR := $(ROOT_DIR)/xi_mocks
DDrppi_mocks_DIR := $(MOCKS_DIR)/DDrppi
WTHETA_DIR := $(MOCKS_DIR)/wtheta
VPF_mocks_DIR := $(MOCKS_DIR)/vpf

DDrppi_mocks_LIB := countpairs_rp_pi_mocks
WTHETA_LIB := countpairs_theta_mocks
VPF_mocks_LIB := countspheres_mocks

include $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk

TARGET := tests_mocks
TARGETS := $(TARGET)
ifneq ($(COMPILE_PYTHON_EXT), 0)
  TARGETS += python_lib
else
  $(warning $(ccmagenta) Skipping python tests since python or numpy is unavailable $(ccreset))
endif

TARGETSRC   := tests_mocks.c $(IO_DIR)/io.c $(IO_DIR)/ftread.c $(UTILS_DIR)/utils.c $(UTILS_DIR)/cosmology_params.c
TARGETOBJS  := $(TARGETSRC:.c=.o)
C_LIBRARIES := $(DDrppi_mocks_DIR)/lib$(DDrppi_mocks_LIB).a $(WTHETA_DIR)/lib$(WTHETA_LIB).a \
             $(VPF_mocks_DIR)/lib$(VPF_mocks_LIB).a
INCL   := $(IO_DIR)/io.h $(IO_DIR)/ftread.h $(UTILS_DIR)/utils.h \
          $(DDrppi_mocks_DIR)/$(DDrppi_mocks_LIB).h $(WTHETA_DIR)/$(WTHETA_LIB).h $(VPF_mocks_DIR)/$(VPF_mocks_LIB).h

EXTRA_INCL:=-DDOUBLE_PREC -I$(DDrppi_mocks_DIR) -I$(WTHETA_DIR) -I$(VPF_mocks_DIR) $(GSL_CFLAGS)
EXTRA_LINK := $(GSL_LINK)

OPT := 

all: tests $(TARGETS) $(INCL) uncompress $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk Makefile

$(TARGET):$(C_LIBRARIES)

$(DDrppi_mocks_DIR)/lib$(DDrppi_mocks_LIB).a: $(DDrppi_mocks_DIR)/*.c $(DDrppi_mocks_DIR)/*.c.src $(DDrppi_mocks_DIR)/*.h.src $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk
	$(MAKE) -C $(DDrppi_mocks_DIR) libs

$(WTHETA_DIR)/lib$(WTHETA_LIB).a: $(WTHETA_DIR)/*.c $(WTHETA_DIR)/*.c.src $(WTHETA_DIR)/*.h.src $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk
	$(MAKE) -C $(WTHETA_DIR) libs

$(VPF_mocks_DIR)/lib$(VPF_mocks_LIB).a: $(VPF_mocks_DIR)/*.c $(VPF_mocks_DIR)/*.c.src $(VPF_mocks_DIR)/*.h.src $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk
	$(MAKE) -C $(VPF_mocks_DIR) libs

python_lib: tests $(TARGETOBJS) $(INCL) $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk Makefile | $(ROOT_DIR)/lib
	@echo 
	@echo "All MOCKS tests are done. Now checking that the C extensions work."
	@echo 
	$(MAKE) -C ../python_bindings tests

tests: $(TARGET)
	./$(TARGET)

uncompress: | data
	@{\
		cd data ; \
		for f in $$(ls -f *.gz 2>/dev/null) ; do \
			gunzip $$f ; \
		done ;\
		cd ..; \
	}
DDrppi: $(TARGET)
	./$(TARGET) 0 3  

wtheta: tests_mocks
	./$(TARGET) 1 4 

vpf: tests_mocks
	./$(TARGET) 2 5

clean:
	$(RM) $(TARGETS) $(TARGETOBJS)
	$(RM) -R *.dSYM


include $(ROOT_DIR)/rules.mk
