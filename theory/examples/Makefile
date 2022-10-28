ROOT_DIR := ../..
UTILS_DIR := $(ROOT_DIR)/utils
IO_DIR := $(ROOT_DIR)/io
INCLUDE_DIR := $(ROOT_DIR)/include
LIB_DIR := $(ROOT_DIR)/lib
DATA_DIR := ../tests/data

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

TARGET := run_correlations
TARGETSRC   := run_correlations.c $(IO_DIR)/ftread.c $(IO_DIR)/io.c  $(UTILS_DIR)/utils.c $(UTILS_DIR)/progressbar.c 
TARGETOBJS  := $(TARGETSRC:.c=.o)
C_LIBRARIES := $(DD_DIR)/lib$(DD_LIB).a $(DDrppi_DIR)/lib$(DDrppi_LIB).a $(DDsmu_DIR)/lib$(DDsmu_LIB).a $(WP_DIR)/lib$(WP_LIB).a \
             $(XI_DIR)/lib$(XI_LIB).a $(VPF_DIR)/lib$(VPF_LIB).a
INCL :=	$(DD_DIR)/$(DD_LIB).h $(DDrppi_DIR)/$(DDrppi_LIB).h $(DDsmu_DIR)/$(DDsmu_LIB).h $(WP_DIR)/$(WP_LIB).h \
        $(XI_DIR)/$(XI_LIB).h $(VPF_DIR)/$(VPF_LIB).h \
        $(UTILS_DIR)/defs.h 

EXTRA_INCL:= -I$(DD_DIR) -I$(DDrppi_DIR) -I$(DDsmu_DIR) -I$(WP_DIR) -I$(XI_DIR) -I$(VPF_DIR) $(GSL_CFLAGS)
EXTRA_LINK:= $(GSL_LINK)

all: $(TARGET) $(TARGETSRC) $(INCL) $(C_LIBRARY) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile uncompress

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

$(TARGET): $(C_LIBRARIES) $(TARGETOBJS)

uncompress: | $(DATA_DIR)
	@{\
		oldpwd=$$(pwd); \
		cd $(DATA_DIR) ; \
		for f in $$(ls -f *.gz 2>/dev/null) ; do \
			gunzip $$f ; \
		done ;\
		cd $$oldpwd; \
	}

.PHONY: clean clena celan uncompress

clena: clean
celan: clean
clean:
	$(RM) $(TARGET) $(TARGETOBJS)
	$(RM) -R *.dSYM

include $(ROOT_DIR)/rules.mk

