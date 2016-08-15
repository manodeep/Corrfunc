ROOT_DIR := ../..
UTILS_DIR := $(ROOT_DIR)/utils
IO_DIR := $(ROOT_DIR)/io
INCLUDE_DIR := $(ROOT_DIR)/include
LIB_DIR := $(ROOT_DIR)/lib
DATA_DIR := ../tests/data

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

TARGET := run_correlations
TARGETSRC   := run_correlations.c $(IO_DIR)/ftread.c $(IO_DIR)/io.c  $(UTILS_DIR)/utils.c $(UTILS_DIR)/progressbar.c 
TARGETOBJS  := $(TARGETSRC:.c=.o)
C_LIBRARY := $(LIB_DIR)/lib$(DD_LIB).a $(LIB_DIR)/lib$(DDrppi_LIB).a $(LIB_DIR)/lib$(WP_LIB).a \
             $(LIB_DIR)/lib$(XI_LIB).a $(LIB_DIR)/lib$(VPF_LIB).a

INCL :=	$(INCLUDE_DIR)/$(DD_LIB).h $(INCLUDE_DIR)/$(DDrppi_LIB).h $(INCLUDE_DIR)/$(WP_LIB).h $(INCLUDE_DIR)/$(XI_LIB).h $(INCLUDE_DIR)/$(VPF_LIB).h \
	$(UTILS_DIR)/defs.h 

all: $(TARGET) $(TARGETSRC) $(INCL) $(C_LIBRARY) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile uncompress

$(TARGET): Makefile $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk $(INCL) $(C_LIBRARY) $(TARGETOBJS) 
	$(CC) $(TARGETOBJS) $(C_LIBRARY) $(CLINK) $(GSL_LINK) -o $@

$(INCLUDE_DIR)/$(DD_LIB).h $(LIB_DIR)/lib$(DD_LIB).a: $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile 
	$(MAKE) -C $(DD_DIR) install

$(INCLUDE_DIR)/$(DDrppi_LIB).h $(LIB_DIR)/lib$(DDrppi_LIB).a: $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile
	$(MAKE) -C $(DDrppi_DIR) install

$(INCLUDE_DIR)/$(WP_LIB).h $(LIB_DIR)/lib$(WP_LIB).a: $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile
	$(MAKE) -C $(WP_DIR) install

$(INCLUDE_DIR)/$(XI_LIB).h $(LIB_DIR)/lib$(XI_LIB).a: $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile
	$(MAKE) -C $(XI_DIR) install

$(INCLUDE_DIR)/$(VPF_LIB).h $(LIB_DIR)/lib$(VPF_LIB).a: $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile
	$(MAKE) -C $(VPF_DIR) install

%.o: %.c $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile
	$(CC) $(CFLAGS) -I$(INCLUDE_DIR) $(OPT) $(INCLUDE) -c $< -o $@

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




