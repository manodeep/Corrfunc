ROOT_DIR := ../..
UTILS_DIR := $(ROOT_DIR)/utils
IO_DIR := $(ROOT_DIR)/io
INCLUDE_DIR := $(ROOT_DIR)/include
LIB_DIR := $(ROOT_DIR)/lib
DATA_DIR := ../tests/data

MOCKS_DIR := $(ROOT_DIR)/xi_mocks
DDrppi_DIR := $(MOCKS_DIR)/DDrppi
WTHETA_DIR := $(MOCKS_DIR)/wtheta
VPF_DIR := $(MOCKS_DIR)/vpf

DDrppi_LIB := countpairs_rp_pi_mocks
WTHETA_LIB := countpairs_theta_mocks
VPF_LIB := countspheres_mocks

include $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk

TARGET := run_correlations_mocks
TARGETSRC   := $(TARGET).c $(IO_DIR)/ftread.c $(IO_DIR)/io.c  $(UTILS_DIR)/utils.c $(UTILS_DIR)/progressbar.c \
               $(UTILS_DIR)/cosmology_params.c 
TARGETOBJS  := $(TARGETSRC:.c=.o)
C_LIBRARIES := $(LIB_DIR)/lib$(DDrppi_LIB).a $(LIB_DIR)/lib$(WTHETA_LIB).a $(LIB_DIR)/lib$(VPF_LIB).a
INCL :=	$(INCLUDE_DIR)/$(DDrppi_LIB).h $(INCLUDE_DIR)/$(WTHETA_LIB).h $(INCLUDE_DIR)/$(VPF_LIB).h \
        $(UTILS_DIR)/defs.h $(IO_DIR)/io.h $(IO_DIR)/ftread.h \
        $(UTILS_DIR)/utils.h $(UTILS_DIR)/gridlink_mocks_impl_double.h $(UTILS_DIR)/gridlink_mocks_impl_float.h \
	$(UTILS_DIR)/function_precision.h $(UTILS_DIR)/cellarray_mocks_double.h $(UTILS_DIR)/cellarray_mocks_float.h \
        $(UTILS_DIR)/progressbar.h $(UTILS_DIR)/cosmology_params.h 

all: $(TARGET) $(TARGETSRC) $(INCL) $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk Makefile $(C_LIBRARIES)

$(TARGET): Makefile $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk $(INCL) $(C_LIBRARIES) $(TARGETOBJS) 
	$(CC) $(OPT) $(TARGETOBJS) $(C_LIBRARIES) $(CLINK) $(GSL_LINK) -o $@

$(INCLUDE_DIR)/$(DDrppi_LIB).h $(LIB_DIR)/lib$(DDrppi_LIB).a: $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk Makefile
	$(MAKE) -C $(DDrppi_DIR) install

$(INCLUDE_DIR)/$(WTHETA_LIB).h $(LIB_DIR)/lib$(WTHETA_LIB).a: $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk Makefile
	$(MAKE) -C $(WTHETA_DIR) install

$(INCLUDE_DIR)/$(VPF_LIB).h $(LIB_DIR)/lib$(VPF_LIB).a: $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk Makefile
	$(MAKE) -C $(VPF_DIR) install

$(TARGET).o: $(TARGET).c $(ROOT_DIR)/common.mk $(ROOT_DIR)/mocks.options Makefile
	$(CC) $(CFLAGS) -I$(ROOT_DIR)/include $(INCLUDE) $(OPT) -c $< -o $@

clean:
	$(RM) $(TARGET) $(TARGETOBJS)
distclean:clean



