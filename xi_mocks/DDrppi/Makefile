ROOT_DIR := ../..
INSTALL_HEADERS_DIR := $(ROOT_DIR)/include
INSTALL_LIB_DIR := $(ROOT_DIR)/lib
INSTALL_BIN_DIR := $(ROOT_DIR)/bin
UTILS_DIR := $(ROOT_DIR)/utils
IO_DIR := $(ROOT_DIR)/io

include $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk

LIBNAME := countpairs_rp_pi_mocks
LIBRARY := lib$(LIBNAME).a
LIBSRC  := countpairs_rp_pi_mocks.c countpairs_rp_pi_mocks_impl_double.c countpairs_rp_pi_mocks_impl_float.c \
           $(UTILS_DIR)/gridlink_mocks_impl_float.c $(UTILS_DIR)/gridlink_mocks_impl_double.c \
           $(UTILS_DIR)/utils.c $(UTILS_DIR)/progressbar.c $(UTILS_DIR)/cpu_features.c \
	   $(UTILS_DIR)/set_cosmo_dist.c $(UTILS_DIR)/cosmology_params.c
LIBRARY_HEADERS := $(LIBNAME).h

GSL_TARGET   := DDrppi_mocks
TARGETS  := $(GSL_TARGET) wprp_mocks
TARGETSRC:= $(GSL_TARGET).c $(IO_DIR)/io.c $(IO_DIR)/ftread.c $(LIBSRC)
INCL     := countpairs_rp_pi_mocks_kernels_float.c countpairs_rp_pi_mocks_kernels_double.c countpairs_rp_pi_mocks_kernels.c.src \
            countpairs_rp_pi_mocks_impl.c.src countpairs_rp_pi_mocks_impl.h.src countpairs_rp_pi_mocks_impl_double.h countpairs_rp_pi_mocks_impl_float.h \
            countpairs_rp_pi_mocks.h \
            $(IO_DIR)/io.h $(IO_DIR)/ftread.h $(IO_DIR)/io.h  \
            $(UTILS_DIR)/gridlink_mocks_impl_double.h $(UTILS_DIR)/gridlink_mocks_impl_float.h $(UTILS_DIR)/gridlink_mocks_impl.h.src \
            $(UTILS_DIR)/cellarray_mocks_float.h $(UTILS_DIR)/cellarray_mocks_double.h $(UTILS_DIR)/cellarray_mocks.h.src \
	    $(UTILS_DIR)/set_cosmo_dist.h $(UTILS_DIR)/cosmology_params.h  $(UTILS_DIR)/progressbar.h $(UTILS_DIR)/cpu_features.h \
	    $(UTILS_DIR)/utils.h $(UTILS_DIR)/function_precision.h $(UTILS_DIR)/avx_calls.h $(UTILS_DIR)/defs.h

TARGETOBJS:=$(TARGETSRC:.c=.o)
LIBOBJS:=$(LIBSRC:.c=.o) 

all: $(TARGETS) $(TARGETSRC) $(INCL) $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk Makefile

$(GSL_TARGET): $(TARGETOBJS) $(INCL) $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk Makefile
	$(CC) $(TARGETOBJS) $(CLINK) $(GSL_LINK) -o $@ 

$(GSL_TARGET).o: $(GSL_TARGET).c $(ROOT_DIR)/common.mk Makefile $(INCL)
	$(CC) $(OPT) $(GSL_CFLAGS) $(CFLAGS) $(INCLUDE) -c $< -o $@

countpairs_rp_pi_mocks_impl_%.o:countpairs_rp_pi_mocks_impl_%.c countpairs_rp_pi_mocks_kernels_%.c $(INCL) $(ROOT_DIR)/common.mk $(ROOT_DIR)/mocks.options Makefile
	$(CC) $(GSL_CFLAGS) $(CFLAGS) $(INCLUDE) -c $< -o $@

wprp_mocks: wprp_mocks.c $(UTILS_DIR)/utils.c 
	$(CC) $(CFLAGS) $(INCLUDE) $^ $(CLINK) -o $@ 

# $(UTILS_DIR)/gridlink_mocks.o:$(UTILS_DIR)/gridlink_mocks.c $(ROOT_DIR)/common.mk $(ROOT_DIR)/mocks.options Makefile
# 	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

libs: lib
lib: $(LIBRARY) $(INSTALL_LIB_DIR)/$(LIBRARY) | $(INSTALL_LIB_DIR)

install: $(INSTALL_BIN_DIR)/DDrppi_mocks $(INSTALL_BIN_DIR)/wprp_mocks $(INSTALL_LIB_DIR)/$(LIBRARY) | $(INSTALL_BIN_DIR)

tests: 
	$(MAKE) -C ../tests DDrppi

clean:
	$(RM) $(TARGETS) $(TARGETOBJS) $(LIBRARY) countpairs_rp_pi_mocks_impl_float.[ch] countpairs_rp_pi_mocks_impl_double.[ch] countpairs_rp_pi_mocks_kernels_double.c countpairs_rp_pi_mocks_kernels_float.c 
	$(RM) -R *.dSYM
	$(MAKE) -C $(UTILS_DIR) clean

distclean:clean | $(INSTALL_LIB_DIR) $(INSTALL_HEADERS_DIR) $(INSTALL_BIN_DIR)
	cd $(INSTALL_LIB_DIR) && $(RM) $(LIBRARY)
	cd $(INSTALL_HEADERS_DIR) && $(RM) $(LIBRARY_HEADERS)
	cd $(INSTALL_BIN_DIR) && $(RM) $(TARGETS)

include $(ROOT_DIR)/rules.mk
