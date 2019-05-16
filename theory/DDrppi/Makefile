ROOT_DIR := ../..
INSTALL_HEADERS_DIR := $(ROOT_DIR)/include
INSTALL_LIB_DIR := $(ROOT_DIR)/lib
INSTALL_BIN_DIR := $(ROOT_DIR)/bin
UTILS_DIR := $(ROOT_DIR)/utils
IO_DIR := $(ROOT_DIR)/io

include $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk

TARGET := DDrppi
TARGETS := $(TARGET) wprp
LIBRARY := libcountpairs_rp_pi.a
LIBSRC  := countpairs_rp_pi.c countpairs_rp_pi_impl_double.c countpairs_rp_pi_impl_float.c \
         $(UTILS_DIR)/gridlink_impl_double.c $(UTILS_DIR)/gridlink_impl_float.c \
         $(UTILS_DIR)/gridlink_utils_double.c $(UTILS_DIR)/gridlink_utils_float.c \
         $(UTILS_DIR)/utils.c $(UTILS_DIR)/progressbar.c $(UTILS_DIR)/cpu_features.c $(UTILS_DIR)/avx512_calls.c
LIBRARY_HEADERS := countpairs_rp_pi.h

TARGETSRC := DDrppi.c $(IO_DIR)/ftread.c $(IO_DIR)/io.c $(LIBSRC)
WPRPSRC := wprp.c $(UTILS_DIR)/utils.c

INCL   := countpairs_rp_pi_kernels_float.c countpairs_rp_pi_kernels_double.c countpairs_rp_pi_kernels.c.src countpairs_rp_pi_impl.c.src countpairs_rp_pi_impl.h.src \
          countpairs_rp_pi.h countpairs_rp_pi_impl_double.h countpairs_rp_pi_impl_float.h \
          $(UTILS_DIR)/gridlink_impl_float.h $(UTILS_DIR)/gridlink_impl_double.h $(UTILS_DIR)/gridlink_impl.h.src \
          $(UTILS_DIR)/gridlink_utils_float.h $(UTILS_DIR)/gridlink_utils_double.h $(UTILS_DIR)/gridlink_utils.h.src \
          $(UTILS_DIR)/cellarray_float.h $(UTILS_DIR)/cellarray_double.h $(UTILS_DIR)/cellarray.h.src \
          $(UTILS_DIR)/cell_pair_float.h $(UTILS_DIR)/cell_pair_double.h $(UTILS_DIR)/cell_pair.h.src \
          $(UTILS_DIR)/function_precision.h  $(UTILS_DIR)/avx512_calls.h $(UTILS_DIR)/avx_calls.h $(UTILS_DIR)/sse_calls.h \
          $(UTILS_DIR)/defs.h $(UTILS_DIR)/cpu_features.h \
          $(IO_DIR)/ftread.h $(IO_DIR)/io.h $(UTILS_DIR)/utils.h $(UTILS_DIR)/progressbar.h \
          $(UTILS_DIR)/weight_functions_double.h $(UTILS_DIR)/weight_functions_float.h $(UTILS_DIR)/weight_functions.h.src \
	  $(UTILS_DIR)/weight_defs_double.h $(UTILS_DIR)/weight_defs_float.h $(UTILS_DIR)/weight_defs.h.src

TARGETOBJS  := $(TARGETSRC:.c=.o)
LIBOBJS := $(LIBSRC:.c=.o)
all: $(TARGETS) $(TARGETSRC) $(WPRPSRC) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile

wprp: $(WPRPSRC) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $(WPRPSRC) $(CLINK)

countpairs_rp_pi_impl_float.o:countpairs_rp_pi_impl_float.c countpairs_rp_pi_impl_float.h countpairs_rp_pi_kernels_float.c $(UTILS_DIR)/gridlink_impl_float.h
countpairs_rp_pi_impl_double.o:countpairs_rp_pi_impl_double.c countpairs_rp_pi_impl_double.h countpairs_rp_pi_kernels_double.c $(UTILS_DIR)/gridlink_impl_double.h
countpairs_rp_pi.o:countpairs_rp_pi.c countpairs_rp_pi_impl_double.h countpairs_rp_pi_impl_float.h countpairs_rp_pi.h
countpairs_rp_pi_impl_float.c countpairs_rp_pi_impl_double.c:countpairs_rp_pi_impl.c.src $(INCL)

libs: lib
lib:  $(LIBRARY)
install:$(INSTALL_BIN_DIR)/$(TARGET) $(INSTALL_LIB_DIR)/$(LIBRARY) $(INSTALL_HEADERS_DIR)/$(LIBRARY_HEADERS)

clean:
	$(RM) $(TARGETOBJS) $(TARGET) $(LIBRARY) countpairs_rp_pi_kernels_float.c countpairs_rp_pi_kernels_double.c countpairs_rp_pi_impl_double.[ch] countpairs_rp_pi_impl_float.[ch]
	$(RM) -R *.dSYM

distclean:clean
	cd $(INSTALL_HEADERS_DIR) && $(RM) $(LIBRARY_HEADERS)
	cd $(INSTALL_LIB_DIR) && $(RM) $(LIBRARY)
	cd $(INSTALL_BIN_DIR) && $(RM) $(TARGET)

tests:
	$(MAKE) -C ../tests DDrppi

include $(ROOT_DIR)/rules.mk
