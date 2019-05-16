ROOT_DIR := ../..
INSTALL_HEADERS_DIR := $(ROOT_DIR)/include
INSTALL_LIB_DIR := $(ROOT_DIR)/lib
INSTALL_BIN_DIR := $(ROOT_DIR)/bin
UTILS_DIR := $(ROOT_DIR)/utils
IO_DIR := $(ROOT_DIR)/io

include $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk
LIBNAME:= countpairs_wp
LIBRARY := lib$(LIBNAME).a
LIBSRC := countpairs_wp.c countpairs_wp_impl_double.c countpairs_wp_impl_float.c \
         $(UTILS_DIR)/gridlink_impl_double.c $(UTILS_DIR)/gridlink_impl_float.c \
         $(UTILS_DIR)/gridlink_utils_float.c $(UTILS_DIR)/gridlink_utils_double.c \
         $(UTILS_DIR)/utils.c $(UTILS_DIR)/progressbar.c $(UTILS_DIR)/cpu_features.c $(UTILS_DIR)/avx512_calls.c
LIBRARY_HEADERS := $(LIBNAME).h

TARGET := wp
TARGETSRC := wp.c $(IO_DIR)/io.c $(IO_DIR)/ftread.c $(LIBSRC)

INCL   := wp_kernels_double.c wp_kernels_float.c wp_kernels.c.src countpairs_wp.h \
          countpairs_wp_impl_float.h countpairs_wp_impl_double.h countpairs_wp_impl.h.src \
          $(UTILS_DIR)/gridlink_impl_float.h $(UTILS_DIR)/gridlink_impl_double.h $(UTILS_DIR)/gridlink_impl.h.src \
          $(UTILS_DIR)/gridlink_utils_float.h $(UTILS_DIR)/gridlink_utils_double.h $(UTILS_DIR)/gridlink_utils.h.src \
          $(UTILS_DIR)/cellarray_double.h $(UTILS_DIR)/cellarray_float.h $(UTILS_DIR)/cellarray.h.src \
          $(UTILS_DIR)/cell_pair_float.h $(UTILS_DIR)/cell_pair_double.h $(UTILS_DIR)/cell_pair.h.src \
          $(UTILS_DIR)/avx512_calls.h $(UTILS_DIR)/avx2_calls.h \
          $(UTILS_DIR)/avx_calls.h $(UTILS_DIR)/sse_calls.h $(UTILS_DIR)/function_precision.h $(UTILS_DIR)/defs.h \
          $(UTILS_DIR)/cpu_features.h \
          $(IO_DIR)/ftread.h $(IO_DIR)/io.h $(UTILS_DIR)/utils.h $(UTILS_DIR)/sglib.h $(UTILS_DIR)/progressbar.h \
	  $(UTILS_DIR)/weight_functions_double.h $(UTILS_DIR)/weight_functions_float.h $(UTILS_DIR)/weight_functions.h.src \
	  $(UTILS_DIR)/weight_defs_double.h $(UTILS_DIR)/weight_defs_float.h $(UTILS_DIR)/weight_defs.h.src


TARGETOBJS  := $(TARGETSRC:.c=.o)
LIBOBJS := $(LIBSRC:.c=.o)

all: $(TARGET) $(TARGETOBJS) $(TARGETSRC) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile

countpairs_wp_impl_float.o:countpairs_wp_impl_float.c countpairs_wp_impl_float.h wp_kernels_float.c $(UTILS_DIR)/gridlink_impl_float.h
countpairs_wp_impl_double.o:countpairs_wp_impl_double.c countpairs_wp_impl_double.h wp_kernels_double.c $(UTILS_DIR)/gridlink_impl_double.h
countpairs_wp.o:countpairs_wp.c countpairs_wp_impl_double.h countpairs_wp_impl_float.h countpairs_wp.h
countpairs_wp_impl_float.c countpairs_wp_impl_double.c:countpairs_wp_impl.c.src $(INCL)



libs: lib
lib: $(LIBRARY)
install:$(INSTALL_BIN_DIR)/$(TARGET) $(INSTALL_LIB_DIR)/$(LIBRARY) $(INSTALL_HEADERS_DIR)/$(LIBRARY_HEADERS)
clean:
	$(RM) $(TARGETOBJS) $(TARGET) $(LIBRARY) wp_kernels_float.c wp_kernels_double.c countpairs_wp_impl_double.[ch] countpairs_wp_impl_float.[ch]
	$(RM) -R *.dSYM

distclean:clean
	cd $(INSTALL_HEADERS_DIR) && $(RM) $(LIBRARY_HEADERS)
	cd $(INSTALL_LIB_DIR) && $(RM) $(LIBRARY)
	cd $(INSTALL_BIN_DIR) && $(RM) $(TARGET)

tests: tests_periodic

tests_periodic:
	$(MAKE) -C ../tests wp

include $(ROOT_DIR)/rules.mk
