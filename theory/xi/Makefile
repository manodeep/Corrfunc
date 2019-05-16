ROOT_DIR := ../..
INSTALL_HEADERS_DIR := $(ROOT_DIR)/include
INSTALL_LIB_DIR := $(ROOT_DIR)/lib
INSTALL_BIN_DIR := $(ROOT_DIR)/bin
UTILS_DIR := $(ROOT_DIR)/utils
IO_DIR := $(ROOT_DIR)/io

include $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk

LIBNAME:= countpairs_xi
LIBRARY := lib$(LIBNAME).a
LIBRARY_HEADERS := $(LIBNAME).h
LIBSRC := countpairs_xi.c countpairs_xi_impl_double.c countpairs_xi_impl_float.c  \
          $(UTILS_DIR)/gridlink_impl_double.c $(UTILS_DIR)/gridlink_impl_float.c \
          $(UTILS_DIR)/gridlink_utils_float.c $(UTILS_DIR)/gridlink_utils_double.c \
          $(UTILS_DIR)/utils.c $(UTILS_DIR)/progressbar.c $(UTILS_DIR)/cpu_features.c $(UTILS_DIR)/avx512_calls.c

TARGET := xi
TARGETSRC := $(TARGET).c $(IO_DIR)/ftread.c $(IO_DIR)/io.c $(LIBSRC)

INCL   := xi_kernels_float.c xi_kernels_double.c xi_kernels.c.src countpairs_xi.h \
          countpairs_xi_impl_float.h countpairs_xi_impl_double.h countpairs_xi_impl.h.src countpairs_xi_impl.c.src \
          $(UTILS_DIR)/gridlink_impl_float.h $(UTILS_DIR)/gridlink_impl_double.h $(UTILS_DIR)/gridlink_impl.h.src \
          $(UTILS_DIR)/gridlink_utils_float.h $(UTILS_DIR)/gridlink_utils_double.h $(UTILS_DIR)/gridlink_utils.h.src \
          $(UTILS_DIR)/cellarray_double.h $(UTILS_DIR)/cellarray_float.h $(UTILS_DIR)/cellarray.h.src \
          $(UTILS_DIR)/cell_pair_float.h $(UTILS_DIR)/cell_pair_double.h $(UTILS_DIR)/cell_pair.h.src \
          $(UTILS_DIR)/avx512_calls.h $(UTILS_DIR)/avx_calls.h $(UTILS_DIR)/sse_calls.h \
          $(IO_DIR)/ftread.h $(IO_DIR)/io.h $(UTILS_DIR)/utils.h \
          $(UTILS_DIR)/function_precision.h $(UTILS_DIR)/defs.h $(UTILS_DIR)/sglib.h $(UTILS_DIR)/progressbar.h \
          $(UTILS_DIR)/weight_functions_double.h $(UTILS_DIR)/weight_functions_float.h $(UTILS_DIR)/weight_functions.h.src \
	  $(UTILS_DIR)/weight_defs_double.h $(UTILS_DIR)/weight_defs_float.h $(UTILS_DIR)/weight_defs.h.src


TARGETOBJS  := $(TARGETSRC:.c=.o)
LIBOBJS := $(LIBSRC:.c=.o)

all: $(TARGET) $(TARGETSRC) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile


countpairs_xi_impl_float.o:countpairs_xi_impl_float.c countpairs_xi_impl_float.h xi_kernels_float.c \
                           $(UTILS_DIR)/gridlink_impl_float.h $(UTILS_DIR)/gridlink_utils_float.h
countpairs_xi_impl_double.o:countpairs_xi_impl_double.c countpairs_xi_impl_double.h xi_kernels_double.c\
                           $(UTILS_DIR)/gridlink_impl_double.h $(UTILS_DIR)/gridlink_utils_float.h
countpairs_xi.o:countpairs_xi.c countpairs_xi_impl_double.h countpairs_xi_impl_float.h countpairs_xi.h
countpairs_xi_impl_float.c countpairs_xi_impl_double.c:countpairs_xi_impl.c.src $(INCL)

libs: lib
lib: $(LIBRARY)
install:$(INSTALL_BIN_DIR)/$(TARGET) $(INSTALL_LIB_DIR)/$(LIBRARY) $(INSTALL_HEADERS_DIR)/$(LIBRARY_HEADERS)

clean:
	$(RM) $(LIBRARY) $(TARGETOBJS) $(TARGET) xi_kernels_float.c xi_kernels_double.c countpairs_xi_impl_float.[ch] countpairs_xi_impl_double.[ch]
	$(RM) -R *.dSYM

distclean:clean
	cd $(INSTALL_HEADERS_DIR) && $(RM) $(LIBRARY_HEADERS)
	cd $(INSTALL_LIB_DIR) && $(RM) $(LIBRARY)
	cd $(INSTALL_BIN_DIR) && $(RM) $(TARGET)

tests: tests_periodic

tests_periodic:
	$(MAKE) -C ../tests xi

include $(ROOT_DIR)/rules.mk
