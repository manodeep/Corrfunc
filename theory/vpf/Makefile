ROOT_DIR := ../..
INSTALL_HEADERS_DIR := $(ROOT_DIR)/include
INSTALL_LIB_DIR := $(ROOT_DIR)/lib
INSTALL_BIN_DIR := $(ROOT_DIR)/bin
UTILS_DIR := $(ROOT_DIR)/utils
IO_DIR := $(ROOT_DIR)/io

GSL_REQUIRED := true
include $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk

LIBRARY := libcountspheres.a
LIBRARY_HEADERS := countspheres.h
LIBSRC := countspheres.c countspheres_impl_double.c countspheres_impl_float.c \
          $(UTILS_DIR)/gridlink_impl_double.c $(UTILS_DIR)/gridlink_impl_float.c \
          $(UTILS_DIR)/gridlink_utils_float.c $(UTILS_DIR)/gridlink_utils_double.c \
          $(UTILS_DIR)/utils.c $(UTILS_DIR)/progressbar.c \
          $(UTILS_DIR)/cpu_features.c $(UTILS_DIR)/avx512_calls.c

TARGET := vpf
TARGETSRC := $(TARGET).c $(IO_DIR)/ftread.c $(IO_DIR)/io.c $(LIBSRC)
INCL   := vpf_kernels_float.c vpf_kernels_double.c vpf_kernels.c.src countspheres.h \
          countspheres_impl_float.h countspheres_impl_double.h countspheres_impl.h.src countspheres_impl.c.src \
          $(UTILS_DIR)/gridlink_impl_float.h $(UTILS_DIR)/gridlink_impl_double.h $(UTILS_DIR)/gridlink_impl.h.src \
          $(UTILS_DIR)/gridlink_utils_float.h $(UTILS_DIR)/gridlink_utils_double.h $(UTILS_DIR)/gridlink_utils.h.src \
          $(UTILS_DIR)/cellarray_double.h $(UTILS_DIR)/cellarray_float.h $(UTILS_DIR)/cellarray.h.src \
          $(IO_DIR)/ftread.h $(IO_DIR)/io.h $(UTILS_DIR)/utils.h \
          $(UTILS_DIR)/avx512_calls.h $(UTILS_DIR)/avx_calls.h $(UTILS_DIR)/sse_calls.h \
	  $(UTILS_DIR)/function_precision.h $(UTILS_DIR)/defs.h $(UTILS_DIR)/sglib.h $(UTILS_DIR)/progressbar.h \
          $(UTILS_DIR)/cpu_features.h

TARGETOBJS  := $(TARGETSRC:.c=.o)
LIBOBJS := $(LIBSRC:.c=.o)

EXTRA_INCL:=$(GSL_CFLAGS)
EXTRA_LINK:=$(GSL_LINK)

all: $(TARGET) $(TARGETSRC) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile
libs: lib
lib: $(LIBRARY)

countspheres_impl_float.o:countspheres_impl_float.c countspheres_impl_float.h vpf_kernels_float.c $(UTILS_DIR)/cellarray_float.h \
                          $(UTILS_DIR)/gridlink_impl_float.h  $(UTILS_DIR)/gridlink_utils_float.h
countspheres_impl_double.o:countspheres_impl_double.c countspheres_impl_double.h vpf_kernels_double.c $(UTILS_DIR)/cellarray_double.h \
                           $(UTILS_DIR)/gridlink_impl_double.h $(UTILS_DIR)/gridlink_utils_double.h
countspheres.o:countspheres.c countspheres_impl_float.h countspheres_impl_double.h

tests:
	$(MAKE) -C ../tests vpf


install:$(INSTALL_BIN_DIR)/$(TARGET) $(INSTALL_LIB_DIR)/$(LIBRARY) $(INSTALL_HEADERS_DIR)/$(LIBRARY_HEADERS)

clean:
	$(RM) $(TARGETOBJS) $(LIBRARY) $(TARGET) vpf_kernels_float.c vpf_kernels_double.c countspheres_impl_double.[ch] countspheres_impl_float.[ch]
	$(RM) -R *.dSYM

distclean:clean
	cd $(INSTALL_HEADERS_DIR) && $(RM) $(LIBRARY_HEADERS)
	cd $(INSTALL_LIB_DIR) && $(RM) $(LIBRARY)
	cd $(INSTALL_BIN_DIR) && $(RM) $(TARGET)

include $(ROOT_DIR)/rules.mk
