ROOT_DIR:=../..
INSTALL_HEADERS_DIR := $(ROOT_DIR)/include
INSTALL_LIB_DIR := $(ROOT_DIR)/lib
INSTALL_BIN_DIR := $(ROOT_DIR)/bin
UTILS_DIR:=$(ROOT_DIR)/utils
IO_DIR:=$(ROOT_DIR)/io

GSL_REQUIRED := true
include $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk

LIBNAME:=countspheres_mocks
LIBRARY:=lib$(LIBNAME).a
LIBSRC:=countspheres_mocks.c countspheres_mocks_impl_float.c countspheres_mocks_impl_double.c \
        $(UTILS_DIR)/gridlink_impl_float.c $(UTILS_DIR)/gridlink_impl_double.c \
        $(UTILS_DIR)/gridlink_utils_float.c $(UTILS_DIR)/gridlink_utils_double.c \
        $(UTILS_DIR)/utils.c $(UTILS_DIR)/progressbar.c $(UTILS_DIR)/cpu_features.c $(UTILS_DIR)/avx512_calls.c \
	$(UTILS_DIR)/set_cosmo_dist.c $(UTILS_DIR)/cosmology_params.c
LIBRARY_HEADERS := $(LIBNAME).h

TARGET := vpf_mocks
TARGETSRC := $(TARGET).c $(IO_DIR)/io.c $(IO_DIR)/ftread.c $(LIBSRC)
INCL   := $(LIBNAME).h vpf_mocks_kernels.c.src countspheres_mocks_impl.c.src countspheres_mocks_impl.h.src \
        countspheres_mocks_impl_float.h countspheres_mocks_impl_double.h \
        $(UTILS_DIR)/gridlink_impl_float.h $(UTILS_DIR)/gridlink_impl_double.h $(UTILS_DIR)/gridlink_impl.h.src \
        $(UTILS_DIR)/gridlink_utils_float.h $(UTILS_DIR)/gridlink_utils_double.h $(UTILS_DIR)/gridlink_utils.h.src \
        $(UTILS_DIR)/cellarray_float.h $(UTILS_DIR)/cellarray_double.h $(UTILS_DIR)/cellarray.h.src \
        $(IO_DIR)/ftread.h $(IO_DIR)/io.h $(UTILS_DIR)/utils.h $(UTILS_DIR)/function_precision.h \
        $(UTILS_DIR)/avx512_calls.h $(UTILS_DIR)/avx_calls.h $(UTILS_DIR)/sse_calls.h \
	$(UTILS_DIR)/defs.h $(UTILS_DIR)/set_cosmo_dist.h $(UTILS_DIR)/cosmology_params.h $(UTILS_DIR)/progressbar.h $(UTILS_DIR)/cpu_features.h

TARGETOBJS  := $(TARGETSRC:.c=.o)
LIBOBJS :=$(LIBSRC:.c=.o)

EXTRA_INCL:=$(GSL_CFLAGS)
EXTRA_LINK:=$(GSL_LINK)

all: $(TARGET) $(TARGETSRC) Makefile $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk

countspheres_mocks_impl_float.o:countspheres_mocks_impl_float.c countspheres_mocks_impl_float.h vpf_mocks_kernels_float.c $(UTILS_DIR)/gridlink_impl_float.h $(UTILS_DIR)/cellarray_float.h $(UTILS_DIR)/gridlink_utils_float.h
countspheres_mocks_impl_double.o:countspheres_mocks_impl_double.c countspheres_mocks_impl_double.h vpf_mocks_kernels_double.c $(UTILS_DIR)/gridlink_impl_double.h $(UTILS_DIR)/cellarray_double.h $(UTILS_DIR)/gridlink_utils_double.h
countspheres_mocks.o: countspheres_mocks.c countspheres_mocks_impl_float.h countspheres_mocks_impl_double.h

libs: lib
lib: $(LIBRARY)

install: $(INSTALL_BIN_DIR)/$(TARGET) $(INSTALL_LIB_DIR)/$(LIBRARY) $(INSTALL_HEADERS_DIR)/$(LIBRARY_HEADERS)

distclean: clean | $(INSTALL_LIB_DIR) $(INSTALL_HEADERS_DIR) $(INSTALL_BIN_DIR)
	$(RM) $(INSTALL_LIB_DIR)/$(LIBRARY)
	cd $(INSTALL_HEADERS_DIR) && $(RM) $(LIBRARY_HEADERS)
	cd $(INSTALL_BIN_DIR) && $(RM) $(TARGET)

clean:
	$(RM) $(TARGET) $(TARGETOBJS) $(LIBRARY) countspheres_mocks_impl_float.[ch] countspheres_mocks_impl_double.[ch] vpf_mocks_kernels_float.c vpf_mocks_kernels_double.c
	$(RM) -R vpf_mocks.dSYM

tests:
	$(MAKE) -C ../tests vpf_mocks


include $(ROOT_DIR)/rules.mk
