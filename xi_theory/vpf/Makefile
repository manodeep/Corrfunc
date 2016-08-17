ROOT_DIR := ../..
INSTALL_HEADERS_DIR := $(ROOT_DIR)/include
INSTALL_LIB_DIR := $(ROOT_DIR)/lib
INSTALL_BIN_DIR := $(ROOT_DIR)/bin
UTILS_DIR := $(ROOT_DIR)/utils
IO_DIR := $(ROOT_DIR)/io

include $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk

LIBRARY := libcountspheres.a
LIBRARY_HEADERS := countspheres.h
LIBSRC := countspheres.c countspheres_impl_double.c countspheres_impl_float.c \
          $(UTILS_DIR)/gridlink_impl_double.c $(UTILS_DIR)/gridlink_impl_float.c \
          $(UTILS_DIR)/utils.c $(UTILS_DIR)/progressbar.c \
          $(UTILS_DIR)/cpu_features.c 

GSL_TARGET := vpf
TARGETSRC := $(GSL_TARGET).c $(IO_DIR)/ftread.c $(IO_DIR)/io.c $(LIBSRC)
INCL   := vpf_kernels_float.c vpf_kernels_double.c vpf_kernels.c.src countspheres.h \
          countspheres_impl_float.h countspheres_impl_double.h countspheres_impl.h.src countspheres_impl.c.src \
          $(UTILS_DIR)/gridlink_impl_float.h $(UTILS_DIR)/gridlink_impl_double.h $(UTILS_DIR)/gridlink_impl.h.src \
          $(UTILS_DIR)/cellarray_double.h $(UTILS_DIR)/cellarray_float.h $(UTILS_DIR)/cellarray.h.src \
          $(IO_DIR)/ftread.h $(IO_DIR)/io.h $(UTILS_DIR)/utils.h $(UTILS_DIR)/avx_calls.h $(UTILS_DIR)/sse_calls.h \
	  $(UTILS_DIR)/function_precision.h $(UTILS_DIR)/defs.h $(UTILS_DIR)/sglib.h $(UTILS_DIR)/progressbar.h \
          $(UTILS_DIR)/cpu_features.h

TARGETOBJS  := $(TARGETSRC:.c=.o)
LIBOBJS := $(LIBSRC:.c=.o)

all: $(GSL_TARGET) $(TARGETSRC) $(INCL) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile 
libs: lib
lib: $(INSTALL_LIB_DIR)/$(LIBRARY) | $(INSTALL_LIB_DIR) $(INSTALL_HEADERS_DIR)

tests: ../tests/tests_periodic ../tests/tests_nonperiodic 
	$(MAKE) -C ../tests vpf

../tests/tests_periodic:
	$(MAKE) -C ../tests

../tests/tests_nonperiodic:
	$(MAKE) -C ../tests

countspheres_%.o: countspheres_%.c vpf_kernels_%.c $(INCL) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile
	$(CC) $(CFLAGS) $(INCLUDE) $(GSL_CFLAGS) -c $< -o $@

$(GSL_TARGET): $(TARGETOBJS) $(INCL) $(ROOT_DIR)/common.mk $(ROOT_DIR)/theory.options Makefile 
	$(CC) $(TARGETOBJS) $(CLINK) $(GSL_LINK) -o $@

$(GSL_TARGET).o: $(GSL_TARGET).c $(ROOT_DIR)/common.mk Makefile $(INCL)
	$(CC) $(OPT) $(GSL_CFLAGS) $(CFLAGS) $(INCLUDE) -c $< -o $@

install: $(INSTALL_BIN_DIR)/$(GSL_TARGET) $(INSTALL_LIB_DIR)/$(LIBRARY)

clean:
	$(RM) $(TARGETOBJS) $(LIBRARY) $(GSL_TARGET) vpf_kernels_float.[ch] vpf_kernels_double.[ch] countspheres_impl_double.[ch] countspheres_impl_float.[ch]
	$(MAKE) -C $(UTILS_DIR) clean
	$(MAKE) -C $(IO_DIR) clean

distclean:clean
	cd $(INSTALL_HEADERS_DIR) && $(RM) $(LIBRARY_HEADERS)
	cd $(INSTALL_LIB_DIR) && $(RM) $(LIBRARY)

include $(ROOT_DIR)/rules.mk
