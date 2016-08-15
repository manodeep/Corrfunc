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
         $(UTILS_DIR)/utils.c $(UTILS_DIR)/progressbar.c $(UTILS_DIR)/cpu_features.c
LIBRARY_HEADERS := $(LIBNAME).h

TARGET := wp
TARGETSRC := wp.c $(IO_DIR)/io.c $(IO_DIR)/ftread.c $(LIBSRC)

INCL   := wp_kernels_double.c wp_kernels_float.c wp_kernels.c.src countpairs_wp.h \
          countpairs_wp_impl_float.h countpairs_wp_impl_double.h countpairs_wp_impl.h.src \
          $(UTILS_DIR)/gridlink_impl_float.h $(UTILS_DIR)/gridlink_impl_double.h $(UTILS_DIR)/gridlink_impl.h.src \
          $(UTILS_DIR)/cellarray_double.h $(UTILS_DIR)/cellarray_float.h $(UTILS_DIR)/cellarray.h.src \
          $(UTILS_DIR)/avx_calls.h $(UTILS_DIR)/sse_calls.h $(UTILS_DIR)/function_precision.h $(UTILS_DIR)/defs.h $(UTILS_DIR)/cpu_features.h \
          $(IO_DIR)/ftread.h $(IO_DIR)/io.h $(UTILS_DIR)/utils.h $(UTILS_DIR)/sglib.h $(UTILS_DIR)/progressbar.h 

TARGETOBJS  := $(TARGETSRC:.c=.o)
LIBOBJS := $(LIBSRC:.c=.o)

all: $(TARGET) $(TARGETOBJS) $(TARGETSRC) $(INCL) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile 
libs: lib
lib: $(LIBRARY) 
install:$(INSTALL_BIN_DIR)/$(TARGET) $(INSTALL_LIB_DIR)/$(LIBRARY)
clean:
	$(RM) $(TARGETOBJS) $(TARGET) wp_kernels_float.[ch] wp_kernels_float_double.[ch] 
	$(MAKE) -C $(UTILS_DIR) clean
	$(MAKE) -C $(IO_DIR) clean

distclean:clean
	cd $(INSTALL_HEADERS_DIR) && $(RM) $(LIBRARY_HEADERS)
	cd $(INSTALL_LIB_DIR) && $(RM) $(LIBRARY)

tests: tests_periodic 

tests_periodic: 
	$(MAKE) -C ../tests wp

include $(ROOT_DIR)/rules.mk


