ROOT_DIR:=../..
INSTALL_HEADERS_DIR := $(ROOT_DIR)/include
INSTALL_LIB_DIR := $(ROOT_DIR)/lib
INSTALL_BIN_DIR := $(ROOT_DIR)/bin
UTILS_DIR:=$(ROOT_DIR)/utils
IO_DIR:=$(ROOT_DIR)/io

include $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk 

LIBNAME:=countspheres_mocks
LIBRARY:=lib$(LIBNAME).a
LIBSRC:=countspheres_mocks.c $(UTILS_DIR)/gridlink_mocks.c $(UTILS_DIR)/utils.c $(UTILS_DIR)/progressbar.c \
	$(UTILS_DIR)/set_cosmo_dist.c $(UTILS_DIR)/cosmology_params.c 
LIBRARY_HEADERS := $(LIBNAME).h

GSL_TARGET := vpf_mocks
TARGETSRC := $(GSL_TARGET).c $(IO_DIR)/io.c $(IO_DIR)/ftread.c $(LIBSRC)
INCL   := countspheres_mocks.h $(IO_DIR)/ftread.h $(IO_DIR)/io.h $(UTILS_DIR)/utils.h $(UTILS_DIR)/gridlink_mocks.h \
	$(UTILS_DIR)/function_precision.h $(UTILS_DIR)/cellarray_mocks.h $(UTILS_DIR)/avx_calls.h \
	$(UTILS_DIR)/defs.h $(UTILS_DIR)/set_cosmo_dist.h $(UTILS_DIR)/cosmology_params.h $(UTILS_DIR)/progressbar.h

TARGETOBJS  := $(TARGETSRC:.c=.o)
LIBOBJS :=$(LIBSRC:.c=.o)

all: $(GSL_TARGET) $(TARGETSRC) $(INCL) Makefile $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk

$(GSL_TARGET): $(TARGETOBJS) $(INCL) $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk Makefile 
	$(CC) $(TARGETOBJS) $(CLINK) $(GSL_LINK) -o $@

$(GSL_TARGET).o: $(GSL_TARGET).c $(INCL) $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk Makefile
	$(CC) $(OPT) $(GSL_CFLAGS) $(CFLAGS) $(INCLUDE) -c $< -o $@

$(LIBNAME).o:$(LIBNAME).c $(INCL) $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) $(GSL_CFLAGS) -c $< -o $@

$(UTILS_DIR)/gridlink_mocks.o:$(UTILS_DIR)/gridlink_mocks.c $(ROOT_DIR)/common.mk $(ROOT_DIR)/mocks.options Makefile
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

libs: lib
lib: $(LIBRARY)

install: $(INSTALL_BIN_DIR)/$(GSL_TARGET) $(INSTALL_LIB_DIR)/$(LIBRARY) | $(INSTALL_BIN_DIR) $(INSTALL_LIB_DIR) $(INSTALL_HEADERS_DIR)

distclean: clean | $(INSTALL_LIB_DIR) $(INSTALL_HEADERS_DIR) $(INSTALL_BIN_DIR)
	$(RM) $(INSTALL_LIB_DIR)/$(LIBRARY)
	cd $(INSTALL_HEADERS_DIR) && $(RM) $(LIBRARY_HEADERS)
	cd $(INSTALL_BIN_DIR) && $(RM) $(GSL_TARGET)

clean:
	$(RM) $(GSL_TARGET) $(TARGETOBJS) $(LIBRARY)
	$(MAKE) -C $(UTILS_DIR) clean
	$(MAKE) -C $(IO_DIR) clean

tests: 
	$(MAKE) -C ../tests vpf


include $(ROOT_DIR)/rules.mk
