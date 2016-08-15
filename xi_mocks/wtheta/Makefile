ROOT_DIR := ../..
INSTALL_HEADERS_DIR := $(ROOT_DIR)/include
INSTALL_LIB_DIR := $(ROOT_DIR)/lib
INSTALL_BIN_DIR := $(ROOT_DIR)/bin
UTILS_DIR := $(ROOT_DIR)/utils
IO_DIR := $(ROOT_DIR)/io

include $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk

LIBNAME:=countpairs_theta_mocks
LIBRARY := lib$(LIBNAME).a
LIBSRC  := countpairs_theta_mocks.c $(UTILS_DIR)/gridlink_mocks.c $(UTILS_DIR)/utils.c $(UTILS_DIR)/progressbar.c 
LIBRARY_HEADERS := $(LIBNAME).h


TARGET := DDtheta_mocks
TARGETS := $(TARGET) wtheta
TARGETSRC := $(LIBSRC) $(TARGET).c $(IO_DIR)/io.c $(IO_DIR)/ftread.c
INCL     := $(LIBNAME).h $(UTILS_DIR)/defs.h \
            $(IO_DIR)/io.h $(IO_DIR)/ftread.h $(IO_DIR)/io.h  \
            $(UTILS_DIR)/gridlink_mocks.h $(UTILS_DIR)/cellarray_mocks.h \
            $(UTILS_DIR)/progressbar.h $(UTILS_DIR)/utils.h $(UTILS_DIR)/function_precision.h $(UTILS_DIR)/avx_calls.h

TARGETOBJS:=$(TARGETSRC:.c=.o)
LIBOBJS:=$(LIBSRC:.c=.o) 

all: $(TARGETS) $(SRC1) $(INCL) $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk Makefile

wtheta: wtheta.c $(UTILS_DIR)/utils.c 
	$(CC) $(CFLAGS) $(INCLUDE) $^ $(CLINK) -o $@ 

countpairs_theta_mocks.o: countpairs_theta_mocks.c $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

$(UTILS_DIR)/gridlink_mocks.o:$(UTILS_DIR)/gridlink_mocks.c $(ROOT_DIR)/common.mk $(ROOT_DIR)/mocks.options Makefile
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

libs:lib
lib: $(LIBRARY) | $(INSTALL_LIB_DIR)

install: $(INSTALL_BIN_DIR)/DDtheta_mocks $(INSTALL_BIN_DIR)/wtheta $(INSTALL_LIB_DIR)/$(LIBRARY) | $(INSTALL_BIN_DIR) $(INSTALL_LIB_DIR)

tests:
	$(MAKE) -C ../tests wtheta

clean:
	$(RM) $(TARGETS) $(TARGETOBJS) $(LIBRARY)
	$(MAKE) -C $(UTILS_DIR) clean
	$(MAKE) -C $(IO_DIR) clean

distclean:clean | $(INSTALL_LIB_DIR) $(INSTALL_HEADERS_DIR) $(INSTALL_BIN_DIR)
	cd $(INSTALL_LIB_DIR) && $(RM) $(LIBRARY)
	cd $(INSTALL_HEADERS_DIR) && $(RM) $(LIBRARY_HEADERS)
	cd $(INSTALL_BIN_DIR) && $(RM) $(TARGETS)

include $(ROOT_DIR)/rules.mk
