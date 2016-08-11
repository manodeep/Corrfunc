ROOT_DIR := ../..
INSTALL_HEADERS_DIR := $(ROOT_DIR)/include
INSTALL_LIB_DIR := $(ROOT_DIR)/lib
INSTALL_BIN_DIR := $(ROOT_DIR)/bin
UTILS_DIR := $(ROOT_DIR)/utils
IO_DIR := $(ROOT_DIR)/io

include $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk

LIBRARY := libcountpairs_xi.a
LIBRARY_HEADERS := countpairs_xi.h 
LIBSRC := countpairs_xi.c countpairs_xi_impl_double.c countpairs_xi_impl_float.c  \
          $(UTILS_DIR)/gridlink_impl_double.c $(UTILS_DIR)/gridlink_impl_float.c \
          $(UTILS_DIR)/utils.c $(UTILS_DIR)/progressbar.c $(UTILS_DIR)/cpu_features.c

TARGET := xi
TARGETSRC := $(TARGET).c $(IO_DIR)/ftread.c $(IO_DIR)/io.c $(LIBSRC)

INCL   := xi_kernels_float.c xi_kernels_double.c xi_kernels.c.src countpairs_xi.h \
          countpairs_xi_impl_float.h countpairs_xi_impl_double.h countpairs_xi_impl.h.src countpairs_xi_impl.c.src \
          $(UTILS_DIR)/gridlink_impl_float.h $(UTILS_DIR)/gridlink_impl_double.h $(UTILS_DIR)/gridlink_impl.h.src \
          $(UTILS_DIR)/cellarray_double.h $(UTILS_DIR)/cellarray_float.h $(UTILS_DIR)/cellarray.h.src \
          $(IO_DIR)/ftread.h $(IO_DIR)/io.h $(UTILS_DIR)/utils.h $(UTILS_DIR)/avx_calls.h $(UTILS_DIR)/sse_calls.h \
	  $(UTILS_DIR)/function_precision.h $(UTILS_DIR)/defs.h $(UTILS_DIR)/sglib.h $(UTILS_DIR)/progressbar.h

TARGETOBJS  := $(TARGETSRC:.c=.o)
LIBOBJS := $(LIBSRC:.c=.o)

all: $(TARGET) $(TARGETSRC) $(INCL) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile 

$(TARGET): $(TARGETOBJS) $(INCL) $(ROOT_DIR)/common.mk $(ROOT_DIR)/theory.options Makefile 
	$(CC) $(TARGETOBJS) $(CLINK) -o $@

%_float.c: %.c.src Makefile
	@$(ECHO_COMMAND) "/* This file is auto-generated from $(notdir $<) */" > $@
	sed -e "s/DOUBLE/float/g" $< >> $@

%_float.h: %.h.src Makefile
	@$(ECHO_COMMAND) "/* This file is auto-generated from $(notdir $<) */" > $@
	sed -e "s/DOUBLE/float/g"  $< >> $@

%_double.c: %.c.src Makefile
	@$(ECHO_COMMAND) "/* This file is auto-generated from $(notdir $<) */" > $@
	sed -e "s/DOUBLE/double/g"  $< >> $@
%_double.h: %.h.src Makefile
	@$(ECHO_COMMAND) "/* This file is auto-generated from $(notdir $<) */" > $@
	sed -e "s/DOUBLE/double/g"  $< >> $@

%_double.o: %_double.c %_double.h $(INCL) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile
	$(CC) -DDOUBLE_PREC $(CFLAGS) $(INCLUDE) -c $< -o $@

%_float.o: %_float.c %_float.h $(INCL)
	$(CC) -DNDOUBLE_PREC $(CFLAGS) $(INCLUDE) -c $< -o $@

countpairs_xi.o: countpairs_xi.c countpairs_xi.h countpairs_xi_impl_double.h countpairs_xi_impl_float.h $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

$(TARGET).o: $(TARGET).c $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile $(INCL)
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

%.o: %.c %.h $(INCL) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

libs: lib
lib:  $(INSTALL_LIB_DIR)/$(LIBRARY) | $(INSTALL_LIB_DIR) $(INSTALL_HEADERS_DIR) 

$(LIBRARY): $(LIBOBJS) $(INCL) $(ROOT_DIR)/common.mk $(ROOT_DIR)/theory.options Makefile | $(INSTALL_LIB_DIR) $(INSTALL_HEADERS_DIR)
	ar rcs $@ $(LIBOBJS)

$(INSTALL_LIB_DIR)/$(LIBRARY): $(LIBRARY)
	cp -p $< $@
	cp -p $(LIBRARY_HEADERS) $(INSTALL_HEADERS_DIR)/

.PHONY: clean clena celan install lib tests distclean

$(INSTALL_LIB_DIR) $(INSTALL_BIN_DIR) $(INSTALL_HEADERS_DIR):
	mkdir -p $@

install: $(INSTALL_BIN_DIR)/$(TARGET) | $(INSTALL_BIN_DIR)

$(INSTALL_BIN_DIR)/$(TARGET): $(TARGET)
	cp -p $< $@

distclean: clean | $(INSTALL_LIB_DIR) $(INSTALL_HEADERS_DIR) $(INSTALL_BIN_DIR)
	$(RM) $(INSTALL_LIB_DIR)/$(LIBRARY) $(INSTALL_HEADERS_DIR)/$(LIBRARY_HEADERS)
	cd $(INSTALL_BIN_DIR) && $(RM) $(TARGET)

clean:
	$(RM) $(TARGET) $(TARGETOBJS) $(LIBOBJS) $(LIBRARY) *_float.[ch] *_double.[ch] $(UTILS_DIR)/*_double.[ch] $(UTILS_DIR)/*_float.[ch]

clena: clean
celan: clean

tests: tests_periodic 

tests_periodic: 
	$(MAKE) -C ../tests xi

