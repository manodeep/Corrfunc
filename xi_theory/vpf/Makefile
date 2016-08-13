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

TARGET := vpf
TARGETSRC := $(TARGET).c $(IO_DIR)/ftread.c $(IO_DIR)/io.c $(LIBSRC)
INCL   := vpf_kernels_float.c vpf_kernels_double.c vpf_kernels.c.src countspheres.h \
          countspheres_impl_float.h countspheres_impl_double.h countspheres_impl.h.src countspheres_impl.c.src \
          $(UTILS_DIR)/gridlink_impl_float.h $(UTILS_DIR)/gridlink_impl_double.h $(UTILS_DIR)/gridlink_impl.h.src \
          $(UTILS_DIR)/cellarray_double.h $(UTILS_DIR)/cellarray_float.h $(UTILS_DIR)/cellarray.h.src \
          $(IO_DIR)/ftread.h $(IO_DIR)/io.h $(UTILS_DIR)/utils.h $(UTILS_DIR)/avx_calls.h $(UTILS_DIR)/sse_calls.h \
	  $(UTILS_DIR)/function_precision.h $(UTILS_DIR)/defs.h $(UTILS_DIR)/sglib.h $(UTILS_DIR)/progressbar.h \
          $(UTILS_DIR)/cpu_features.h

TARGETOBJS  := $(TARGETSRC:.c=.o)
LIBOBJS := $(LIBSRC:.c=.o)

all: $(TARGET) $(TARGETSRC) $(INCL) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile 

$(TARGET): $(TARGETOBJS) $(INCL) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile 
	$(CC) $(TARGETOBJS) $(CLINK) $(GSL_LINK) -o $@

$(TARGET).o: $(TARGET).c $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile $(INCL)
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

$(UTILS_DIR)/%.o: $(ROOT_DIR)/common.mk Makefile $(UTILS_DIR)/%.c
	$(MAKE) -C $(UTILS_DIR)

$(UTILS_DIR)/%.h: $(ROOT_DIR)/common.mk Makefile $(UTILS_DIR)/%.h.src
	$(MAKE) -C $(UTILS_DIR)

$(IO_DIR)/%.h: $(ROOT_DIR)/common.mk Makefile $(IO_DIR)/%.h.src
	$(MAKE) -C $(IO_DIR)

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
	$(CC) -DDOUBLE_PREC $(CFLAGS) $(INCLUDE) $(GSL_CFLAGS) -c $< -o $@

%_float.o: %_float.c %_float.h $(INCL)
	$(CC) -DNDOUBLE_PREC $(CFLAGS) $(INCLUDE) $(GSL_CFLAGS) -c $< -o $@

%.o: %.c %.h $(INCL) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile
	$(CC) $(CFLAGS) $(INCLUDE) $(GSL_CFLAGS) -c $< -o $@

libs: lib
lib: $(INSTALL_LIB_DIR)/$(LIBRARY) | $(INSTALL_LIB_DIR) $(INSTALL_HEADERS_DIR)

$(LIBRARY): $(LIBOBJS) $(INCL) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile | $(INSTALL_LIB_DIR) $(INSTALL_HEADERS_DIR)
	ar rcs $@ $(LIBOBJS)

$(INSTALL_LIB_DIR)/$(LIBRARY): $(LIBRARY) $(INSTALL_HEADERS_DIR)/defs.h | $(INSTALL_LIB_DIR)
	cp -p $< $@
	cp -p $(LIBRARY_HEADERS) $(INSTALL_HEADERS_DIR)/

$(INSTALL_HEADERS_DIR)/%.h: $(UTILS_DIR)/%.h | $(INSTALL_HEADERS_DIR)
	cp -p $< $(INSTALL_HEADERS_DIR)/

$(INSTALL_LIB_DIR) $(INSTALL_BIN_DIR) $(INSTALL_HEADERS_DIR):
	mkdir -p $@

install: $(INSTALL_BIN_DIR)/$(TARGET) | $(INSTALL_BIN_DIR)

$(INSTALL_BIN_DIR)/$(TARGET): $(TARGET)
	cp -p $< $@
$(INSTALL_BIN_DIR)/%: % | $(INSTALL_BIN_DIR)
	cp -p $< $(INSTALL_BIN_DIR)/

distclean: clean | $(INSTALL_LIB_DIR) $(INSTALL_HEADERS_DIR) $(INSTALL_BIN_DIR)
	$(RM) $(INSTALL_LIB_DIR)/$(LIBRARY)
	cd $(INSTALL_HEADERS_DIR) && $(RM) $(LIBRARY_HEADERS)
	cd $(INSTALL_BIN_DIR) && $(RM) $(TARGET)

clean:
	$(RM) $(TARGET) $(TARGETOBJS) $(LIBOBJS) $(LIBRARY) $(UTILS_DIR)/*_float.[ch] $(UTILS_DIR)/*_double.[ch] *_float.[ch] *_double.[ch]

clena: clean
celan: clean

.PHONY: clean clena celan install lib tests distclean

tests: ../tests/tests_periodic ../tests/tests_nonperiodic 
	$(MAKE) -C ../tests vpf

../tests/tests_periodic:
	$(MAKE) -C ../tests

../tests/tests_nonperiodic:
	$(MAKE) -C ../tests

