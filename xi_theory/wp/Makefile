include ../../theory.options ../../common.mk


INSTALL_HEADERS_DIR := ../../include
INSTALL_LIB_DIR := ../../lib
INSTALL_BIN_DIR := ../../bin
UTILS_DIR := ../../utils
IO_DIR := ../../io

TARGET := wp
TARGETSRC := wp.c countpairs_wp.c countpairs_wp_impl_double.c countpairs_wp_impl_float.c \
          $(IO_DIR)/io.c $(IO_DIR)/ftread.c \
          $(UTILS_DIR)/utils.c $(UTILS_DIR)/gridlink_impl_float.c $(UTILS_DIR)/gridlink_impl_double.c \
          $(UTILS_DIR)/progressbar.c $(UTILS_DIR)/cpu_features.c

LIBRARY := libcountpairs_wp.a 
LIBSRC := countpairs_wp.c countpairs_wp_impl_double.c countpairs_wp_impl_float.c \
         $(UTILS_DIR)/gridlink_impl_double.c $(UTILS_DIR)/gridlink_impl_float.c \
         $(UTILS_DIR)/utils.c $(UTILS_DIR)/progressbar.c $(UTILS_DIR)/cpu_features.c
LIBRARY_HEADER := countpairs_wp.h 

INCL   := wp_kernels_double.c wp_kernels_float.c countpairs_wp.h countpairs_wp_impl_float.h countpairs_wp_impl_double.h \
          $(IO_DIR)/ftread.h $(IO_DIR)/io.h $(UTILS_DIR)/utils.h \
          $(UTILS_DIR)/gridlink_impl_float.h $(UTILS_DIR)/gridlink_impl_double.h \
          $(UTILS_DIR)/function_precision.h $(UTILS_DIR)/cellarray_double.h $(UTILS_DIR)/cellarray_float.h $(UTILS_DIR)/avx_calls.h \
          $(UTILS_DIR)/defs.h $(UTILS_DIR)/sglib.h $(UTILS_DIR)/progressbar.h $(UTILS_DIR)/cpu_features.h

TARGETOBJS  := $(TARGETSRC:.c=.o)
LIBOBJS :=$(LIBSRC:.c=.o)

all: $(TARGET) $(TARGETSRC) $(INCL) ../../theory.options ../../common.mk Makefile 

$(TARGET): $(TARGETOBJS) $(INCL) ../../theory.options ../../common.mk Makefile 
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

%_double.o: %_double.c $(INCL) ../../theory.options ../../common.mk Makefile
	$(CC) -DDOUBLE_PREC $(CFLAGS) $(INCLUDE) -c $< -o $@

%_float.o: %_float.c %_float.h $(INCL)
	$(CC) -DNDOUBLE_PREC $(CFLAGS) $(INCLUDE) -c $< -o $@

countpairs_wp.o: countpairs_wp.c countpairs_wp.h countpairs_wp_impl_double.h countpairs_wp_impl_float.h ../../theory.options ../../common.mk Makefile
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

$(TARGET).o: $(TARGET).c ../../theory.options ../../common.mk Makefile $(INCL)
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

%.o: %.c $(INCL) ../../theory.options ../../common.mk Makefile
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@


libs: lib
lib: $(INSTALL_LIB_DIR)/$(LIBRARY) | $(INSTALL_LIB_DIR) $(INSTALL_HEADERS_DIR)

$(LIBRARY): $(LIBOBJS) $(INCL) ../../theory.options ../../common.mk Makefile | $(INSTALL_LIB_DIR) $(INSTALL_HEADERS_DIR)
	ar rcs $@ $(LIBOBJS)

$(INSTALL_LIB_DIR)/$(LIBRARY): $(LIBRARY)
	cp -p $< $@
	cp -p $(LIBRARY_HEADER) $(INSTALL_HEADERS_DIR)/

.PHONY: clean clena celan install lib tests distclean all

$(INSTALL_LIB_DIR) $(INSTALL_BIN_DIR) $(INSTALL_HEADERS_DIR):
	mkdir -p $@

install: $(INSTALL_BIN_DIR)/$(TARGET) | $(INSTALL_BIN_DIR)

$(INSTALL_BIN_DIR)/$(TARGET): $(TARGET)
	cp -p $< $@

distclean: clean | $(INSTALL_LIB_DIR) $(INSTALL_HEADERS_DIR) $(INSTALL_BIN_DIR)
	$(RM) $(INSTALL_LIB_DIR)/$(LIBRARY) $(INSTALL_HEADERS_DIR)/$(LIBRARY_HEADER)
	cd $(INSTALL_BIN_DIR) && $(RM) $(TARGET)

clean:
	$(RM) $(TARGET) $(TARGETOBJS) $(LIBOBJS) $(LIBRARY) *_float.[ch] *_double.[ch]

clena: clean

celan: celan

tests: tests_periodic 

tests_periodic: 
	$(MAKE) -C ../tests wp

