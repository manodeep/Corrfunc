ROOT_DIR := ../..
INSTALL_HEADERS_DIR := $(ROOT_DIR)/include
INSTALL_LIB_DIR := $(ROOT_DIR)/lib
INSTALL_BIN_DIR := $(ROOT_DIR)/bin
UTILS_DIR := $(ROOT_DIR)/utils
IO_DIR := $(ROOT_DIR)/io

include $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk

TARGET := DDrppi
TARGETS := $(TARGET) wprp
LIBRARY := libcountpairs_rp_pi.a
LIBSRC  := countpairs_rp_pi.c countpairs_rp_pi_impl_double.c countpairs_rp_pi_impl_float.c \
         $(UTILS_DIR)/gridlink_impl_double.c $(UTILS_DIR)/gridlink_impl_float.c \
         $(UTILS_DIR)/utils.c $(UTILS_DIR)/progressbar.c $(UTILS_DIR)/cpu_features.c
LIBRARY_HEADERS := countpairs_rp_pi.h

TARGETSRC := DDrppi.c $(IO_DIR)/ftread.c $(IO_DIR)/io.c $(LIBSRC)
WPRPSRC := wprp.c $(UTILS_DIR)/utils.c

INCL   := countpairs_rp_pi_kernels_float.c countpairs_rp_pi_kernels_double.c countpairs_rp_pi_kernels.c.src countpairs_rp_pi_impl.c.src countpairs_rp_pi_impl.h.src \
          countpairs_rp_pi.h countpairs_rp_pi_impl_double.h countpairs_rp_pi_impl_float.h \
          $(UTILS_DIR)/gridlink_impl_float.h $(UTILS_DIR)/gridlink_impl_double.h $(UTILS_DIR)/gridlink_impl.h.src \
          $(UTILS_DIR)/cellarray_float.h $(UTILS_DIR)/cellarray_double.h $(UTILS_DIR)/cellarray.h.src \
	  $(UTILS_DIR)/function_precision.h  $(UTILS_DIR)/avx_calls.h $(UTILS_DIR)/sse_calls.h \
	  $(UTILS_DIR)/defs.h $(UTILS_DIR)/cpu_features.h \
          $(IO_DIR)/ftread.h $(IO_DIR)/io.h $(UTILS_DIR)/utils.h $(UTILS_DIR)/progressbar.h 

TARGETOBJS  := $(TARGETSRC:.c=.o)
LIBOBJS  := $(LIBSRC:.c=.o)

all: $(TARGETS) $(TARGETSRC) $(WPRPSRC) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile 

$(TARGET): $(TARGETOBJS) $(INCL) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile
	$(CC) $(TARGETOBJS) $(CLINK) -o $@

wprp: $(WPRPOBJS) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $(WPRPSRC) $(CLINK)

$(UTILS_DIR)/%.o: $(ROOT_DIR)/common.mk Makefile $(UTILS_DIR)/%.c
	$(MAKE) -C $(UTILS_DIR)

$(UTILS_DIR)/%.h: $(ROOT_DIR)/common.mk Makefile $(UTILS_DIR)/%.h.src
	$(MAKE) -C $(UTILS_DIR)

$(IO_DIR)/%.h: $(ROOT_DIR)/common.mk Makefile $(IO_DIR)/%.h.src
	$(MAKE) -C $(IO_DIR)

$(TARGET).o: $(TARGET).c $(INCL) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

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

%.o: %.c %.h $(INCL) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

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

$(INSTALL_BIN_DIR)/%: % | $(INSTALL_BIN_DIR)
	cp -p $< $(INSTALL_BIN_DIR)/

install: $(INSTALL_BIN_DIR)/DDrppi $(INSTALL_BIN_DIR)/wprp | $(INSTALL_BIN_DIR)

.PHONY: clean clena celan install lib tests distclean

clean:
	$(RM) $(TARGETS) $(LIBOBJS) $(WPRPOBJS) $(LIBRARY) *_double.[ch] *_float.[ch] $(UTILS_DIR)/*_double.[ch] $(UTILS_DIR)/*_float.[ch]
	$(RM) -R *.dSYM

clena: clean
celan: celan

distclean:clean | $(INSTALL_LIB_DIR) $(INSTALL_HEADERS_DIR) $(INSTALL_BIN_DIR)
	$(RM) $(INSTALL_LIB_DIR)/$(LIBRARY)
	cd $(INSTALL_HEADERS_DIR) && $(RM) $(LIBRARY_HEADERS)
	cd $(INSTALL_BIN_DIR) && $(RM) $(TARGETS)

tests: 
	$(MAKE) -C ../tests xi_rp_pi
