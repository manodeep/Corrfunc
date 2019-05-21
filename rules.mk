TARGETOBJS := $(TARGETSRC:.c=.o)
LIBOBJS :=$(LIBSRC:.c=.o)

gridlink_impl_double.o:gridlink_impl_double.c gridlink_impl_double.h gridlink_utils_double.h cell_pair_double.h cellarray_double.h gridlink_utils_double.c
gridlink_impl_float.o:gridlink_impl_float.c gridlink_impl_float.h gridlink_utils_float.h cell_pair_float.h cellarray_float.h gridlink_utils_float.c
gridlink_mocks_impl_double.o:gridlink_mocks_impl_double.c gridlink_mocks_impl_double.h gridlink_utils_double.h cell_pair_double.h cellarray_mocks_double.h
gridlink_mocks_impl_float.o:gridlink_mocks_impl_float.c gridlink_mocks_impl_float.h gridlink_utils_float.h cell_pair_float.h cellarray_mocks_float.h

cellarray_double.h:weight_functions_double.h
cellarray_float.h:weight_functions_float.h
weight_functions_double.h:weight_defs_double.h
weight_functions_float.h:weight_defs_float.h
gridlink_mocks_impl_double.h:cellarray_mocks_double.h
gridlink_mocks_impl_float.h:cellarray_mocks_float.h
gridlink_impl_double.h:cellarray_double.h
gridlink_impl_float.h:cellarray_float.h

.SUFFIXES:

$(TARGET): $(TARGETOBJS) $(ROOT_DIR)/common.mk
	$(CC) $(TARGETOBJS) $(C_LIBRARIES) $(CLINK) $(EXTRA_LINK) -o $@

$(TARGET).o: $(TARGET).c $(ROOT_DIR)/common.mk Makefile $(ROOT_DIR)/theory.options $(ROOT_DIR)/mocks.options
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) $(EXTRA_INCL) -c $< -o $@

%_float.c: %.c.src Makefile
	@echo "/* This file is auto-generated from $(notdir $<) */" > $@
	@echo "#ifdef DOUBLE_PREC" >> $@
	@echo "#undef DOUBLE_PREC" >> $@
	@echo "#endif" >> $@
	sed -e "/DOUBLE_PREC/!s/DOUBLE/float/g" $< >> $@

%_float.h: %.h.src Makefile
	@echo "/* This file is auto-generated from $(notdir $<) */" > $@
	@echo "#ifdef DOUBLE_PREC" >> $@
	@echo "#undef DOUBLE_PREC" >> $@
	@echo "#endif" >> $@
	sed -e "/DOUBLE_PREC/!s/DOUBLE/float/g"  $< >> $@

%_double.c: %.c.src Makefile
	@echo "/* This file is auto-generated from $(notdir $<) */" > $@
	@echo "#ifndef DOUBLE_PREC" >> $@
	@echo "#define DOUBLE_PREC" >> $@
	@echo "#endif" >> $@
	sed -e "/DOUBLE_PREC/!s/DOUBLE/double/g"  $< >> $@

%_double.h: %.h.src Makefile
	@echo "/* This file is auto-generated from $(notdir $<) */" > $@
	@echo "#ifndef DOUBLE_PREC" >> $@
	@echo "#define DOUBLE_PREC" >> $@
	@echo "#endif" >> $@
	sed -e "/DOUBLE_PREC/!s/DOUBLE/double/g"  $< >> $@

$(LIBNAME).o:$(LIBNAME).c $(LIBNAME).h $(LIBNAME)_impl_double.h $(LIBNAME)_impl_float.h
$(LIBNAME)_impl_double.h: $(UTILS_DIR)/weight_defs_double.h
$(LIBNAME)_impl_float.h:  $(UTILS_DIR)/weight_defs_float.h
$(LIBNAME)_impl_double.o: $(UTILS_DIR)/weight_functions_double.h
$(LIBNAME)_impl_float.o: $(UTILS_DIR)/weight_functions_float.h

%_double.o: %_double.c %_double.h
	$(CC) -DDOUBLE_PREC $(CFLAGS) $(INCLUDE) $(EXTRA_INCL) -c $< -o $@

%_float.o: %_float.c %_float.h
	$(CC) -DNDOUBLE_PREC $(CFLAGS) $(INCLUDE) $(EXTRA_INCL) -c $< -o $@

%.o: %.c $(ROOT_DIR)/common.mk $(ROOT_DIR)/utils/defs.h Makefile
	$(CC) $(CFLAGS) $(INCLUDE) $(EXTRA_INCL) -c $< -o $@

$(LIBRARY): $(LIBOBJS) $(ROOT_DIR)/mocks.options $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile
	ar rcs $@ $(LIBOBJS)

$(INSTALL_LIB_DIR)/%.a: %.a | $(INSTALL_LIB_DIR)
	cp -p $(LIBRARY) $(INSTALL_LIB_DIR)/

$(INSTALL_HEADERS_DIR)/%.h: %.h $(INSTALL_HEADERS_DIR)/defs.h | $(INSTALL_HEADERS_DIR)
	cp -p $< $@

$(INSTALL_HEADERS_DIR)/defs.h:$(UTILS_DIR)/defs.h | $(INSTALL_HEADERS_DIR)
	cp -p $(UTILS_DIR)/defs.h $(INSTALL_HEADERS_DIR)/

$(INSTALL_BIN_DIR)/%: %
	cp -p $< $(INSTALL_BIN_DIR)/

$(INSTALL_LIB_DIR) $(INSTALL_BIN_DIR) $(INSTALL_HEADERS_DIR):
	mkdir -p $@

ifdef PROJECT
  ifeq ($(UNAME), Darwin)
    PYTHON_LINK += -Xlinker -install_name -Xlinker "$(PROJECT)$(PYTHON_SOABI).so"
  else
    PYTHON_LINK += -Xlinker -soname -Xlinker "$(PROJECT)$(PYTHON_SOABI).so.$(MAJOR)"
  endif

endif



.PHONY: clean clena celan install lib tests distclean all realclean libs

celna: clean
celan: clean
clena: clean
