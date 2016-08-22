TARGETOBJS := $(TARGETSRC:.c=.o)
LIBOBJS :=$(LIBSRC:.c=.o)

.SUFFIXES:

$(TARGET): $(TARGETOBJS) $(ROOT_DIR)/common.mk 
	$(CC) $(TARGETOBJS) $(C_LIBRARIES) $(CLINK) $(EXTRA_LINK) -o $@

$(TARGET).o: $(TARGET).c $(ROOT_DIR)/common.mk Makefile $(INCL) $(ROOT_DIR)/theory.options $(ROOT_DIR)/mocks.options 
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) $(EXTRA_INCL) -c $< -o $@

$(UTILS_DIR)/%.o: $(ROOT_DIR)/common.mk Makefile $(UTILS_DIR)/%.c
	$(MAKE) -C $(UTILS_DIR)

$(UTILS_DIR)/%.h: $(UTILS_DIR)/%.h.src
	$(MAKE) -C $(UTILS_DIR)

$(IO_DIR)/%.h: $(IO_DIR)/%.h.src
	$(MAKE) -C $(IO_DIR)

$(IO_DIR)/%.o: $(IO_DIR)/%.c $(ROOT_DIR)/common.mk Makefile
	$(MAKE) -C $(IO_DIR)

%_float.c: %.c.src Makefile
	@echo "/* This file is auto-generated from $(notdir $<) */" > $@
	@echo "#ifdef DOUBLE_PREC" >> $@
	@echo "#undef DOUBLE_PREC" >> $@
	@echo "#endif" >> $@
	sed -e "s/DOUBLE/float/g" $< >> $@

%_float.h: %.h.src Makefile
	@echo "/* This file is auto-generated from $(notdir $<) */" > $@
	@echo "#ifdef DOUBLE_PREC" >> $@
	@echo "#undef DOUBLE_PREC" >> $@
	@echo "#endif" >> $@
	sed -e "s/DOUBLE/float/g"  $< >> $@

%_double.c: %.c.src Makefile
	@echo "/* This file is auto-generated from $(notdir $<) */" > $@
	@echo "#ifndef DOUBLE_PREC" >> $@
	@echo "#define DOUBLE_PREC" >> $@
	@echo "#endif" >> $@
	sed -e "s/DOUBLE/double/g"  $< >> $@

%_double.h: %.h.src Makefile
	@echo "/* This file is auto-generated from $(notdir $<) */" > $@
	@echo "#ifndef DOUBLE_PREC" >> $@
	@echo "#define DOUBLE_PREC" >> $@
	@echo "#endif" >> $@
	sed -e "s/DOUBLE/double/g"  $< >> $@

%_double.o: %_double.c %.c.src Makefile
	$(CC) -DDOUBLE_PREC $(CFLAGS) $(INCLUDE) $(EXTRA_INCL) -c $< -o $@

%_float.o: %_float.c %.c.src Makefile
	$(CC) -DNDOUBLE_PREC $(CFLAGS) $(INCLUDE) $(EXTRA_INCL) -c $< -o $@

%.o: %.c $(INCL) $(ROOT_DIR)/common.mk Makefile
	$(CC) $(CFLAGS) $(INCLUDE) $(EXTRA_INCL) -c $< -o $@

$(LIBRARY): $(LIBOBJS) $(ROOT_DIR)/mocks.options $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile 
	ar rcs $@ $(LIBOBJS)


$(INSTALL_LIB_DIR)/%.a: %.a $(UTILS_DIR)/defs.h | $(INSTALL_LIB_DIR) $(INSTALL_HEADERS_DIR)
	cp -p $(LIBRARY) $(INSTALL_LIB_DIR)/
	cp -p $(LIBRARY_HEADERS) $(INSTALL_HEADERS_DIR)/
	cp -p $(UTILS_DIR)/defs.h $(INSTALL_HEADERS_DIR)/

$(INSTALL_BIN_DIR)/%: %
	cp -p $< $(INSTALL_BIN_DIR)/

$(INSTALL_LIB_DIR) $(INSTALL_BIN_DIR) $(INSTALL_HEADERS_DIR):
	mkdir -p $@

.PHONY: clean clena celan install lib tests distclean all realclean libs

celna: clean
celan: clean
clena: clean

