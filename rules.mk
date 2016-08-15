TARGETOBJS := $(TARGETSRC:.c=.o)
LIBOBJS :=$(LIBSRC:.c=.o)

$(TARGET): $(TARGETOBJS) $(ROOT_DIR)/common.mk $(ROOT_DIR)/theory.options Makefile 
	$(CC) $(TARGETOBJS) $(CLINK) -o $@

$(TARGET).o: $(TARGET).c $(ROOT_DIR)/common.mk Makefile $(INCL)
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

$(UTILS_DIR)/%.o: $(ROOT_DIR)/common.mk Makefile $(UTILS_DIR)/%.c
	$(MAKE) -C $(UTILS_DIR)

$(UTILS_DIR)/%.h: $(UTILS_DIR)/%.h.src
	$(MAKE) -C $(UTILS_DIR)

$(IO_DIR)/%.h: $(IO_DIR)/%.h.src
	$(MAKE) -C $(IO_DIR)

$(IO_DIR)/%.o: $(IO_DIR)/%.c $(ROOT_DIR)/common.mk Makefile
	$(MAKE) -C $(IO_DIR)

%_float.c: %.c.src Makefile
	@$(ECHO_COMMAND) "/* This file is auto-generated from $(notdir $<) */" > $@
	@$(ECHO_COMMAND) "#ifdef DOUBLE_PREC" >> $@
	@$(ECHO_COMMAND) "#undef DOUBLE_PREC" >> $@
	@$(ECHO_COMMAND) "#endif" >> $@
	sed -e "s/DOUBLE/float/g" $< >> $@

%_float.h: %.h.src Makefile
	@$(ECHO_COMMAND) "/* This file is auto-generated from $(notdir $<) */" > $@
	@$(ECHO_COMMAND) "#ifdef DOUBLE_PREC" >> $@
	@$(ECHO_COMMAND) "#undef DOUBLE_PREC" >> $@
	@$(ECHO_COMMAND) "#endif" >> $@
	sed -e "s/DOUBLE/float/g"  $< >> $@

%_double.c: %.c.src Makefile
	@$(ECHO_COMMAND) "/* This file is auto-generated from $(notdir $<) */" > $@
	@$(ECHO_COMMAND) "#ifndef DOUBLE_PREC" >> $@
	@$(ECHO_COMMAND) "#define DOUBLE_PREC" >> $@
	@$(ECHO_COMMAND) "#endif" >> $@
	sed -e "s/DOUBLE/double/g"  $< >> $@

%_double.h: %.h.src Makefile
	@$(ECHO_COMMAND) "/* This file is auto-generated from $(notdir $<) */" > $@
	@$(ECHO_COMMAND) "#ifndef DOUBLE_PREC" >> $@
	@$(ECHO_COMMAND) "#define DOUBLE_PREC" >> $@
	@$(ECHO_COMMAND) "#endif" >> $@
	sed -e "s/DOUBLE/double/g"  $< >> $@

%.o: %.c $(INCL) $(ROOT_DIR)/common.mk Makefile
	@echo running default make rule
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

$(LIBRARY): $(LIBOBJS) $(ROOT_DIR)/mocks.options $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile 
	ar rcs $@ $(LIBOBJS)

$(INSTALL_LIB_DIR)/%.a: %.a $(UTILS_DIR)/defs.h | $(INSTALL_LIB_DIR) $(INSTALL_HEADERS_DIR)
	cp -p $(LIBRARY) $(INSTALL_LIB_DIR)/
	sed -e "s/DOUBLE/$(VECTOR_TYPE)/g" $(LIBRARY_HEADERS) > $(INSTALL_HEADERS_DIR)/$(LIBRARY_HEADERS)
	cp -p $(UTILS_DIR)/defs.h $(INSTALL_HEADERS_DIR)/

$(INSTALL_BIN_DIR)/%: %
	cp -p $< $(INSTALL_BIN_DIR)/

$(INSTALL_LIB_DIR) $(INSTALL_BIN_DIR) $(INSTALL_HEADERS_DIR):
	mkdir -p $@

.PHONY: clean clena celan install lib tests distclean all realclean libs

celna: clean
celan: clean
clena: clean

