TARGETOBJS := $(TARGETSRC:.c=.o)
LIBOBJS :=$(LIBSRC:.c=.o)

$(TARGET): $(TARGETOBJS) $(INCL) $(ROOT_DIR)/common.mk $(ROOT_DIR)/theory.options Makefile 
	$(CC) $(TARGETOBJS) $(CLINK) $(GSL_LINK) -o $@

$(TARGET).o: $(TARGET).c $(ROOT_DIR)/common.mk Makefile $(INCL)
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) $(GSL_CFLAGS) -c $< -o $@

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
	$(CC) -DDOUBLE_PREC $(CFLAGS) $(INCLUDE) -c $< -o $@

%_float.o: %_float.c %_float.h $(INCL)
	$(CC) -DNDOUBLE_PREC $(CFLAGS) $(INCLUDE) -c $< -o $@

%.o: %.c %.h $(INCL) $(ROOT_DIR)/common.mk Makefile
	$(CC) $(CFLAGS) $(INCLUDE) $(GSL_CFLAGS) -c $< -o $@


$(LIBRARY): $(LIBOBJS) $(INCL) $(ROOT_DIR)/mocks.options $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile 
	ar rcs $@ $(LIBOBJS)

$(INSTALL_LIB_DIR)/%.a: %.a $(UTILS_DIR)/defs.h | $(INSTALL_LIB_DIR) $(INSTALL_HEADERS_DIR)
	cp -p $(LIBRARY) $(INSTALL_LIB_DIR)/
	sed -e "s/DOUBLE/$(VECTOR_TYPE)/g" $(LIBRARY_HEADERS) > $(INSTALL_HEADERS_DIR)/$(LIBRARY_HEADERS)
	cp -p $(UTILS_DIR)/defs.h $(INSTALL_HEADERS_DIR)/


$(INSTALL_BIN_DIR)/%: %
	cp -p $< $(INSTALL_BIN_DIR)/
$(INSTALL_LIB_DIR) $(INSTALL_BIN_DIR) $(INSTALL_HEADERS_DIR):
	mkdir -p $@

.PHONY: clean clena celan install lib tests distclean all realclean

celna: clean
celan: clean
clena: clean

