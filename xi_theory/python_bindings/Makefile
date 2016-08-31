ROOT_DIR := ../..
UTILS_DIR:=$(ROOT_DIR)/utils
INSTALL_HEADERS_DIR := $(ROOT_DIR)/include
INSTALL_LIB_DIR := $(ROOT_DIR)/lib
INSTALL_BIN_DIR := $(ROOT_DIR)/bin

THEORY_DIR := $(ROOT_DIR)/xi_theory
DD_DIR := $(THEORY_DIR)/xi_of_r
DDrppi_DIR := $(THEORY_DIR)/xi_rp_pi
WP_DIR := $(THEORY_DIR)/wp
XI_DIR := $(THEORY_DIR)/xi
VPF_DIR := $(THEORY_DIR)/vpf

DD_LIB := countpairs
DDrppi_LIB := countpairs_rp_pi
WP_LIB := countpairs_wp
XI_LIB := countpairs_xi
VPF_LIB := countspheres

include $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk

PROJECT := _countpairs
SOURCES := $(PROJECT).c
OBJECTS := $(SOURCES:.c=.o)
PYTHON_LIBRARY := $(PROJECT).so.$(MAJOR).$(MINOR).$(PATCHLEVEL)
C_LIBRARIES := $(DD_DIR)/lib$(DD_LIB).a $(DDrppi_DIR)/lib$(DDrppi_LIB).a $(WP_DIR)/lib$(WP_LIB).a \
             $(XI_DIR)/lib$(XI_LIB).a $(VPF_DIR)/lib$(VPF_LIB).a
INCL := $(DD_DIR)/$(DD_LIB).h $(DDrppi_DIR)/$(DDrppi_LIB).h $(WP_DIR)/$(WP_LIB).h \
        $(XI_DIR)/$(XI_LIB).h $(VPF_DIR)/$(VPF_LIB).h 
LIB_INCLUDE:=-I$(DD_DIR) -I$(DDrppi_DIR) -I$(WP_DIR) -I$(XI_DIR) -I$(VPF_DIR)


all: sharedlib $(SOURCES) $(C_LIBRARIES) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile 


$(PROJECT).so: $(PYTHON_LIBRARY)
	$(RM) $(PROJECT).so
	ln -s $(PYTHON_LIBRARY) $(PROJECT).so

$(PROJECT).o: $(PROJECT).c $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile $(C_LIBRARIES) $(INCL)
	$(CC) $(CFLAGS) $(LIB_INCLUDE) $(PYTHON_CFLAGS) $(INCLUDE) $< -c -o $@

$(PYTHON_LIBRARY): $(OBJECTS) $(C_LIBRARIES) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile 
	$(CC) $(OBJECTS) $(C_LIBRARIES) $(LINKER_FLAGS) $(CLINK) $(GSL_LINK) $(PYTHON_LINK) -shared -o $@


$(DD_DIR)/lib$(DD_LIB).a: $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk
	$(MAKE) -C $(DD_DIR) libs

$(DDrppi_DIR)/lib$(DDrppi_LIB).a: $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk
	$(MAKE) -C $(DDrppi_DIR) libs

$(WP_DIR)/lib$(WP_LIB).a: $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk
	$(MAKE) -C $(WP_DIR) libs

$(XI_DIR)/lib$(XI_LIB).a: $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk
	$(MAKE) -C $(XI_DIR) libs

$(VPF_DIR)/lib$(VPF_LIB).a: $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk
	$(MAKE) -C $(VPF_DIR) libs


install: sharedlib $(INSTALL_LIB_DIR)/$(PYTHON_LIBRARY) 

sharedlib: $(PROJECT).so
ifeq ($(FIX_PYTHON_LINK), 1)
	@{ \
		CURRENT_PYTHON_LIB=`otool -L $(PYTHON_LIBRARY) | grep -i python | cut -d " " -f1 | xargs` ; \
		PYTHON_LIB_NAME=$(PYTHON_LIB_BASE).dylib ; \
		LINK_PYTHON_LIB=$(PYTHON_LIBDIR)/$$PYTHON_LIB_NAME ;\
		if [[ "$$CURRENT_PYTHON_LIB" != "$$LINK_PYTHON_LIB" ]] ; then \
		install_name_tool -change $$CURRENT_PYTHON_LIB $$LINK_PYTHON_LIB $(PYTHON_LIBRARY); \
		fi ;\
	}
endif

$(INSTALL_LIB_DIR)/$(PYTHON_LIBRARY): $(PYTHON_LIBRARY) | $(INSTALL_LIB_DIR)
	cp -p $(PYTHON_LIBRARY) $(INSTALL_LIB_DIR)/

$(INSTALL_LIB_DIR)/$(PROJECT).so:$(INSTALL_LIB_DIR)/$(PYTHON_LIBRARY)
	$(RM) $(INSTALL_LIB_DIR)/$(PROJECT).so
	cd $(INSTALL_LIB_DIR) && ln -s $(PYTHON_LIBRARY) $(PROJECT).so

tests: sharedlib $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile $(C_LIBRARIES) $(OBJECTS) | $(ROOT_DIR)/lib
	python call_correlation_functions.py

clean:
	$(RM) $(OBJECTS) $(PYTHON_LIBRARY) $(PROJECT).so
	$(RM) -R *.dSYM

distclean:clean
	cd $(INSTALL_LIB_DIR) && $(RM) $(PYTHON_LIBRARY) $(PROJECT).so

include $(ROOT_DIR)/rules.mk

