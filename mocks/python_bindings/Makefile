ROOT_DIR := ../..
UTILS_DIR := $(ROOT_DIR)/utils
IO_DIR := $(ROOT_DIR)/io
LIB_DIR := $(ROOT_DIR)/lib
DATA_DIR := ../tests/data
INSTALL_HEADERS_DIR := $(ROOT_DIR)/include
INSTALL_LIB_DIR := $(ROOT_DIR)/lib
INSTALL_BIN_DIR := $(ROOT_DIR)/bin

MOCKS_DIR := $(ROOT_DIR)/mocks
DDrppi_mocks_DIR := $(MOCKS_DIR)/DDrppi
WTHETA_DIR := $(MOCKS_DIR)/wtheta
VPF_mocks_DIR := $(MOCKS_DIR)/vpf

DDrppi_mocks_LIB := countpairs_rp_pi_mocks
WTHETA_LIB := countpairs_theta_mocks
VPF_mocks_LIB := countspheres_mocks

include $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk

PROJECT := _countpairs_mocks
PYTHON_EXTN := $(PROJECT).so.$(MAJOR).$(MINOR).$(PATCHLEVEL)
SOURCES := $(PROJECT).c
OBJECTS := $(SOURCES:.c=.o)
C_LIBRARIES := $(DDrppi_mocks_DIR)/lib$(DDrppi_mocks_LIB).a $(WTHETA_DIR)/lib$(WTHETA_LIB).a $(VPF_mocks_DIR)/lib$(VPF_mocks_LIB).a
INCL :=	$(DDrppi_mocks_DIR)/$(DDrppi_mocks_LIB).h $(WTHETA_DIR)/$(WTHETA_LIB).h $(VPF_mocks_DIR)/$(VPF_mocks_LIB).h \
        $(UTILS_DIR)/defs.h $(IO_DIR)/io.h $(IO_DIR)/ftread.h \
        $(UTILS_DIR)/utils.h \
	$(UTILS_DIR)/function_precision.h \
        $(UTILS_DIR)/progressbar.h $(UTILS_DIR)/cosmology_params.h \
        $(UTILS_DIR)/cpu_features.h $(UTILS_DIR)/macros.h

LIB_INCLUDE:=-I$(DDrppi_mocks_DIR) -I$(WTHETA_DIR) -I$(VPF_mocks_DIR)


ifeq ($(UNAME), Darwin)
  LINKER_FLAGS := -install_name "$(PYTHON_EXTN)"
else
  LINKER_FLAGS := -Wl,-soname=$(PROJECT).$(MAJOR)
endif

all: sharedlib $(LIBRARY) $(SOURCES) $(C_LIBRARIES) $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk Makefile 

$(PROJECT).so: $(PYTHON_EXTN)
	$(RM) $(PROJECT).so
	ln -s $(PYTHON_EXTN) $(PROJECT).so

$(PROJECT).o: $(PROJECT).c $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk Makefile $(C_LIBRARIES) $(INCL)
	$(CC) $(CFLAGS) $(LIB_INCLUDE) $(PYTHON_CFLAGS)  $(INCLUDE) $< -c -o $@

tests: sharedlib $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk Makefile $(C_LIBRARIES) $(OBJECTS) | $(ROOT_DIR)/lib
	python call_correlation_functions_mocks.py

$(PYTHON_EXTN): $(OBJECTS) $(C_LIBRARIES) $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk Makefile 
	$(CC) $(OBJECTS) $(C_LIBRARIES) $(LINKER_FLAGS) $(CLINK) $(GSL_LINK) $(PYTHON_LINK) -shared -o $@

$(DDrppi_mocks_DIR)/lib$(DDrppi_mocks_LIB).a: $(DDrppi_mocks_DIR)/*.c $(DDrppi_mocks_DIR)/*.c.src $(DDrppi_mocks_DIR)/*.h.src $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk
	$(MAKE) -C $(DDrppi_mocks_DIR) libs

$(WTHETA_DIR)/lib$(WTHETA_LIB).a: $(WTHETA_DIR)/*.c $(WTHETA_DIR)/*.c.src $(WTHETA_DIR)/*.h.src $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk
	$(MAKE) -C $(WTHETA_DIR) libs

$(VPF_mocks_DIR)/lib$(VPF_mocks_LIB).a: $(VPF_mocks_DIR)/*.c $(VPF_mocks_DIR)/*.c.src $(VPF_mocks_DIR)/*.h.src $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk
	$(MAKE) -C $(VPF_mocks_DIR) libs

sharedlib: $(PROJECT).so
ifeq ($(FIX_PYTHON_LINK), 1)
	@{ \
		CURRENT_PYTHON_LIB=`otool -L $(PYTHON_EXTN) | grep -i python | cut -d " " -f1 | grep -v $(PROJECT) | xargs` ; \
		PYTHON_LIB_NAME=$(PYTHON_LIB_BASE).dylib ; \
		LINK_PYTHON_LIB=$(PYTHON_LIBDIR)/$$PYTHON_LIB_NAME ;\
		if [[ "$$CURRENT_PYTHON_LIB" != "$$LINK_PYTHON_LIB" ]] ; then \
		install_name_tool -change $$CURRENT_PYTHON_LIB $$LINK_PYTHON_LIB $(PYTHON_EXTN); \
		fi ;\
	}
endif

install: sharedlib $(LIB_DIR)/$(PYTHON_EXTN) 

.PHONY: sharedlib

$(INSTALL_LIB_DIR)/$(PYTHON_EXTN): $(PYTHON_EXTN) | $(INSTALL_LIB_DIR)
	cp -p $(PYTHON_EXTN) $(INSTALL_LIB_DIR)/

$(INSTALL_LIB_DIR)/$(PROJECT).so:$(INSTALL_LIB_DIR)/$(PYTHON_EXTN)
	$(RM) $(INSTALL_LIB_DIR)/$(PROJECT).so
	cd $(INSTALL_LIB_DIR) && ln -s $(PYTHON_EXTN) $(PROJECT).so


distclean: clean | $(LIB_DIR)
	cd $(INSTALL_LIB_DIR) && $(RM) $(PYTHON_EXTN) $(PROJECT).so

clean:
	$(RM) $(PYTHON_EXTN) $(OBJECTS) $(PROJECT).so


