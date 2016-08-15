ROOT_DIR := ../..
UTILS_DIR := $(ROOT_DIR)/utils
IO_DIR := $(ROOT_DIR)/io
INCLUDE_DIR := $(ROOT_DIR)/include
LIB_DIR := $(ROOT_DIR)/lib
DATA_DIR := ../tests/data

MOCKS_DIR := $(ROOT_DIR)/xi_mocks
DDrppi_DIR := $(MOCKS_DIR)/DDrppi
WTHETA_DIR := $(MOCKS_DIR)/wtheta
VPF_DIR := $(MOCKS_DIR)/vpf

DDrppi_LIB := countpairs_rp_pi_mocks
WTHETA_LIB := countpairs_theta_mocks
VPF_LIB := countspheres_mocks

include $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk

PROJECT := _countpairs_mocks
PYTHON_LIBRARY := $(PROJECT).so.$(MAJOR).$(MINOR).$(PATCHLEVEL)
SOURCES := $(PROJECT).c
OBJECTS := $(SOURCES:.c=.o)
C_LIBRARIES := $(LIB_DIR)/libcountpairs_rp_pi_mocks.a $(LIB_DIR)/libcountpairs_theta_mocks.a $(LIB_DIR)/libcountspheres_mocks.a 
INCL :=	$(INCLUDE_DIR)/$(DDrppi_LIB).h $(INCLUDE_DIR)/$(WTHETA_LIB).h $(INCLUDE_DIR)/$(VPF_LIB).h \
        $(UTILS_DIR)/defs.h $(IO_DIR)/io.h $(IO_DIR)/ftread.h \
        $(UTILS_DIR)/utils.h $(UTILS_DIR)/gridlink_mocks.h \
	$(UTILS_DIR)/function_precision.h $(UTILS_DIR)/cellarray_mocks.h \
        $(UTILS_DIR)/progressbar.h $(UTILS_DIR)/cosmology_params.h 

ifeq ($(UNAME), Darwin)
  LINKER_FLAGS := -install_name "$(PYTHON_LIBRARY)"
else
  LINKER_FLAGS := -Wl,-soname=$(PROJECT).$(MAJOR)
endif

all: sharedlib $(LIBRARY) $(SOURCES) $(C_LIBRARIES) ../../mocks.options ../../common.mk Makefile | $(LIB_DIR)

$(PROJECT).so: $(PYTHON_LIBRARY)
	$(RM) $(PROJECT).so
	ln -s $(PYTHON_LIBRARY) $(PROJECT).so

$(PROJECT).o: $(PROJECT).c $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk Makefile $(C_LIBRARIES)
	$(CC) $(OPT) $(CFLAGS) -I$(ROOT_DIR)/include $(PYTHON_CFLAGS)  $(INCLUDE) $< -c -o $@

tests: sharedlib $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk Makefile $(C_LIBRARIES) $(OBJECTS) | $(ROOT_DIR)/lib
	python call_correlation_functions_mocks.py

$(LIB_DIR)/$(PYTHON_LIBRARY): $(PYTHON_LIBRARY) | $(LIB_DIR)
	cp -p $(PYTHON_LIBRARY) $(LIB_DIR)/

$(LIB_DIR)/$(PROJECT).so:$(LIB_DIR)/$(PYTHON_LIBRARY)
	$(RM) $(LIB_DIR)/$(PROJECT).so
	cd $(LIB_DIR) && ln -s $(PYTHON_LIBRARY) $(PROJECT).so

$(PYTHON_LIBRARY): $(OBJECTS) $(C_LIBRARIES) $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk Makefile 
	$(CC) $(OBJECTS) $(C_LIBRARIES) $(LINKER_FLAGS) $(CLINK) $(GSL_LINK) $(PYTHON_LINK) -shared -o $@

sharedlib: $(PROJECT).so | $(LIB_DIR)
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

install: sharedlib $(LIB_DIR)/$(PYTHON_LIBRARY) 

$(INCLUDE_DIR)/$(DDrppi_LIB).h $(LIB_DIR)/lib$(DDrppi_LIB).a: $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk Makefile
	$(MAKE) -C $(DDrppi_DIR) install

$(INCLUDE_DIR)/$(WTHETA_LIB).h $(LIB_DIR)/lib$(WTHETA_LIB).a: $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk Makefile
	$(MAKE) -C $(WTHETA_DIR) install

$(INCLUDE_DIR)/$(VPF_LIB).h $(LIB_DIR)/lib$(VPF_LIB).a: $(ROOT_DIR)/mocks.options $(ROOT_DIR)/common.mk Makefile
	$(MAKE) -C $(VPF_DIR) install


.PHONY: sharedlib

distclean: clean | $(LIB_DIR)
	cd $(LIB_DIR) && $(RM) $(PYTHON_LIBRARY)
	cd $(LIB_DIR) && $(RM) $(PROJECT).so

clean:
	$(RM) $(PYTHON_LIBRARY) $(OBJECTS) $(PROJECT).so


