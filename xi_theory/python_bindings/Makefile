ROOT_DIR := ../..
INSTALL_HEADERS_DIR := $(ROOT_DIR)/include
INSTALL_LIB_DIR := $(ROOT_DIR)/lib
INSTALL_BIN_DIR := $(ROOT_DIR)/bin
UTILS_DIR := $(ROOT_DIR)/utils
IO_DIR := $(ROOT_DIR)/io

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
LIBRARY := $(PROJECT).so.$(MAJOR).$(MINOR).$(PATCHLEVEL)
SOURCES := $(PROJECT).c 
OBJECTS := $(SOURCES:.c=.o)
C_LIBRARY := $(DD_DIR)/lib$(DD_LIB).a $(DDrppi_DIR)/lib$(DDrppi_LIB).a $(WP_DIR)/lib$(WP_LIB).a \
             $(XI_DIR)/lib$(XI_LIB).a $(VPF_DIR)/lib$(VPF_LIB).a
PHONY_FIX_PY_LINK := 

ifeq ($(FIX_PYTHON_LINK), 1)
  PHONY_FIX_PY_LINK := sharedlib
endif

ifeq ($(UNAME), Darwin)
  LINKER_FLAGS := -install_name "$(LIBRARY)"
else
  LINKER_FLAGS := -Wl,-soname=$(PROJECT).$(MAJOR)
endif



all: $(PHONY_FIX_PY_LINK) $(LIBRARY) $(SOURCES) $(C_LIBRARY) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile | $(ROOT_DIR)/lib

$(LIBRARY): $(SOURCES) $(OBJECTS) $(C_LIBRARY) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile $(C_LIBRARY) $(OBJECTS) | $(ROOT_DIR)/lib
	$(CC) $(OBJECTS) $(C_LIBRARY) $(LINKER_FLAGS) $(CLINK) $(GSL_LINK) $(PYTHON_LINK) -shared -o $@
	$(RM) $(PROJECT).so
	ln -s $(LIBRARY) $(PROJECT).so

sharedlib: $(LIBRARY) | $(ROOT_DIR)/lib
	@{ \
		CURRENT_PYTHON_LIB=`otool -L $(LIBRARY) | grep -i python | cut -d " " -f1 | xargs` ; \
		PYTHON_LIB_NAME=$(PYTHON_LIB_BASE).dylib ; \
		LINK_PYTHON_LIB=$(PYTHON_LIBDIR)/$$PYTHON_LIB_NAME ;\
		if [[ "$$CURRENT_PYTHON_LIB" != "$$LINK_PYTHON_LIB" ]] ; then \
		install_name_tool -change $$CURRENT_PYTHON_LIB $$LINK_PYTHON_LIB $(LIBRARY); \
		fi ;\
	}


$(ROOT_DIR)/lib:
	mkdir -p $(ROOT_DIR)/lib

$(C_LIBRARY):
	$(MAKE) -C ../ libs

%.o: %.c $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile $(C_LIBRARY)
	$(CC) $(OPT) $(CFLAGS) -I$(ROOT_DIR)/include $(PYTHON_CFLAGS)  $(INCLUDE) $< -c -o $@

install: $(PHONY_FIX_PY_LINK) $(ROOT_DIR)/lib/$(LIBRARY) | $(ROOT_DIR)/lib

$(ROOT_DIR)/lib/$(LIBRARY): $(LIBRARY)
	cp -p $(LIBRARY) $(ROOT_DIR)/lib/
	$(RM) $(ROOT_DIR)/lib/$(PROJECT).so
	ln -s $(LIBRARY) $(ROOT_DIR)/lib/$(PROJECT).so


tests: $(PHONY_FIX_PY_LINK) $(LIBRARY) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile $(C_LIBRARY) $(OBJECTS) | $(ROOT_DIR)/lib
	python call_correlation_functions.py

distclean: clean | $(ROOT_DIR)/lib
	$(RM) $(ROOT_DIR)/lib/$(LIBRARY)
	$(RM) $(ROOT_DIR)/lib/$(PROJECT).so

.PHONY: clean celna clena celan install distclean sharedlib 

clean:
	$(RM) $(LIBRARY) $(OBJECTS) $(PROJECT).so

celan: clean
celna: clean
clena: clean

