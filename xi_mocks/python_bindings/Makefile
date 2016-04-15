include ../../mocks.options ../../common.mk

PROJECT := _countpairs_mocks
LIBRARY := $(PROJECT).so.$(MAJOR).$(MINOR).$(PATCHLEVEL)
SOURCES := _countpairs_mocks.c 
OBJECTS := $(SOURCES:.c=.o)
C_LIBRARY := ../../lib/libcountpairs_rp_pi_mocks.a ../../lib/libcountpairs_theta_mocks.a ../../lib/libcountspheres_mocks.a 
PHONY_FIX_PY_LINK := 

ifeq ($(FIX_PYTHON_LINK), 1)
  PHONY_FIX_PY_LINK = sharedlib
endif

ifeq ($(UNAME), Darwin)
  LINKER_FLAGS := -install_name "$(LIBRARY)"
else
  LINKER_FLAGS := -Wl,-soname=$(PROJECT).$(MAJOR)
endif

all: $(PHONY_FIX_PY_LINK) $(LIBRARY) $(SOURCES) $(C_LIBRARY) ../../mocks.options ../../common.mk Makefile | ../../lib

$(LIBRARY): Makefile ../../mocks.options ../../common.mk $(C_LIBRARY) $(OBJECTS) | ../../lib
	$(CC) $(OBJECTS) $(C_LIBRARY) $(LINKER_FLAGS) $(CLINK) $(PYTHON_LINK)  $(GSL_LINK) -shared -o $@
	$(RM) $(PROJECT).so
	ln -s $(LIBRARY) $(PROJECT).so

sharedlib: $(LIBRARY) | ../../lib
	@{ \
		CURRENT_PYTHON_LIB=`otool -L $(LIBRARY) | grep -i python | cut -d " " -f1 | xargs` ; \
    PYTHON_LIB_NAME=$(PYTHON_LIB_BASE).dylib ; \
		LINK_PYTHON_LIB=$(PYTHON_LIBDIR)/$$PYTHON_LIB_NAME ;\
		if [[ "$$CURRENT_PYTHON_LIB" != "$$LINK_PYTHON_LIB" ]] ; then \
		install_name_tool -change $$CURRENT_PYTHON_LIB $$LINK_PYTHON_LIB $(LIBRARY); \
		fi ;\
	}

../../lib:
	mkdir -p ../../lib

../../lib/libcountpairs_rp_pi_mocks.a: ../../mocks.options ../../common.mk Makefile | ../../lib
	$(MAKE) -C ../DDrppi lib

../../lib/libcountpairs_theta_mocks.a: ../../mocks.options ../../common.mk Makefile | ../../lib
	$(MAKE) -C ../wtheta lib

../../lib/libcountspheres_mocks.a: ../../mocks.options ../../common.mk Makefile | ../../lib
	$(MAKE) -C ../vpf lib


install: $(PHONY_FIX_PY_LINK) ../../lib/$(LIBRARY) | ../../lib

../../lib/$(LIBRARY): $(LIBRARY)
	cp -p $(LIBRARY) $@
	$(RM) ../../lib/$(PROJECT).so 
	ln -s $(LIBRARY) ../../lib/$(PROJECT).so

tests: $(PHONY_FIX_PY_LINK) $(LIBRARY) ../../theory.options ../../common.mk Makefile $(C_LIBRARY) $(OBJECTS) | ../../lib
	python call_correlation_functions_mocks.py

%.o: %.c Makefile ../../common.mk $(C_LIBRARY)
	$(CC) $(OPT) $(CFLAGS) -I../../include $(PYTHON_CFLAGS)  $(INCLUDE) $< -c -o $@

.PHONY: clean clena celan install sharedlib

distclean: clean | ../../lib
	$(RM) ../../lib/$(LIBRARY)
	$(RM) ../../lib/$(PROJECT).so

realclean: distclean

clean:
	$(RM) $(LIBRARY) $(OBJECTS) $(PROJECT).so


celan: clean
clena: clean

