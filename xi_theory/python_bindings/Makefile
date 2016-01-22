include ../../theory.options ../../common.mk

PROJECT := _countpairs
LIBRARY := $(PROJECT).so.$(MAJOR).$(MINOR).$(PATCHLEVEL)
SOURCES := _countpairs.c 
OBJECTS := $(SOURCES:.c=.o)
C_LIBRARY := ../../lib/libcountpairs.a ../../lib/libcountpairs_rp_pi.a ../../lib/libcountpairs_wp.a ../../lib/libcountpairs_xi.a ../../lib/libcountspheres.a
PHONY_FIX_PY_LINK := 

ifeq ($(FIX_PYTHON_LINK), 1)
PHONY_FIX_PY_LINK := sharedlib
endif

all: $(PHONY_FIX_PY_LINK) $(LIBRARY) $(SOURCES) $(C_LIBRARY) ../../theory.options ../../common.mk Makefile | ../../lib

$(LIBRARY): $(SOURCES) ../../theory.options ../../common.mk Makefile $(C_LIBRARY) $(OBJECTS) | ../../lib
	$(CC) $(OBJECTS) $(C_LIBRARY) $(CLINK) $(GSL_LINK) $(PYTHON_LINK) -shared -o $@
	$(RM) $(PROJECT).so
	ln -s $(LIBRARY) $(PROJECT).so

sharedlib: $(LIBRARY) | ../../lib
	@{ \
		CURRENT_PYTHON_LIB=`otool -L $(LIBRARY) | grep -i python | cut -d " " -f1 | xargs` ; \
		PYTHON_LIB_NAME=libpython$(PYTHON_VERSION_MAJOR).$(PYTHON_VERSION_MINOR).dylib ; \
		install_name_tool -change $$CURRENT_PYTHON_LIB @loader_path/../lib/$$PYTHON_LIB_NAME $(LIBRARY); \
	}

../../lib:
	mkdir -p ../../lib

../../lib/libcountpairs.a: ../../theory.options ../../common.mk Makefile | ../../lib
	$(MAKE) -C ../xi_of_r/ lib

../../lib/libcountpairs_rp_pi.a: ../../theory.options ../../common.mk Makefile | ../../lib
	$(MAKE) -C ../xi_rp_pi/ lib

../../lib/libcountpairs_wp.a: ../../theory.options ../../common.mk Makefile | ../../lib
	$(MAKE) -C ../wp/ lib

../../lib/libcountpairs_xi.a: ../../theory.options ../../common.mk Makefile | ../../lib
	$(MAKE) -C ../xi/ lib

../../lib/libcountspheres.a: ../../theory.options ../../common.mk Makefile | ../../lib
	$(MAKE) -C ../vpf/ lib

%.o: %.c ../../theory.options ../../common.mk Makefile $(C_LIBRARY)
	$(CC) $(OPT) $(CFLAGS) -I../../include $(PYTHON_CFLAGS)  $(INCLUDE) $< -c -o $@

install: $(PHONY_FIX_PY_LINK) $(LIBRARY) | ../../lib
	cp -p $(LIBRARY) ../../lib/
	$(RM) ../../lib/$(PROJECT).so 
	ln -s $(LIBRARY) ../../lib/$(PROJECT).so

distclean: clean | ../../lib
	$(RM) ../../lib/$(LIBRARY)
	$(RM) ../../lib/$(PROJECT).so

.PHONY: clean celna clena celan install distclean sharedlib

clean:
	$(RM) $(LIBRARY) $(OBJECTS) $(PROJECT).so 

celan: clean
celna: clean
clena: clean

