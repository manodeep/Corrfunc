include ../../theory.options ../../common.mk

PYTHON_CFLAGS := $(shell python-config --includes) $(shell python -c "from __future__ import print_function; import numpy; print('-I' + numpy.__path__[0] + '/core/include/numpy/')")
PYTHON_LIBDIR := $(shell python-config --prefix)/lib
PYTHON_LINK   := -L$(PYTHON_LIBDIR) $(shell python-config --ldflags) -Xlinker -rpath -Xlinker $(PYTHON_LIBDIR)

PROJECT := _countpairs
LIBRARY := $(PROJECT).so.$(MAJOR).0.$(MINOR)
SOURCES := _countpairs.c 
OBJECTS := $(SOURCES:.c=.o)
C_LIBRARY := ../../lib/libcountpairs.a ../../lib/libcountpairs_rp_pi.a ../../lib/libcountpairs_wp.a ../../lib/libcountpairs_xi.a ../../lib/libcountspheres.a

all: $(LIBRARY) $(SOURCES) $(C_LIBRARY) ../../theory.options ../../common.mk Makefile | ../../lib

$(LIBRARY): $(SOURCES) ../../theory.options ../../common.mk Makefile $(C_LIBRARY) $(OBJECTS) | ../../lib
	$(CC) $(OBJECTS) $(C_LIBRARY) $(CLINK) $(GSL_LINK) $(PYTHON_LINK) -shared -o $@
	rm -f $(PROJECT).so 
	ln -s $(LIBRARY) $(PROJECT).so
	cp -p $(LIBRARY) ../../lib/
	rm -f ../../lib/$(PROJECT).so 
	ln -s $(LIBRARY) ../../lib/$(PROJECT).so

../../lib:
	mkdir -p ../../lib

%.o: %.c ../../theory.options ../../common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) -I../../include $(PYTHON_CFLAGS)  $(INCLUDE) $< -c -o $@

../../lib/libcountpairs.a: 
	$(MAKE) -C ../xi_of_r/ lib

../../lib/libcountpairs_rp_pi.a: 
	$(MAKE) -C ../xi_rp_pi/ lib

../../lib/libcountpairs_wp.a: 
	$(MAKE) -C ../wp/ lib

../../lib/libcountpairs_xi.a: 
	$(MAKE) -C ../xi/ lib

../../lib/libcountspheres.a: 
	$(MAKE) -C ../vpf/ lib

install: $(LIBRARY) | ../../lib
	cp -p $(LIBRARY) ../../lib/
	$(RM) ../../lib/$(PROJECT).so 
	ln -s $(LIBRARY) ../../lib/$(PROJECT).so

distclean: clean | ../../lib
	$(RM) ../../lib/$(LIBRARY)
	$(RM) ../../lib/$(PROJECT).so

.PHONY: clean clena celan install distclean

clean:
	$(RM) $(LIBRARY) $(OBJECTS) $(PROJECT).so 

celan: clean

clena: clean

