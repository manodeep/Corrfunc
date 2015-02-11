include ../../mocks.options ../../common.mk

PYTHON_CFLAGS := $(shell python-config --includes) $(shell python -c "import numpy; print '-I' + numpy.__path__[0] + '/core/include/numpy/'")
PYTHON_LIBDIR := $(shell python-config --prefix)/lib 
PYTHON_LINK   := $(shell python-config --ldflags) -Xlinker -rpath -Xlinker $(PYTHON_LIBDIR)

PROJECT := _countpairs_mocks
LIBRARY := $(PROJECT).so.$(MAJOR).0.$(MINOR)
SOURCES := _countpairs_mocks.c 
OBJECTS := $(SOURCES:.c=.o)
C_LIBRARY := ../../lib/libcountpairs_rp_pi_mocks.a ../../lib/libcountpairs_theta_mocks.a ../../lib/libcountspheres_mocks.a 

all: $(LIBRARY) | ../../lib

$(LIBRARY): Makefile ../../mocks.options ../../common.mk $(C_LIBRARY) $(OBJECTS) | ../../lib
	$(CC) $(OBJECTS) $(C_LIBRARY) $(CLINK) -L$(PYTHON_LIBDIR) $(PYTHON_LINK)  $(GSL_LINK) -shared -o $@
	rm -f $(PROJECT).so 
	ln -s $(LIBRARY) $(PROJECT).so
	cp -p $(LIBRARY) ../../lib/
	rm -f ../../lib/$(PROJECT).so
	ln -s $(LIBRARY) ../../lib/$(PROJECT).so

../../lib:
	mkdir -p ../../lib

%.o: %.c Makefile ../../common.mk 
	$(CC) $(OPT) $(CFLAGS) -I../../include $(PYTHON_CFLAGS)  $(INCLUDE) $< -c -o $@

../../lib/libcountpairs_rp_pi_mocks.a: 
	$(MAKE) -C ../DDrppi lib

../../lib/libcountpairs_theta_mocks.a: 
	$(MAKE) -C ../wtheta lib

../../lib/libcountspheres_mocks.a:
	$(MAKE) -C ../vpf lib

install: $(LIBRARY) | ../../lib
	cp -p $(LIBRARY) ../../lib/
	$(RM) ../../lib/$(PROJECT).so
	ln -s $(LIBRARY) ../../lib/$(PROJECT).so

.PHONY: clean clena celan install 

distclean: clean | ../../lib
	$(RM) ../../lib/$(LIBRARY)
	$(RM) ../../lib/$(PROJECT).so

realclean: distclean

clean:
	$(RM) $(LIBRARY) $(OBJECTS) $(PROJECT).so

celan: clean
clena: clean

