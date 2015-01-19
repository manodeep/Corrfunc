include ../common.mk

PYTHON_CFLAGS := $(shell python-config --includes) $(shell python -c "import numpy; print '-I' + numpy.__path__[0] + '/core/include/numpy/'")
PYTHON_LINK   := $(shell python-config --ldflags) 
PYTHON_LIBDIR := $(shell python-config --prefix)/lib

PROJECT := _countpairs
LIBRARY := $(PROJECT).so.$(MAJOR).0.$(MINOR)
SOURCES := _countpairs.c 
OBJECTS = $(SOURCES:.c=.o)
C_LIBRARY := ../lib/libcountpairs.a ../lib/libcountpairs_rp_pi.a ../lib/libcountpairs_wp.a
INCL := ../include/defs.h ../include/countpairs.h ../include/countpairs_rp_pi.h ../include/countpairs_wp.h

all: $(LIBRARY)

$(LIBRARY): libs Makefile ../common.mk $(C_LIBRARY) $(INCL) $(OBJECTS)
	$(CC) $(OBJECTS) $(C_LIBRARY) $(CLINK) -L$(PYTHON_LIBDIR) $(PYTHON_LINK) -Xlinker -rpath -Xlinker $(PYTHON_LIBDIR) -shared -o $@
	@rm -f _countpairs.so 
	@ln -s $(LIBRARY) _countpairs.so
	@cp -p $(LIBRARY) ../lib/
	@rm -f ../lib/_countpairs.so 
	@ln -s $(LIBRARY) ../lib/_countpairs.so

libs:
	$(MAKE) -C ../ libs

%.o: %.c Makefile ../common.mk $(INCL)
	$(CC) $(OPT) $(CFLAGS) $(PYTHON_CFLAGS)  $(INCLUDE) $< -c -o $@

.PHONY: clean clena celan install libs 

install: $(LIBRARY)
	cp -p $(LIBRARY) ../lib/
	$(RM) ../lib/_countpairs.so 
	ln -s $(LIBRARY) ../lib/_countpairs.so

clean:
	$(RM) $(LIBRARY) $(OBJECTS) _countpairs.so 

celan: clean

clena: clean

