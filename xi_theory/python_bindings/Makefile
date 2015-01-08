include ../common.mk

PYTHON_CFLAGS = $(shell python-config --includes) $(shell python -c "import numpy; print '-I' + numpy.__path__[0] + '/core/include/numpy/'")
PYTHON_LINK   = $(shell python-config --ldflags) 
PYTHON_LIBDIR = $(shell python-config --prefix)/lib
CFLAGS += $(PYTHON_CFLAGS)

PROJECT = _countpairs
LIBRARY = $(PROJECT).so.$(MAJOR).0.$(MINOR)
SOURCES = _countpairs.c 
OBJECTS = $(SOURCES:.c=.o)
C_LIBRARY = ../lib/libcountpairs.a ../lib/libcountpairs_rp_pi.a ../lib/libcountpairs_wp.a
INCL   = ../include/defs.h ../include/countpairs.h ../include/countpairs_rp_pi.h ../include/countpairs_wp.h

all: $(LIBRARY)

$(LIBRARY): ../lib $(OBJECTS) lib Makefile ../common.mk 
	$(CC) $(OBJECTS) $(C_LIBRARY) $(CFLAGS) $(CLINK) -L$(PYTHON_LIBDIR) $(PYTHON_LINK) -shared -o $@
	rm -f _countpairs.so 
	ln -s $(LIBRARY) _countpairs.so
	cp -p $(LIBRARY) ../lib/
	rm -f ../lib/_countpairs.so 
	ln -s $(LIBRARY) ../lib/_countpairs.so


../lib:
	mkdir -p ../lib 

lib:
	$(MAKE) -C ../ lib

%.o: %.c $(INCL) Makefile ../common.mk 
	$(CC) $(OPT) $(CFLAGS) $(PYTHON_CFLAGS) $(INCLUDE) $< -c -o $@

.PHONY: clean clena celan install

install: $(LIBRARY)
	cp -p $(LIBRARY) ../lib/
	rm -f ../lib/_countpairs.so 
	ln -s $(LIBRARY) ../lib/_countpairs.so

clean:
	rm -f $(LIBRARY) $(OBJECTS)

celan:
	rm -f $(LIBRARY) $(OBJECTS)

clena:
	rm -f $(LIBRARY) $(OBJECTS) 
