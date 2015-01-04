include ../common.mk

PYTHON_CFLAGS = $(shell python-config --includes) $(shell python -c "import numpy; print '-I' + numpy.__path__[0] + '/core/include/numpy/'")
PYTHON_LINK = $(shell python-config --ldflags) 
CFLAGS += $(PYTHON_CFLAGS)
CLINK  += $(PYTHON_LINK)  -Xlinker -rpath -Xlinker ../lib

PROJECT = _countpairs
LIBRARY = $(PROJECT).so.$(MAJOR).0.$(MINOR)
SOURCES = _countpairs.c 
OBJECTS = $(SOURCES:.c=.o)
C_LIBRARY = ../lib/libcountpairs.so ../lib/libcountpairs_rp_pi.so ../lib/libcountpairs_wp.so 
INCL   = ../include/defs.h ../include/countpairs.h ../include/countpairs_rp_pi.h ../include/countpairs_wp.h

all: $(LIBRARY)

$(LIBRARY): dir $(OBJECTS) lib Makefile ../common.mk 
	$(CC) $(CFLAGS) $(CLINK) $(C_LIBRARY) $(OBJECTS) -shared -o $@
	cp -p $(LIBRARY) ../lib/
	rm -f ../lib/_countpairs.so 
	ln -s ../lib/$(LIBRARY) ../lib/_countpairs.so

dir:
	mkdir -p ../lib 

lib:
	$(MAKE) -C ../ lib

%.o: %.c $(INCL) Makefile ../common.mk 
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) $< -c -o $@

.PHONY: clean clena celan install

install: $(LIBRARY)
	cp -p $(LIBRARY) ../lib/

clean:
	rm -f $(LIBRARY) $(OBJECTS)
	rm -f ../lib/$(LIBRARY)

celan:
	rm -f $(LIBRARY) $(OBJECTS)
	rm -f ../lib/$(LIBRARY)

clena:
	rm -f $(LIBRARY) $(OBJECTS) 
	rm -f ../lib/$(LIBRARY)
