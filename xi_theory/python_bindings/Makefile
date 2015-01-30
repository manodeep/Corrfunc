include ../../common.mk

PYTHON_CFLAGS := $(shell python-config --includes) $(shell python -c "import numpy; print '-I' + numpy.__path__[0] + '/core/include/numpy/'")
PYTHON_LINK   := $(shell python-config --ldflags) 
PYTHON_LIBDIR := $(shell python-config --prefix)/lib

PROJECT := _countpairs
LIBRARY := $(PROJECT).so.$(MAJOR).0.$(MINOR)
SOURCES := _countpairs.c 
OBJECTS := $(SOURCES:.c=.o)
C_LIBRARY := ../lib/libcountpairs.a ../lib/libcountpairs_rp_pi.a ../lib/libcountpairs_wp.a

all: $(LIBRARY) | ../lib

$(LIBRARY): Makefile ../../common.mk $(C_LIBRARY) $(OBJECTS) | ../lib
	$(CC) $(OBJECTS) $(C_LIBRARY) $(CLINK) -L$(PYTHON_LIBDIR) $(PYTHON_LINK) -Xlinker -rpath -Xlinker $(PYTHON_LIBDIR) -shared -o $@
	rm -f _countpairs.so 
	ln -s $(LIBRARY) _countpairs.so
	cp -p $(LIBRARY) ../lib/
	rm -f ../lib/_countpairs.so 
	ln -s $(LIBRARY) ../lib/_countpairs.so

../lib:
	mkdir -p ../lib

%.o: %.c Makefile ../../common.mk 
	$(CC) $(OPT) $(CFLAGS) -I../include $(PYTHON_CFLAGS)  $(INCLUDE) $< -c -o $@

../lib/libcountpairs.a: 
	$(MAKE) -C ../xi_of_r/ lib

../lib/libcountpairs_rp_pi.a: 
	$(MAKE) -C ../xi_rp_pi/ lib

../lib/libcountpairs_wp.a: 
	$(MAKE) -C ../wp/ lib

install: $(LIBRARY) | ../lib
	cp -p $(LIBRARY) ../lib/
	$(RM) ../lib/_countpairs.so 
	ln -s $(LIBRARY) ../lib/_countpairs.so

.PHONY: clean clena celan install 

clean:
	$(RM) $(LIBRARY) $(OBJECTS) _countpairs.so 

celan: clean

clena: clean

