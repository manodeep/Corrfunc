include ../common.mk

targets:= dirs DDrppi wtheta vpf examples
ifneq ($(COMPILE_PYTHON_EXT), 0)
  targets += python_bindings
else
  $(warning $(ccmagenta) Not compiling C extensions for mocks. Either python or numpy not available $(ccreset))
endif

all: $(targets) 

dirs: | ../lib ../bin ../include

../lib ../bin ../include:
	mkdir -p $@

.PHONY: clean celna clena celan $(targets) tests distclean realclean distclena realclena dirs test  python_bindings libs all 

DDrppi:
	$(MAKE) -C DDrppi

wtheta:
	$(MAKE) -C wtheta

vpf:
	$(MAKE) -C vpf

examples: libs
	$(MAKE) -C examples

python_bindings: libs 
	$(MAKE) -C python_bindings

distclean:realclean
distclena:realclean
realclena:realclean

realclean:
	$(MAKE) clean
	$(MAKE) -C DDrppi distclean
	$(MAKE) -C wtheta distclean
	$(MAKE) -C vpf distclean
	$(MAKE) -C python_bindings distclean

clean:
	$(MAKE) -C DDrppi clean
	$(MAKE) -C wtheta clean
	$(MAKE) -C vpf clean
	$(MAKE) -C examples clean
	$(MAKE) -C tests clean
	$(MAKE) -C python_bindings clean

clena: clean
celan: clean
celna: clean

install: examples | dirs 
	$(MAKE) -C DDrppi install
	$(MAKE) -C wtheta install
	$(MAKE) -C vpf install
	$(MAKE) -C python_bindings install

libs: | dirs
	$(MAKE) -C DDrppi lib
	$(MAKE) -C wtheta lib
	$(MAKE) -C vpf lib

test: tests
tests:
	$(MAKE) -C tests

