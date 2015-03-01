include ../mocks.options ../common.mk

all: dirs DDrppi wtheta vpf libs python_bindings examples

dirs:
	mkdir -p ../lib .,/bin ../include 

.PHONY: clean celna clena celan python_bindings tests DDrppi wtheta distclean realclean distclena realclena dirs test

DDrppi:
	$(MAKE) -C DDrppi

wtheta:
	$(MAKE) -C wtheta

vpf:
	$(MAKE) -C vpf

examples:
	$(MAKE) -C examples

python_bindings:
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
	$(RM) $(target) $(OBJS1) 
	$(MAKE) -C DDrppi clean
	$(MAKE) -C wtheta clean
	$(MAKE) -C vpf clean
	$(MAKE) -C examples clean
	$(MAKE) -C tests clean
	$(MAKE) -C python_bindings clean

clena: clean
celan: clean
celna: clean

install: | dirs 
	$(MAKE) -C DDrppi install
	$(MAKE) -C wtheta install
	$(MAKE) -C vpf install
	$(MAKE) -C python_bindings install 

libs: | dirs
	$(MAKE) -C DDrppi lib
	$(MAKE) -C wtheta lib
	$(MAKE) -C vpf lib

tar: dist

dist:
	hg archive $(DISTNAME).$(MAJOR).0.$(MINOR).tar.gz 

test: tests
tests:
	$(MAKE) -C tests
	cd tests && ./tests_mocks

