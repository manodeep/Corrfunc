include ../mocks.options ../common.mk

all: dirs DDrppi wtheta libs

dirs:
	mkdir -p ../lib .,/bin ../include 

.PHONY: clean celna clena celan DDrppi wtheta distclean realclean distclena realclena dirs 

DDrppi:
	$(MAKE) -C DDrppi

wtheta:
	$(MAKE) -C wtheta

distclean:realclean

distclena:realclean

realclena:realclean

realclean:
	$(MAKE) clean
	$(MAKE) -C DDrppi distclean
	$(MAKE) -C wtheta distclean

clean:
	$(RM) $(target) $(OBJS1) 
	$(MAKE) -C DDrppi clean
	$(MAKE) -C wtheta clean

clena: clean
celan: clean
celna: clean

install: | dirs 
	$(MAKE) -C DDrppi install
	$(MAKE) -C wtheta install

libs: | dirs
	$(MAKE) -C DDrppi lib
	$(MAKE) -C wtheta lib

tar: dist

dist:
	hg archive $(DISTNAME).$(MAJOR).0.$(MINOR).tar.gz 

tests:
	$(MAKE) -C tests
	cd tests && ./test_nonperiodic

