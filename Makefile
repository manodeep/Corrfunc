include theory.options common.mk

all: dirs xi_theory xi_mocks

# logbins: logbins.c  | dirs 
# 	$(CC) logbins.c  $(CLINK) -o $@

dirs:
	mkdir -p lib bin include 

xi_theory: | dirs
	$(MAKE) -C xi_theory

xi_mocks: | dirs
	$(MAKE) -C xi_mocks

install: | dirs
	$(MAKE) -C xi_theory install
	$(MAKE) -C xi_mocks install

libs: | dirs
	$(MAKE) -C xi_theory libs
	$(MAKE) -C xi_mocks libs

tar: 
	hg archive $(DISTNAME).$(MAJOR).0.$(MINOR).no_version_control.tar.gz -X ".hg*"

dist:
	hg archive $(DISTNAME).$(MAJOR).0.$(MINOR).tar.gz 

tests:
	$(MAKE) -C xi_theory tests
	$(MAKE) -C xi_mocks tests


.PHONY: clean celna clena celan xi_theory xi_mocks install distclean realclean

distclean:realclean
distclena:realclean
realclena:realclean

realclean:
	$(MAKE) -C xi_theory distclean
	$(MAKE) -C xi_mocks distclean

clean:
	$(RM) $(DISTNAME).$(MAJOR).0.$(MINOR).no_version_control.tar.gz $(DISTNAME).$(MAJOR).0.$(MINOR).tar.gz
	$(MAKE) -C xi_theory clean
	$(MAKE) -C xi_mocks clean

clena: clean
celan: clean
celna: clean


