include theory.options common.mk

all: $(target) dirs xi_theory xi_mocks

logbins: $(OBJS1)  | dirs 
	$(CC) $(OBJS1)  $(CLINK) -o $@

dirs:
	mkdir -p lib bin include 

.PHONY: clean celna clena celan xi_theory xi_mocks

xi_theory: | dirs
	$(MAKE) -C xi_theory

xi_mocks: | dirs
	$(MAKE) -C xi_mocks

distclean:realclean

distclena:realclean

realclena:realclean

realclean:
	$(MAKE) -C xi_theory distclean
	$(MAKE) -C xi_mocks distclean

clean:
	$(MAKE) -C xi_theory clean
	$(MAKE) -C xi_mocks clean

clena: clean

celan: clean

celna: clean

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

