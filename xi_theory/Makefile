include ../theory.options ../common.mk

target=logbins 
all: $(target) dirs DD DDrppi wp xi libs examples python_bindings 

OBJS1= logbins.o ../utils/utils.o 

logbins: $(OBJS1)  | dirs 
	$(CC) $(OBJS1)  $(CLINK) -o $@

%.o: %.c 
	$(CC) $(CFLAGS) -c $< -o $@ 

dirs:
	mkdir -p ../lib .,/bin ../include 

.PHONY: clean celna clena celan DD DDrppi libs install wp tests distclean realclean distclena realclena examples dirs xi

DD:
	$(MAKE) -C xi_of_r

DDrppi:
	$(MAKE) -C xi_rp_pi

wp:
	$(MAKE) -C wp

xi:
	$(MAKE) -C xi

python_bindings: libs
	$(MAKE) -C python_bindings

### Must create the libraries first before making the examples.
examples: libs
	$(MAKE) -C examples

distclean:realclean

distclena:realclean

realclena:realclean

realclean:
	$(RM) ../bin/$(target)
	$(MAKE) clean
	$(MAKE) -C xi_of_r distclean
	$(MAKE) -C xi_rp_pi distclean
	$(MAKE) -C wp distclean
	$(MAKE) -C xi distclean
	$(MAKE) -C vpf distclean
	$(MAKE) -C python_bindings distclean

clean:
	$(RM) $(target) $(OBJS1) 
	$(MAKE) -C xi_of_r clean
	$(MAKE) -C xi_rp_pi clean
	$(MAKE) -C wp clean
	$(MAKE) -C xi clean
	$(MAKE) -C vpf clean
	$(MAKE) -C examples clean
	$(MAKE) -C python_bindings clean
	$(MAKE) -C tests clean

clena: clean
celan: clean
celna: clean
install: $(target) | dirs
	cp -p $(target) ../bin
	$(MAKE) -C xi_of_r install
	$(MAKE) -C xi_rp_pi install
	$(MAKE) -C wp install
	$(MAKE) -C xi install
	$(MAKE) -C vpf install
	$(MAKE) -C python_bindings install

lib: libs

libs: | dirs
	$(MAKE) -C xi_of_r lib
	$(MAKE) -C xi_rp_pi lib
	$(MAKE) -C wp lib
	$(MAKE) -C xi lib
	$(MAKE) -C vpf lib	

tar: dist

dist:
	hg archive $(DISTNAME).$(MAJOR).0.$(MINOR).tar.gz 

tests:
	$(MAKE) -C tests
	cd tests && ./test_periodic
	cd tests && ./test_nonperiodic

