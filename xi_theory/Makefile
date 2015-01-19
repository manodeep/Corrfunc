include common.mk

target=logbins 
all: $(target) dirs DD DDrppi wp libs examples python_bindings 

OBJS1= logbins.o utils/utils.o 

logbins: $(OBJS1) Makefile 
	$(CC) $(CFLAGS) -o $@ $(OBJS1)  $(CLINK)

dirs:
	mkdir -p lib bin

.PHONY: clean celna clena celan DD DDrppi libs install wp tests realclean examples dirs

DD:
	$(MAKE) -C xi_of_r

DDrppi:
	$(MAKE) -C xi_rp_pi

wp:
	$(MAKE) -C wp

python_bindings: libs
	$(MAKE) -C python_bindings

### Must create the libraries first before making the examples.
examples: libs
	$(MAKE) -C examples

realclean:
	$(MAKE) clean
	$(RM) bin/$(target)
	$(RM) lib/libcountpairs.a 
	$(RM) lib/libcountpairs_rp_pi.a
	$(RM) lib/libcountpairs_wp.a 
	$(RM) bin/DD bin/DDrppi bin/wprp bin/wp
	$(RM) include/countpairs.h include/countpairs_rp_pi.h include/countpairs_wp.h 
	$(RM) lib/_countpairs.so lib/_countpairs.so.$(MAJOR).0.$(MINOR)

clean:
	$(RM) $(target) $(OBJS1) 
	$(MAKE) -C xi_of_r clean
	$(MAKE) -C xi_rp_pi clean
	$(MAKE) -C wp clean
	$(MAKE) -C examples clean
	$(MAKE) -C python_bindings clean

clena: clean

celan: clean

celna: clean

install: 
	cp -p $(target) bin/
	$(MAKE) -C xi_of_r install
	$(MAKE) -C xi_rp_pi install
	$(MAKE) -C wp install
	$(MAKE) -C python_bindings install

libs: 
	$(MAKE) -C xi_of_r lib
	$(MAKE) -C xi_rp_pi lib
	$(MAKE) -C wp lib


dist:
	hg archive $(DISTNAME).$(MAJOR).0.$(MINOR).tar.gz 

tests:
	$(MAKE) -C xi_of_r tests
	# $(MAKE) -C xi_rp_pi tests
	# $(MAKE) -C wp tests
