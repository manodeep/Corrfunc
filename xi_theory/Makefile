include common.mk

target=logbins
all: $(target)  DD DDrppi wp lib examples python_bindings 

OBJS1= logbins.o utils/utils.o 

logbins: $(OBJS1) Makefile 
	$(CC) $(CFLAGS) -o $@ $(OBJS1)  $(CLINK)


.PHONY: clean clena celan DD DDrppi lib install wp tests realclean examples

DD:
	$(MAKE) -C xi_of_r

DDrppi:
	$(MAKE) -C xi_rp_pi

wp:
	$(MAKE) -C wp

python_bindings: lib
	$(MAKE) -C python_bindings

### Must create the libraries first before making the examples.
examples: lib
	$(MAKE) -C examples

realclean:
	$(MAKE) clean
	rm -f lib/libcountpairs.a 
	rm -f lib/libcountpairs_rp_pi.a
	rm -f lib/libcountpairs_wp.a 
	rm -f bin/DD bin/DDrppi bin/wprp bin/wp
	rm -f include/countpairs.h include/countpairs_rp_pi.h include/countpairs_wp.h 
	rm -f lib/_countpairs.a


clean:
	rm -f $(target) $(OBJS1) 
	$(MAKE) -C xi_of_r clean
	$(MAKE) -C xi_rp_pi clean
	$(MAKE) -C wp clean
	$(MAKE) -C examples clean
	$(MAKE) -C python_bindings clean

clena:
	rm -f $(target) $(OBJS1) 
	$(MAKE) -C xi_of_r clean
	$(MAKE) -C xi_rp_pi clean
	$(MAKE) -C wp clean
	$(MAKE) -C examples clean
	$(MAKE) -C python_bindings clean

celan:
	rm -f $(target) $(OBJS1)
	$(MAKE) -C xi_of_r clean
	$(MAKE) -C xi_rp_pi clean
	$(MAKE) -C wp clean
	$(MAKE) -C examples clean
	$(MAKE) -C python_bindings clean

install:
	mkdir -p bin
	cp -p $(target) bin/
	$(MAKE) -C xi_of_r install
	$(MAKE) -C xi_rp_pi install
	$(MAKE) -C wp install
	$(MAKE) -C python_bindings install

lib:
	mkdir -p lib
	$(MAKE) -C xi_of_r lib
	$(MAKE) -C xi_rp_pi lib
	$(MAKE) -C wp lib

tests:
	$(MAKE) -C xi_of_r tests
	# $(MAKE) -C xi_rp_pi tests
	# $(MAKE) -C wp tests
