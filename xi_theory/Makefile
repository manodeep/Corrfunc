include common.mk

target=logbins
all: $(target)  DD DDrppi wp

OBJS1= logbins.o utils/utils.o 

logbins: $(OBJS1) Makefile 
	$(CC) $(CFLAGS) -o $@ $(OBJS1)  $(CLINK)


.PHONY: clean clena celan DD DDrppi lib install wp tests realclean

DD:
	$(MAKE) -C xi_of_r

DDrppi:
	$(MAKE) -C xi_rp_pi

wp:
	$(MAKE) -C wp

realclean:
	$(MAKE) clean
	rm -f bin/*
	rm -f lib/*

clean:
	rm -f $(target) $(OBJS1) 
	$(MAKE) -C xi_of_r clean
	$(MAKE) -C xi_rp_pi clean
	$(MAKE) -C wp clean

clena:
	rm -f $(target) $(OBJS1) 
	$(MAKE) -C xi_of_r clean
	$(MAKE) -C xi_rp_pi clean
	$(MAKE) -C wp clean

celan:
	rm -f $(target) $(OBJS1)
	$(MAKE) -C xi_of_r clean
	$(MAKE) -C xi_rp_pi clean
	$(MAKE) -C wp clean

install:
	cp -p $(target) bin/
	$(MAKE) -C xi_of_r install
	$(MAKE) -C xi_rp_pi install
	$(MAKE) -C wp install

lib:
	$(MAKE) -C xi_of_r lib
	$(MAKE) -C xi_rp_pi lib
	$(MAKE) -C wp lib

tests:
	$(MAKE) -C xi_of_r tests
	# $(MAKE) -C xi_rp_pi tests
	# $(MAKE) -C wp tests
