include common.mk

target=logbins
all: $(target)  DD DDrppi  

OBJS1= logbins.o utils/utils.o 

logbins: $(OBJS1) Makefile 
	$(CC) $(CFLAGS) -o $@ $(OBJS1)  $(CLINK)


.PHONY: clean clena celan DD DDrppi lib install

DD:
	$(MAKE) -C xi_of_r

DDrppi:
	$(MAKE) -C xi_rp_pi

clean:
	rm -f $(target) $(OBJS1) 
	$(MAKE) -C xi_of_r clean
	$(MAKE) -C xi_rp_pi clean

clena:
	rm -f $(target) $(OBJS1) 
	$(MAKE) -C xi_of_r clean
	$(MAKE) -C xi_rp_pi clean

celan:
	rm -f $(target) $(OBJS1)
	$(MAKE) -C xi_of_r clean
	$(MAKE) -C xi_rp_pi clean


install:
	cp -p $(target) bin/
	$(MAKE) -C xi_of_r install
	$(MAKE) -C xi_rp_pi install

lib:
	$(MAKE) -C xi_of_r lib
	$(MAKE) -C xi_rp_pi lib
