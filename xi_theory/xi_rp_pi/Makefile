include ../common.mk

target=DDrppi wprp
all: $(target) 

OBJS1= DDrppi.o countpairs_rp_pi.o ../utils/utils.o ../utils/gridlink.o ../io/io.o ../io/ftread.o 


DDrppi: $(OBJS1) Makefile
	$(CC) $(CFLAGS) -o $@ $(OBJS1) $(CLINK)

OBJS2= wprp.o ../utils/utils.o
wprp: $(OBJS2) Makefile 
	$(CC) $(CFLAGS) -o $@ $(OBJS2) $(CLINK)

OBJS3=countpairs_rp_pi.o ../utils/gridlink.o ../utils/utils.o 
lib:
	 $(CC) $(CFLAGS) $(LDFLAGS) $(OBJS3) -shared -o libcountpairs_rp_pi.so 


.PHONY: clean clena celan

clean:
	rm -f $(target) $(OBJS1) $(OBJS2) $(OBJS3) libcountpairs_rp_pi.so 

clena:
	rm -f $(target) $(OBJS1) $(OBJS2) $(OBJS3) libcountpairs_rp_pi.so 

celan:
	rm -f $(target) $(OBJS1) $(OBJS2) $(OBJS3) libcountpairs_rp_pi.so

install:
	cp -p $(target) ../bin/
	cp -p libcountpairs_rp_pi.so ../lib/

