include ../common.mk

target=wp
all: $(target) 

OBJS1= wp.o countpairs_wp.o ../utils/utils.o ../utils/gridlink.o ../io/io.o ../io/ftread.o 

wp: $(OBJS1) Makefile
	$(CC) $(CFLAGS) -o $@ $(OBJS1) $(CLINK)


OBJS2=countpairs_wp.o ../utils/gridlink.o ../utils/utils.o 
lib:
	 $(CC) $(CFLAGS) $(LDFLAGS) $(OBJS2) -shared -o libcountpairs_wp.so 


.PHONY: clean clena celan

clean:
	rm -f $(target) $(OBJS1) $(OBJS2) libcountpairs_wp.so 

clena:
	rm -f $(target) $(OBJS1) $(OBJS2) libcountpairs_wp.so 

celan:
	rm -f $(target) $(OBJS1) $(OBJS2) libcountpairs_wp.so

install:
	cp -p $(target) ../bin/
	cp -p libcountpairs_wp.so ../lib/

