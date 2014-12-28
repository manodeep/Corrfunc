include ../common.mk

target=DD 
all: $(target) 

OBJS1= DD.o countpairs.o ../utils/gridlink.o ../utils/utils.o ../io/ftread.o ../io/io.o 

DD: $(OBJS1) Makefile 
	$(CC) $(CFLAGS) -o $@ $(OBJS1)  $(CLINK)

OBJS2=countpairs.o ../utils/gridlink.o ../utils/utils.o 
lib:
	 $(CC) $(CFLAGS) $(LDFLAGS) $(OBJS2) -shared -o libcountpairs.so 

.PHONY: clean clena celan

clean:
	rm -f $(target) $(OBJS1) $(OBJS2) libcountpairs.so 

clena:
	rm -f $(target) $(OBJS1) $(OBJS2) libcountpairs.so 

celan:
	rm -f $(target) $(OBJS1) $(OBJS2) libcountpairs.so 

install:
	cp -p $(target) ../bin
	cp -p libcountpairs.so ../lib
