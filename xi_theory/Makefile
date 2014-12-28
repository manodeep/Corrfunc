include common.mk

target=logbins
all: $(target) 

OBJS1= logbins.o utils/utils.o 

logbins: $(OBJS1) Makefile 
	$(CC) $(CFLAGS) -o $@ $(OBJS1)  $(CLINK)


.PHONY: clean clena celan

clean:
	rm -f $(target) $(OBJS1) 

clena:
	rm -f $(target) $(OBJS1) 

celan:
	rm -f $(target) $(OBJS1) 

install:
	cp -p $(target) bin/
