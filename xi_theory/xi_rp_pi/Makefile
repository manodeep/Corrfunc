include ../common.mk

target=DDrppi wprp
all: $(target) 

OBJS1= DDrppi.o countpairs.o ../utils/utils.o ../utils/gridlink.o ../io/io.o ../io/ftread.o 


DDrppi: $(OBJS1) Makefile
	$(CC) $(CFLAGS) -o $@ $(OBJS1) $(CLINK)

OBJS2= wprp.o ../utils/utils.o
wprp: $(OBJS2) Makefile 
	$(CC) $(CFLAGS) -o $@ $(OBJS2) $(CLINK)

.PHONY: clean clena celan

clean:
	rm -f $(target) $(OBJS1) $(OBJS2) 

clena:
	rm -f $(target) $(OBJS1) $(OBJS2) 

celan:
	rm -f $(target) $(OBJS1) $(OBJS2) 

