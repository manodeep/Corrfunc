include ../common.mk

target=wp
SRC1= wp.c countpairs_wp.c ../utils/utils.c ../utils/gridlink.c ../io/io.c ../io/ftread.c 
OBJS1  = $(SRC1:.c=.o)
INCL   = countpairs_wp.h ../io/ftread.h ../io/io.h ../utils/utils.h ../utils/gridlink.h \
					../include/function_precision.h ../include/cellarray.h ../include/avx_calls.h \
					../include/defs.h ../include/sglib.h 

all: $(target) $(SRC1) $(INCL)

$(target): $(OBJS1) $(INCL) ../common.mk Makefile 
	$(CC) $(CLINK) $(OBJS1) -o $@

%.o: %.c $(INCL) ../common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

OBJS2=countpairs_wp.o ../utils/gridlink.o ../utils/utils.o 
lib: libcountpairs_wp.a | ../lib

libcountpairs_wp.a: $(OBJS2) | ../lib
	ar rcs $@ $(OBJS2)
	cp -p $@ ../lib/
	sed -e "s/DOUBLE/$(VECTOR_TYPE)/g" countpairs_wp.h > ../include/countpairs_wp.h

../lib:
	mkdir -p ../lib

../bin:
	mkdir -p ../bin

install: $(target) | ../bin
	cp -p $(target) ../bin/

.PHONY: clean clena celan install lib

clean:
	$(RM) $(target) $(OBJS1) $(OBJS2) libcountpairs_wp.a

clena: clean

celan: celan



