include ../common.mk

CC=icc
CLINK=
CFLAGS=  $(OPT) -Wall -Wextra -Wshadow -Wunused -std=c99 -g 

ifeq (icc,$(findstring icc,$(CC)))
CFLAGS += -xhost -ipo -opt-prefetch
else
CFLAGS += -mavx -march=native -funroll-loops -O3
CFLAGS  +=  -Wformat=3  -Wpacked  -Wnested-externs -Wpointer-arith  -Wredundant-decls  -Wfloat-equal -Wcast-qual  
CFLAGS  +=  -Wcast-align -Wmissing-declarations -Wmissing-prototypes  -Wnested-externs -fno-strict-aliasing #-D_POSIX_C_SOURCE=2 -Wpadded -Wstrict-prototypes -Wconversion
CLINK  += -lm
endif

target=DDrppi wprp
all: $(target) 

OBJS1= DDrppi.o gridlink.o countpairs.o \
	utils.o ftread.o 


DDrppi: $(OBJS1) Makefile
	$(CC) $(CFLAGS) -o $@ $(OBJS1) $(CLINK)

OBJS2= wprp.o utils.o
wprp: $(OBJS2) Makefile 
	$(CC) $(CFLAGS) -o $@ $(OBJS2) $(CLINK)

.PHONY: clean clena celan

clean:
	rm -f $(target) $(OBJS1) $(OBJS2) $(OBJS3) $(OBJS4)

clena:
	rm -f $(target) $(OBJS1) $(OBJS2) $(OBJS3) $(OBJS4)

celan:
	rm -f $(target) $(OBJS1) $(OBJS2) $(OBJS3) $(OBJS4)

