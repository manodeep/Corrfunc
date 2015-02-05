include ../../mocks.options ../../common.mk

target=DDtheta wtheta
all: $(target) 

ifeq (USE_MKL,$(findstring USE_MKL,$(OPT)))
MKL_CHECK=mkl
endif

INCL   = countpairs_theta.h ../../io/io.h ../../io/ftread.h ../../io/io.h ../../utils/utils.h ../../utils/gridlink_mocks.h \
        ../../utils/function_precision.h ../../utils/cellarray_mocks.h ../../utils/avx_calls.h \
        ../../utils/defs.h ../../utils/progressbar.h

##countpairs_theta_brute_force.c 
SRC1= DDtheta.c  countpairs_theta.c ../../utils/gridlink_mocks.c \
			../../utils/utils.c ../../io/io.c ../../io/ftread.c ../../utils/progressbar.c

OBJS1=$(SRC1:.c=.o) 
DDtheta: $(OBJS1) $(INCL) ../../mocks.options ../../common.mk Makefile
	$(CC) $(OBJS1) $(CLINK) -o $@ 

SRC2=wtheta.c ../../utils/utils.c
OBJS2=$(SRC2:.c=.o)
wtheta: $(OBJS2) $(INCL) ../../mocks.options ../../common.mk Makefile
	$(CC) $(OBJS2) $(CLINK) -o $@ 

%.o: %.c $(INCL) ../../mocks.options ../../common.mk Makefile
	$(CC) $(OPT) $(CFLAGS) $(INCLUDE) -c $< -o $@

# countpairs_theta_brute_force.o: countpairs_theta_brute_force.c $(INCL) ../../mocks.options ../../common.mk Makefile $(MKL_CHECK)
# 	$(CC) $(OPT) $(CFLAGS) $(BLAS_INCLUDE) $(INCLUDE) -c $< -o $@

mkl:
ifndef MKLROOT
	$(error You have requested MKL but MKLROOT environment variable is not defined)
endif

.PHONY: clean clena celan

clean:
	rm -f $(target) $(OBJS1) $(OBJS2) 

clena: clean
celan: clean
