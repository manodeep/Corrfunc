include ../../mocks.options
ifneq (clean,$(findstring clean,$(MAKECMDGOALS)))
  include ../../common.mk
endif

target = run_correlations_mocks

SRC1   := run_correlations_mocks.c ../../io/ftread.c ../../io/io.c  ../../utils/utils.c ../../utils/progressbar.c \
					../../utils/cosmology_params.c 
OBJS1  := $(SRC1:.c=.o)
C_LIBRARIES = ../../lib/libcountpairs_rp_pi_mocks.a ../../lib/libcountpairs_theta_mocks.a  ../../lib/libcountspheres_mocks.a 

INCL :=	../../include/countpairs_rp_pi_mocks.h ../../include/countpairs_theta_mocks.h ../../include/countspheres_mocks.h\
		../../io/ftread.h ../../io/io.h ../../utils/utils.h ../../utils/gridlink_mocks.h \
		../../utils/function_precision.h ../../utils/cellarray_mocks.h ../../utils/avx_calls.h \
		../../utils/defs.h ../../utils/progressbar.h ../../utils/cosmology_params.h 

all: $(target) $(SRC1) $(INCL) ../../mocks.options ../../common.mk Makefile $(C_LIBRARIES)

$(target): Makefile ../../mocks.options ../../common.mk $(INCL) $(C_LIBRARIES) $(OBJS1) 
	$(CC) $(OBJS1) $(C_LIBRARIES) $(CLINK) $(GSL_LINK) -o $@

../../include/countpairs_rp_pi_mocks.h ../../lib/libcountpairs_rp_pi_mocks.a: ../../common.mk ../../mocks.options Makefile
	$(MAKE) -C ../DDrppi lib

../../include/countpairs_theta_mocks.h ../../lib/libcountpairs_theta_mocks.a: ../../common.mk ../../mocks.options Makefile
	$(MAKE) -C ../wtheta lib

../../include/countspheres_mocks.h ../../lib/libcountspheres_mocks.a: ../../common.mk ../../mocks.options Makefile
	$(MAKE) -C ../vpf lib

%.o: %.c ../../common.mk ../../mocks.options Makefile
	$(CC) $(CFLAGS) -I../../include $(INCLUDE) $(OPT) -c $< -o $@

.PHONY: clean clena celan 

clean:
	$(RM) $(target) $(OBJS1) 

clena: clean
celan: clean



