ROOT_DIR := ..
include $(ROOT_DIR)/common.mk
TARGETSRC   := cosmology_params.c gridlink_impl_double.c gridlink_impl_float.c progressbar.c set_cosmo_dist.c utils.c cpu_features.c
TARGETOBJS  := $(TARGETSRC:.c=.o)
INCL  := avx_calls.h cellarray_float.h cellarray_double.h cellarray_mocks.h cosmology_params.h defs.h function_precision.h \
         gridlink_impl_double.h gridlink_impl_float.h gridlink_impl_float.c gridlink_impl_double.c \
         progressbar.h set_cosmo_dist.h set_cosmology.h sglib.h utils.h gridlink_impl.c.src gridlink_impl.h.src cellarray.h.src

all: $(TARGETOBJS) Makefile $(ROOT_DIR)/common.mk $(ROOT_DIR)/theory.options $(ROOT_DIR)/mocks.options

set_cosmo_dist.o: set_cosmo_dist.c $(INCL) $(ROOT_DIR)/mocks.options 
	$(CC) $(CFLAGS) $(GSL_CFLAGS) -c $< -o $@

clean:
	$(RM) $(TARGETOBJS) cellarray_float.h cellarray_double.h gridlink_impl_float.[ch] gridlink_impl_double.[ch]

include $(ROOT_DIR)/rules.mk
