ROOT_DIR := ..

## While GSL is required in one object file, GSL_REQUIRED is intentionally 
## *NOT* set because that flag will be set/unset by the actual pair-counting. 
## If we set the GSL_REQUIRED flag in this Makefile, then *ALL* directories
## effectively inherit GSL_REQUIRED.

include $(ROOT_DIR)/common.mk
TARGETSRC   := cosmology_params.c gridlink_impl_double.c gridlink_impl_float.c gridlink_mocks_impl_float.c \
               gridlink_mocks_impl_double.c gridink_utils_double.c gridlink_utils_float.c \
               progressbar.c set_cosmo_dist.c utils.c cpu_features.c avx512_calls.c
TARGETOBJS  := $(TARGETSRC:.c=.o)
INCL  := avx512_calls.h avx2_calls.h avx_calls.h sse_calls.h defs.h defs.h function_precision.h cosmology_params.h \
         cellarray_float.h cellarray_double.h cellarray.h.src \
         cellarray_mocks_float.h cellarray_mocks_double.h cellarray_mocks.h.src \
         gridlink_utils_float.h gridlink_utils_double.h gridlink_utils_impl.h.src \
         gridlink_utils_float.c gridlink_utils_double.c gridlink_utils_impl.c.src \
         gridlink_impl_float.c gridlink_impl_double.c \
         gridlink_mocks_impl_float.c gridlink_mocks_impl_double.c \
         gridlink_impl_double.h gridlink_impl_float.h gridlink_impl.c.src gridlink_impl.h.src \
         gridlink_mocks_impl_float.h gridlink_mocks_impl_double.h gridlink_mocks_impl.h.src gridlink_mocks_impl.c.src \
         progressbar.h set_cosmo_dist.h set_cosmology.h sglib.h utils.h \
         weight_functions_double.h weight_functions_float.h weight_functions.h.src \
	 weight_defs_double.h weight_defs_float.h weight_defs.h.src \
         cell_pair_float.h cell_pair_double.h cell_pair.h.src

all: $(TARGETOBJS) Makefile $(ROOT_DIR)/common.mk $(ROOT_DIR)/theory.options $(ROOT_DIR)/mocks.options

set_cosmo_dist.o: set_cosmo_dist.c $(ROOT_DIR)/mocks.options
	$(CC) $(CFLAGS) $(GSL_CFLAGS) -c $< -o $@

clean:
	$(RM) $(TARGETOBJS) cellarray_float.h cellarray_double.h gridlink_impl_float.[ch] gridlink_impl_double.[ch] cellarray_mocks_float.h cellarray_mocks_double.h gridlink_mocks_impl_float.[ch] gridlink_utils_double.[ch] gridlink_utils_float.[ch] gridlink_mocks_impl_double.[ch] weight_functions_double.h weight_functions_float.h weight_defs_double.h weight_defs_float.h cell_pair_double.h cell_pair_float.h

include $(ROOT_DIR)/rules.mk
