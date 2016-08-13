ROOT_DIR := ..

include $(ROOT_DIR)/common.mk
SRC   := cosmology_params.c gridlink_impl_double.c gridlink_impl_float.c gridlink_mocks.c progressbar.c set_cosmo_dist.c utils.c cpu_features.c
OBJS  := $(SRC:.c=.o)
INCL  := avx_calls.h cellarray_float.h cellarray_double.h cellarray_mocks.h cosmology_params.h defs.h function_precision.h \
         gridlink_impl_double.h gridlink_impl_float.h gridlink_mocks.h gridlink_impl_float.c gridlink_impl_double.c \
         progressbar.h set_cosmo_dist.h set_cosmology.h sglib.h utils.h gridlink_impl.c.src gridlink_impl.h.src cellarray.h.src

all: $(OBJS) Makefile $(ROOT_DIR)/common.mk $(ROOT_DIR)/theory.options $(ROOT_DIR)/mocks.options

%_float.c: %.c.src Makefile
	@$(ECHO_COMMAND) "/* This file is auto-generated from $(notdir $<) */" > $@
	sed -e "s/DOUBLE/float/g" $< >> $@

%_float.h: %.h.src Makefile
	@$(ECHO_COMMAND) "/* This file is auto-generated from $(notdir $<) */" > $@
	sed -e "s/DOUBLE/float/g"  $< >> $@

%_double.c: %.c.src Makefile
	@$(ECHO_COMMAND) "/* This file is auto-generated from $(notdir $<) */" > $@
	sed -e "s/DOUBLE/double/g"  $< >> $@
%_double.h: %.h.src Makefile
	@$(ECHO_COMMAND) "/* This file is auto-generated from $(notdir $<) */" > $@
	sed -e "s/DOUBLE/double/g"  $< >> $@

%.o: %.c $(INCL) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk $(ROOT_DIR)/mocks.options Makefile
	$(CC) $(CFLAGS) -c $< -o $@

set_cosmo_dist.o: set_cosmo_dist.c $(INCL) $(ROOT_DIR)/mocks.options 
	$(CC) $(CFLAGS) $(GSL_CFLAGS) -c $< -o $@

clean:
	$(RM) $(OBJS) *_float.[ch] *_double.[ch]

clena: clean
celan: celan


