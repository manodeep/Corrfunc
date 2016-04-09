include ../theory.options ../common.mk

targets:=dirs DD DDrppi wp xi vpf examples 

ifneq ($(COMPILE_PYTHON_EXT), 0)
  PYTHON_DIR := python_bindings
else
  PYTHON_DIR := 
  $(warning $(ccmagenta) Not compiling C extensions for theory. Either python or numpy not available $(ccreset))
endif

all: $(targets) $(PYTHON_DIR) logbins

OBJS1= logbins.o ../utils/utils.o 
logbins: $(OBJS1)  | dirs 
	$(CC) $(OBJS1) $(CLINK) -o $@

%.o: %.c Makefile ../common.mk ../theory.options
	$(CC) $(CFLAGS) -c $< -o $@ 

dirs:  | ../lib ../bin ../include

../lib ../bin ../include:
	mkdir -p $@

.PHONY: clean celna clena celan $(targets) libs install tests distclean realclean distclena realclena

DD:
	$(MAKE) -C xi_of_r

DDrppi:
	$(MAKE) -C xi_rp_pi

wp:
	$(MAKE) -C wp

xi:
	$(MAKE) -C xi

vpf:
	$(MAKE) -C vpf

python_bindings: libs

	$(MAKE) -C python_bindings

### Must create the libraries first before making the examples.
examples: examples/run_correlations 

## Not a phony target
examples/run_correlations: libs
	$(MAKE) -C examples

distclean:realclean
distclena:realclean
realclena:realclean

realclean:
	$(RM) ../bin/logbins
	$(MAKE) clean
	$(MAKE) -C xi_of_r distclean
	$(MAKE) -C xi_rp_pi distclean
	$(MAKE) -C wp distclean
	$(MAKE) -C xi distclean
	$(MAKE) -C vpf distclean
	$(MAKE) -C python_bindings distclean
	$(MAKE) -C tests clean

clean:
	$(RM) logbins $(OBJS1) 
	$(MAKE) -C xi_of_r clean
	$(MAKE) -C xi_rp_pi clean
	$(MAKE) -C wp clean
	$(MAKE) -C xi clean
	$(MAKE) -C vpf clean
	$(MAKE) -C examples clean
	$(MAKE) -C python_bindings clean
	$(MAKE) -C tests clean

clena: clean
celan: clean
celna: clean
install: ../bin/logbins examples $(targets) $(PYTHON_DIR) | dirs
	$(MAKE) -C xi_of_r install
	$(MAKE) -C xi_rp_pi install
	$(MAKE) -C wp install
	$(MAKE) -C xi install
	$(MAKE) -C vpf install
	$(MAKE) -C python_bindings install

../bin/logbins: logbins
	cp -p logbins ../bin

lib: libs

libs: | dirs
	$(MAKE) -C xi_of_r lib
	$(MAKE) -C xi_rp_pi lib
	$(MAKE) -C wp lib
	$(MAKE) -C xi lib
	$(MAKE) -C vpf lib	

test: tests
tests:  
	$(MAKE) -C tests

