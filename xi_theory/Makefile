ROOT_DIR:=..

include $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk

PYTHON_EXT_DIR:=python_bindings
TARGETS:=DD DDrppi wp xi vpf examples $(PYTHON_EXT_DIR)

all: $(TARGETS) $(PYTHON_EXT_DIR) logbins

logbins: logbins.c 
	$(CC) $(CFLAGS) $< $(CLINK) -o $@

.PHONY: $(TARGETS) libs install tests distclean realclean distclena realclena clean celna 

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

realclean: clean
	$(RM) $(ROOT_DIR)/bin/logbins
	$(RM) -R *.dSYM
	$(MAKE) -C xi_of_r distclean
	$(MAKE) -C xi_rp_pi distclean
	$(MAKE) -C wp distclean
	$(MAKE) -C xi distclean
	$(MAKE) -C vpf distclean
	$(MAKE) -C python_bindings distclean
	$(MAKE) -C tests clean

clean:
	$(RM) logbins
	$(RM) -R *.dSYM
	$(MAKE) -C xi_of_r clean
	$(MAKE) -C xi_rp_pi clean
	$(MAKE) -C wp clean
	$(MAKE) -C xi clean
	$(MAKE) -C vpf clean
	$(MAKE) -C examples clean
	$(MAKE) -C python_bindings clean
	$(MAKE) -C tests clean
	$(MAKE) -C ../io clean
	$(MAKE) -C ../utils clean

install: $(ROOT_DIR)/bin/logbins examples $(TARGETS) $(PYTHON_EXT_DIR) 
	$(MAKE) -C xi_of_r install
	$(MAKE) -C xi_rp_pi install
	$(MAKE) -C wp install
	$(MAKE) -C xi install
	$(MAKE) -C vpf install
	$(MAKE) -C python_bindings install

$(ROOT_DIR)/bin/logbins: logbins
	cp -p logbins $(ROOT_DIR)/bin

lib: libs

libs: 
	$(MAKE) -C xi_of_r lib
	$(MAKE) -C xi_rp_pi lib
	$(MAKE) -C wp lib
	$(MAKE) -C xi lib
	$(MAKE) -C vpf lib	

test: tests
tests:  
	$(MAKE) -C tests

include $(ROOT_DIR)/rules.mk
