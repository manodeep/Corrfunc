ROOT_DIR:=..

include $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk

PYTHON_EXT_DIR:=
ifneq ($(COMPILE_PYTHON_EXT), 0)
  PYTHON_EXT_DIR:=python_bindings
else
  $(warning $(ccmagenta) Not compiling C extensions for mocks. Either python or numpy not available $(ccreset))
endif

TARGETS:=DD DDrppi DDsmu wp xi vpf examples $(PYTHON_EXT_DIR)

all: $(TARGETS) $(PYTHON_EXT_DIR) logbins

logbins: logbins.c
	$(CC) $(CFLAGS) $< $(CLINK) -o $@

.PHONY: $(TARGETS) libs install tests distclean realclean distclena realclena clean celna

DD DDrppi DDsmu xi wp vpf:
	$(MAKE) -C $@

python_bindings: libs
	$(MAKE) -C $@

### Must create the libraries first before making the examples.
examples: examples/run_correlations

## Not a phony target
examples/run_correlations: libs
	$(MAKE) -C examples

distclena distclean realclena:realclean

realclean:
	$(RM) $(ROOT_DIR)/bin/logbins
	$(RM) -R *.dSYM
	for dir in $(TARGETS); do \
		$(MAKE) -C $$dir distclean; \
	done
	$(MAKE) -C python_bindings distclean
	$(MAKE) -C tests clean

celan clena celna: clean

clean:
	$(RM) logbins
	$(RM) -R *.dSYM
	for dir in $(TARGETS); do \
		$(MAKE) -C $$dir clean; \
	done
	$(MAKE) -C python_bindings clean
	$(MAKE) -C tests clean
	$(MAKE) -C ../io clean
	$(MAKE) -C ../utils clean

install: $(ROOT_DIR)/bin/logbins examples $(TARGETS) $(PYTHON_EXT_DIR)
	for dir in $(TARGETS); do \
		$(MAKE) -C $$dir install; \
	done
ifneq ($(COMPILE_PYTHON_EXT), 0)
	$(MAKE) -C python_bindings install
endif

$(ROOT_DIR)/bin/logbins: logbins
	cp -p logbins $(ROOT_DIR)/bin

lib: libs

libs:
	for dir in $(TARGETS); do \
		$(MAKE) -C $$dir lib; \
	done

test: tests
tests:
	$(MAKE) -C tests

include $(ROOT_DIR)/rules.mk
