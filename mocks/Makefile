include ../mocks.options ../common.mk

TARGETS:= DDrppi_mocks DDtheta_mocks DDsmu_mocks vpf_mocks examples
ifneq ($(COMPILE_PYTHON_EXT), 0)
  TARGETS += python_bindings
else
  $(warning $(ccmagenta) Not compiling C extensions for mocks. Either python or numpy not available $(ccreset))
endif

all: $(TARGETS)

dirs: | ../lib ../bin ../include

../lib ../bin ../include:
	mkdir -p $@

.PHONY: clean celna clena celan $(TARGETS) tests distclean realclean distclena realclena dirs test  python_bindings libs all

DDrppi_mocks DDtheta_mocks vpf_mocks DDsmu_mocks: | dirs
	$(MAKE) -C $@

examples: libs
	$(MAKE) -C examples

python_bindings: libs
	$(MAKE) -C $@

distclean:realclean
distclena:realclean
realclena:realclean

realclean:
	$(MAKE) clean
	for dir in $(TARGETS); do \
		$(MAKE) -C $$dir distclean; \
	done
	$(MAKE) -C ../utils clean
	$(MAKE) -C ../io clean

clean:
	for dir in $(TARGETS); do \
		$(MAKE) -C $$dir clean; \
	done
	$(MAKE) -C tests clean
	$(MAKE) -C python_bindings clean

clena: clean
celan: clean
celna: clean

install: examples | dirs
	for dir in $(TARGETS); do \
		$(MAKE) -C $$dir install; \
	done
ifneq ($(COMPILE_PYTHON_EXT), 0)
	$(MAKE) -C python_bindings install
endif

libs:  | dirs
	for dir in $(TARGETS); do \
		$(MAKE) -C $$dir lib; \
	done

test: tests
tests:
	$(MAKE) -C tests
