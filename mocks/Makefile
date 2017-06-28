include ../mocks.options ../common.mk

TARGETS:= dirs DDrppi_mocks DDtheta_mocks DDsmu_mocks vpf_mocks examples
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

DDrppi_mocks DDtheta_mocks vpf_mocks DDsmu_mocks:
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
	$(MAKE) -C DDrppi_mocks distclean
	$(MAKE) -C DDsmu_mocks distclean
	$(MAKE) -C DDtheta_mocks distclean
	$(MAKE) -C vpf_mocks distclean
	$(MAKE) -C python_bindings distclean
	$(MAKE) -C ../utils clean
	$(MAKE) -C ../io clean

clean:
	$(MAKE) -C DDrppi_mocks clean
	$(MAKE) -C DDsmu_mocks clean
	$(MAKE) -C DDtheta_mocks clean
	$(MAKE) -C vpf_mocks clean
	$(MAKE) -C examples clean
	$(MAKE) -C tests clean
	$(MAKE) -C python_bindings clean

clena: clean
celan: clean
celna: clean

install: examples | dirs
	$(MAKE) -C DDrppi_mocks install
	$(MAKE) -C DDsmu_mocks install
	$(MAKE) -C DDtheta_mocks install
	$(MAKE) -C vpf_mocks install
	$(MAKE) -C python_bindings install

libs:  | dirs
	$(MAKE) -C DDrppi_mocks lib
	$(MAKE) -C DDsmu_mocks lib
	$(MAKE) -C DDtheta_mocks lib
	$(MAKE) -C vpf_mocks lib

test: tests
tests:
	$(MAKE) -C tests
