all: dirs xi_theory xi_mocks

dirs: | lib bin include

lib bin include: 
	mkdir -p $@

xi_theory: | dirs
	$(MAKE) -C xi_theory

xi_mocks: | dirs
	$(MAKE) -C xi_mocks

install: | dirs
	$(MAKE) -C xi_theory install
	$(MAKE) -C xi_mocks install

libs: | dirs
	$(MAKE) -C xi_theory libs
	$(MAKE) -C xi_mocks libs

tests:
	$(MAKE) -C xi_theory tests
	$(MAKE) -C xi_mocks tests


.PHONY: clean celna clena celan xi_theory xi_mocks install distclean realclean libs lib

distclean:realclean
distclena:realclean
realclena:realclean

realclean:|dirs
	$(MAKE) -C xi_theory distclean
	$(MAKE) -C xi_mocks distclean
	@{\
		if [ 0 -eq $$(ls -1 lib/lib*.a 2>/dev/null | wc -l) ]; then \
			echo "No static libs in lib/. Removing defs.h " ;\
			rm -f include/defs.h;\
		fi;\
	}
	$(MAKE) -C utils clean
	$(MAKE) -C io clean

clean:
	$(MAKE) -C xi_theory clean
	$(MAKE) -C xi_mocks clean

clena: clean
celan: clean
celna: clean


