all: dirs theory mocks

dirs: | lib bin include

lib bin include: 
	mkdir -p $@

theory: | dirs
	$(MAKE) -C theory

mocks: | dirs
	$(MAKE) -C mocks

install: | dirs
	$(MAKE) -C theory install
	$(MAKE) -C mocks install

libs: | dirs
	$(MAKE) -C theory libs
	$(MAKE) -C mocks libs

tests:
	$(MAKE) -C theory tests
	$(MAKE) -C mocks tests


.PHONY: clean celna clena celan theory mocks install distclean realclean libs lib

distclean:realclean
distclena:realclean
realclena:realclean

realclean:|dirs
	$(MAKE) -C theory distclean
	$(MAKE) -C mocks distclean
	@{\
		if [ 0 -eq $$(ls -1 lib/lib*.a 2>/dev/null | wc -l) ]; then \
			echo "No static libs in lib/. Removing defs.h " ;\
			rm -f include/defs.h;\
		fi;\
	}
	$(MAKE) -C utils clean
	$(MAKE) -C io clean

clean:
	$(MAKE) -C theory clean
	$(MAKE) -C mocks clean

clena: clean
celan: clean
celna: clean


