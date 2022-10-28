ROOT_DIR := ../..
UTILS_DIR:=$(ROOT_DIR)/utils
IO_DIR:=$(ROOT_DIR)/io
INSTALL_HEADERS_DIR := $(ROOT_DIR)/include
INSTALL_LIB_DIR := $(ROOT_DIR)/lib
INSTALL_BIN_DIR := $(ROOT_DIR)/bin

THEORY_DIR := $(ROOT_DIR)/theory
DD_DIR := $(THEORY_DIR)/DD
DDrppi_DIR := $(THEORY_DIR)/DDrppi
WP_DIR := $(THEORY_DIR)/wp
XI_DIR := $(THEORY_DIR)/xi
DDSMU_DIR := $(THEORY_DIR)/DDsmu
VPF_DIR := $(THEORY_DIR)/vpf

DD_LIB := countpairs
DDrppi_LIB := countpairs_rp_pi
WP_LIB := countpairs_wp
XI_LIB := countpairs_xi
DDSMU_LIB := countpairs_s_mu
VPF_LIB := countspheres

GSL_REQUIRED := true
include $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk

PROJECT := _countpairs
SOURCES := $(PROJECT).c
OBJECTS := $(SOURCES:.c=.o)
PYTHON_EXTN := $(PROJECT).so.$(MAJOR).$(MINOR).$(PATCHLEVEL)
C_LIBRARIES := $(DD_DIR)/lib$(DD_LIB).a $(DDrppi_DIR)/lib$(DDrppi_LIB).a $(WP_DIR)/lib$(WP_LIB).a \
             $(XI_DIR)/lib$(XI_LIB).a $(DDSMU_DIR)/lib$(DDSMU_LIB).a $(VPF_DIR)/lib$(VPF_LIB).a
INCL := $(DD_DIR)/$(DD_LIB).h $(DDrppi_DIR)/$(DDrppi_LIB).h $(WP_DIR)/$(WP_LIB).h \
        $(XI_DIR)/$(XI_LIB).h $(DDSMU_DIR)/$(DDSMU_LIB).h $(VPF_DIR)/$(VPF_LIB).h \
        $(UTILS_DIR)/defs.h $(IO_DIR)/io.h $(IO_DIR)/ftread.h \
        $(UTILS_DIR)/utils.h \
	$(UTILS_DIR)/function_precision.h $(UTILS_DIR)/progressbar.h \
        $(UTILS_DIR)/cpu_features.h $(UTILS_DIR)/macros.h
LIB_INCLUDE:=-I$(DD_DIR) -I$(DDrppi_DIR) -I$(WP_DIR) -I$(XI_DIR) -I$(DDSMU_DIR) -I$(VPF_DIR)


all: sharedlib $(SOURCES) $(C_LIBRARIES) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile


$(PROJECT).so: $(PYTHON_EXTN)
	$(RM) $(PROJECT).so
	ln -s $(PYTHON_EXTN) $(PROJECT).so

$(PROJECT).o: $(PROJECT).c $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile $(C_LIBRARIES) $(INCL)
	$(CC) $(CFLAGS) $(LIB_INCLUDE) $(PYTHON_CFLAGS) $(INCLUDE) $< -c -o $@

$(PYTHON_EXTN): $(OBJECTS) $(C_LIBRARIES) $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile
	$(CC) $(OBJECTS) $(C_LIBRARIES) $(LINKER_FLAGS) $(CLINK) $(GSL_LINK) $(PYTHON_LINK) -shared -o $@

$(DD_DIR)/lib$(DD_LIB).a: $(DD_DIR)/*.c $(DD_DIR)/*.c.src $(DD_DIR)/*.h.src $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk
	$(MAKE) -C $(DD_DIR) libs

$(DDrppi_DIR)/lib$(DDrppi_LIB).a: $(DDrppi_DIR)/*.c $(DDrppi_DIR)/*.c.src $(DDrppi_DIR)/*.h.src $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk
	$(MAKE) -C $(DDrppi_DIR) libs

$(WP_DIR)/lib$(WP_LIB).a: $(WP_DIR)/*.c $(WP_DIR)/*.c.src $(WP_DIR)/*.h.src $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk
	$(MAKE) -C $(WP_DIR) libs

$(XI_DIR)/lib$(XI_LIB).a: $(XI_DIR)/*.c $(XI_DIR)/*.c.src $(XI_DIR)/*.h.src $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk
	$(MAKE) -C $(XI_DIR) libs

$(DDSMU_DIR)/lib$(DDSMU_LIB).a: $(DDSMU_DIR)/*.c $(DDSMU_DIR)/*.c.src $(DDSMU_DIR)/*.h.src $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk
	$(MAKE) -C $(DDSMU_DIR) libs

$(VPF_DIR)/lib$(VPF_LIB).a: $(VPF_DIR)/*.c $(VPF_DIR)/*.c.src $(VPF_DIR)/*.h.src $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk
	$(MAKE) -C $(VPF_DIR) libs


install: sharedlib $(INSTALL_LIB_DIR)/$(PYTHON_EXTN)
sharedlib: $(PROJECT).so


$(INSTALL_LIB_DIR)/$(PYTHON_EXTN): $(PYTHON_EXTN) | $(INSTALL_LIB_DIR)
	cp -p $(PYTHON_EXTN) $(INSTALL_LIB_DIR)/

$(INSTALL_LIB_DIR)/$(PROJECT).so:$(INSTALL_LIB_DIR)/$(PYTHON_EXTN)
	$(RM) $(INSTALL_LIB_DIR)/$(PROJECT).so
	cd $(INSTALL_LIB_DIR) && ln -s $(PYTHON_EXTN) $(PROJECT).so

tests: sharedlib $(ROOT_DIR)/theory.options $(ROOT_DIR)/common.mk Makefile $(C_LIBRARIES) $(OBJECTS) | $(ROOT_DIR)/lib
	$(PYTHON) call_correlation_functions.py

clean:
	$(RM) $(OBJECTS) $(PYTHON_EXTN) $(PROJECT).so
	$(RM) -R *.dSYM

distclean:clean
	cd $(INSTALL_LIB_DIR) && $(RM) $(PYTHON_EXTN) $(PROJECT).so
	cd ../../Corrfunc && $(RM) $(PROJECT).so

include $(ROOT_DIR)/rules.mk
