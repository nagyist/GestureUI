SUBDIRS= _GClasses _GSup _wizard _generate _plot _learn _transform _test

DBG_SUBDIRS= $(SUBDIRS:_%=DBG_%)
OPT_SUBDIRS= $(SUBDIRS:_%=OPT_%)
CLEAN_SUBDIRS= $(SUBDIRS:_%=CLEAN_%)
INSTALL_SUBDIRS= INSTALL_GClasses INSTALL_GSup INSTALL_wizard INSTALL_generate INSTALL_plot INSTALL_learn INSTALL_transform INSTALL_test
UNINSTALL_SUBDIRS= UNINSTALL_GClasses UNINSTALL_GSup UNINSTALL_wizard UNINSTALL_generate UNINSTALL_plot UNINSTALL_learn UNINSTALL_transform UNINSTALL_test
INSTALL_UNIX_LIB_SUBDIRS= INSTALL_UNIX_LIB_GClasses

.PHONY: dbg opt install uninstall $(DBG_SUBDIRS) $(OPT_SUBDIRS)

dbg: $(DBG_SUBDIRS)

opt: $(OPT_SUBDIRS)

clean: $(CLEAN_SUBDIRS)

configure:
	# There is no need to configure. Just run "make dbg" or "make opt"

install: $(INSTALL_SUBDIRS)

uninstall: $(UNINSTALL_SUBDIRS)

#NOTE: you may want to change the UNIX_LIB_DIR and UNIX_INCLUDE_DIR
#directories in GClasses since they determine where the libraries will
#be installed
install_unix_lib: $(INSTALL_UNIX_LIB_SUBDIRS)

$(INSTALL_UNIX_LIB_SUBDIRS): INSTALL_UNIX_LIB_% :
	$(MAKE) -C $* install_unix_lib

$(DBG_SUBDIRS): DBG_% :
	$(MAKE) -C $* dbg

$(OPT_SUBDIRS): OPT_% :
	$(MAKE) -C $* opt

$(CLEAN_SUBDIRS): CLEAN_% :
	$(MAKE) -C $* clean

$(INSTALL_SUBDIRS): INSTALL_% :
	$(MAKE) -C $* install

$(UNINSTALL_SUBDIRS): UNINSTALL_% :
	$(MAKE) -C $* uninstall
