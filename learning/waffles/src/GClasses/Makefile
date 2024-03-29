################
# Paths and Flags
################
SHELL = /bin/bash
TARGET_PATH = ../../lib
TARGET_NAME_OPT = GClasses.a
TARGET_NAME_DBG = GClassesDbg.a
OBJ_PATH = ../../obj/GClasses
UNAME = $(shell uname -s)
ifeq ($(UNAME),Darwin)
#CFLAGS = -D_THREAD_SAFE -DDARWIN -I/sw/include -DBYTE_ORDER_BIG_ENDIAN -no-cpp-precomp
CFLAGS = -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp
else
CFLAGS = -Wall
endif
DBG_CFLAGS = $(CFLAGS) -g -D_DEBUG
OPT_CFLAGS = $(CFLAGS) -O3
UNIX_LIB_DIR = /usr/local/lib
UNIX_INCLUDE_DIR = /usr/local/include

################
# Source
################

CPP_FILES = \
                   GActivation.cpp\
                   GApp.cpp\
                   GBayesianNetwork.cpp\
                   GBits.cpp\
                   GBitTable.cpp\
                   GCluster.cpp\
                   GData.cpp\
                   GDecisionTree.cpp\
                   GDistribution.cpp\
                   GEnsemble.cpp\
                   GFile.cpp\
                   GFourier.cpp\
                   GEvolutionary.cpp\
                   GFunction.cpp\
                   GGraph.cpp\
                   GHashTable.cpp\
                   GHeap.cpp\
                   GHiddenMarkovModel.cpp\
                   GHillClimber.cpp\
                   GHistogram.cpp\
                   GHolders.cpp\
                   GHtml.cpp\
                   GImage.cpp\
                   GKalman.cpp\
                   GKernelTrick.cpp\
                   GKNN.cpp\
                   GLearner.cpp\
                   GLinear.cpp\
                   GMacros.cpp\
                   GManifold.cpp\
                   GMath.cpp\
                   GMixtureOfGaussians.cpp\
                   GNaiveBayes.cpp\
                   GNaiveInstance.cpp\
                   GNeighborFinder.cpp\
                   GNeuralNet.cpp\
                   GOptimizer.cpp\
                   GParticleSwarm.cpp\
                   GPlot.cpp\
                   GPolicyLearner.cpp\
                   GPolynomial.cpp\
                   GPriorityQueue.cpp\
                   GRect.cpp\
                   GRand.cpp\
                   GRegion.cpp\
                   GReinforcement.cpp\
                   GSelfOrganizingMap.cpp\
                   GSparseMatrix.cpp\
                   GSpinLock.cpp\
                   GStabSearch.cpp\
                   GStemmer.cpp\
                   GString.cpp\
                   GSystemLearner.cpp\
                   GText.cpp\
                   GThread.cpp\
                   GTime.cpp\
                   GTransform.cpp\
                   GTwt.cpp\
                   GVec.cpp\
                   sha1.cpp\
                   sha2.cpp

################
# Lists
################

TEMP_LIST_OPT = $(CPP_FILES:%=$(OBJ_PATH)/opt/%)
TEMP_LIST_DBG = $(CPP_FILES:%=$(OBJ_PATH)/dbg/%)
OBJECTS_OPT = $(TEMP_LIST_OPT:%.cpp=%.o)
OBJECTS_DBG = $(TEMP_LIST_DBG:%.cpp=%.o)
DEPS_OPT = $(TEMP_LIST_OPT:%.cpp=%.d)
DEPS_DBG = $(TEMP_LIST_DBG:%.cpp=%.d)

################
# Rules
################

.DELETE_ON_ERROR:

dbg : titledbg $(TARGET_PATH)/$(TARGET_NAME_DBG)

opt : titleopt $(TARGET_PATH)/$(TARGET_NAME_OPT)

titleopt :
	@for ((i=0;i<10;i+=1));do echo "";done
	@echo "Building $(TARGET_NAME_OPT)..."
	@echo ""

titledbg :
	@for ((i=0;i<10;i+=1));do echo "";done
	@echo "Building $(TARGET_NAME_DBG)..."
	@echo ""

usage:
	#
	# Usage:
	#  make usage   (to see this info)
	#  make clean   (to delete all the .o files)
	#  make dbg     (to build a debug version)
	#  make opt     (to build an optimized version)
	#

# This rule makes the optimized ".a" file by using the "ar" tool with all of the optimized ".o" files
$(TARGET_PATH)/$(TARGET_NAME_OPT) : partialcleanopt $(OBJECTS_OPT)
	ar -cr $(TARGET_PATH)/$(TARGET_NAME_OPT) $(OBJECTS_OPT)

# This rule makes the debug ".a" file by using the "ar" tool with all of the debug ".o" files
$(TARGET_PATH)/$(TARGET_NAME_DBG) : partialcleandbg $(OBJECTS_DBG)
	ar -cr $(TARGET_PATH)/$(TARGET_NAME_DBG) $(OBJECTS_DBG)

# This includes all of the ".d" files. Each ".d" file contains a
# generated rule that tells it how to make .o files. (The reason these are generated is so that
# dependencies for these rules can be generated.)
-include $(DEPS_OPT)

-include $(DEPS_DBG)

# This rule makes the optimized ".d" files by using "g++ -MM" with the corresponding ".cpp" file
# The ".d" file will contain a rule that says how to make an optimized ".o" file.
# "$<" refers to the ".cpp" file, and "$@" refers to the ".d" file
$(DEPS_OPT) : $(OBJ_PATH)/opt/%.d : %.cpp
	@echo -e "Computing opt dependencies for $<"
	@-rm -f $$(dirname $@)/$$(basename $@ .d).o
	@if [ ! -d "$$(dirname $@)" ]; then umask 0;mkdir -p "$$(dirname $@)"; fi
	@umask 0;echo -en "$$(dirname $@)/" > $@
	@g++ $(OPT_CFLAGS) -MM $< >> $@
	@echo -e "	umask 0;g++ $(OPT_CFLAGS) -c $< -o $$(dirname $@)/$$(basename $@ .d).o" >> $@

# This rule makes the debug ".d" files by using "g++ -MM" with the corresponding ".cpp" file
# The ".d" file will contain a rule that says how to make a debug ".o" file.
# "$<" refers to the ".cpp" file, and "$@" refers to the ".d" file
$(DEPS_DBG) : $(OBJ_PATH)/dbg/%.d : %.cpp
	@echo -e "Computing dbg dependencies for $<"
	@-rm -f $$(dirname $@)/$$(basename $@ .d).o
	@if [ ! -d "$$(dirname $@)" ]; then umask 0;mkdir -p "$$(dirname $@)"; fi
	@umask 0;echo -en "$$(dirname $@)/" > $@
	@g++ $(DBG_CFLAGS) -MM $< >> $@
	@echo -e "	umask 0;g++ $(DBG_CFLAGS) -c $< -o $$(dirname $@)/$$(basename $@ .d).o" >> $@

partialcleandbg :
	rm -f $(TARGET_PATH)/$(TARGET_NAME_DBG)

partialcleanopt :
	rm -f $(TARGET_PATH)/$(TARGET_NAME_OPT)

clean : partialcleandbg partialcleanopt
	rm -f $(OBJECTS_OPT)
	rm -f $(OBJECTS_DBG)
	rm -f $(DEPS_OPT)
	rm -f $(DEPS_DBG)

install_unix_lib: $(TARGET_PATH)/GClasses.a
	mkdir -m 0755 -p $(UNIX_LIB_DIR)
	mkdir -m 0755 -p $(UNIX_INCLUDE_DIR)/waffles
	install $(TARGET_PATH)/GClasses.a $(UNIX_LIB_DIR)/libGClasses.a
	install *.h $(UNIX_INCLUDE_DIR)/waffles

install: $(TARGET_PATH)/$(TARGET_NAME_OPT) $(TARGET_PATH)/$(TARGET_NAME_DBG)
	mkdir -m 0755 -p $(UNIX_LIB_DIR)
	mkdir -m 0755 -p $(UNIX_INCLUDE_DIR)/GClasses
	install $(TARGET_PATH)/GClasses.a $(UNIX_LIB_DIR)/libGClasses.a
	install $(TARGET_PATH)/GClassesDbg.a $(UNIX_LIB_DIR)/libGClassesDbg.a
	install *.h $(UNIX_INCLUDE_DIR)/GClasses

uninstall:
	rm -f $(UNIX_LIB_DIR)/libGClasses.a
	rm -f $(UNIX_LIB_DIR)/libGClassesDbg.a
	rm -rf $(UNIX_INCLUDE_DIR)/GClasses

.PHONY: clean partialcleandbg partialcleanopt install_unix_lib install uninstall dbg opt
