../../obj/GClasses/dbg/GHillClimber.o: GHillClimber.cpp GHillClimber.h GOptimizer.h GMacros.h \
  GData.h GHolders.h GVec.h GRand.h GImage.h GBitTable.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GHillClimber.cpp -o ../../obj/GClasses/dbg/GHillClimber.o
