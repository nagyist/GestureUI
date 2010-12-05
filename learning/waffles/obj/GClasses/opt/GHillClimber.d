../../obj/GClasses/opt/GHillClimber.o: GHillClimber.cpp GHillClimber.h GOptimizer.h GMacros.h \
  GData.h GHolders.h GVec.h GRand.h GImage.h GBitTable.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GHillClimber.cpp -o ../../obj/GClasses/opt/GHillClimber.o
