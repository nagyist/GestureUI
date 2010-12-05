../../obj/GClasses/dbg/GKNN.o: GKNN.cpp GKNN.h GLearner.h GData.h GMacros.h GHolders.h GTwt.h \
  GHeap.h GDistribution.h GRand.h GNeighborFinder.h GVec.h GHillClimber.h \
  GOptimizer.h GCluster.h GTransform.h GBitTable.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GKNN.cpp -o ../../obj/GClasses/dbg/GKNN.o
