../../obj/GClasses/dbg/GCluster.o: GCluster.cpp GCluster.h GTransform.h GLearner.h GData.h \
  GMacros.h GHolders.h GNeighborFinder.h GBitTable.h GHeap.h GMath.h \
  GVec.h GNeuralNet.h GHillClimber.h GOptimizer.h GRand.h GKNN.h GTime.h \
  GGraph.h GTwt.h GImage.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GCluster.cpp -o ../../obj/GClasses/dbg/GCluster.o
