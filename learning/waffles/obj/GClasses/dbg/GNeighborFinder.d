../../obj/GClasses/dbg/GNeighborFinder.o: GNeighborFinder.cpp GNeighborFinder.h GData.h \
  GMacros.h GHolders.h GVec.h GRand.h GPlot.h GRect.h GMath.h \
  GOptimizer.h GHillClimber.h GGraph.h GBitTable.h GTwt.h GHeap.h GKNN.h \
  GLearner.h GTransform.h GNeuralNet.h GImage.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GNeighborFinder.cpp -o ../../obj/GClasses/dbg/GNeighborFinder.o
