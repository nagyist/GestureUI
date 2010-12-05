../../obj/GClasses/dbg/GManifold.o: GManifold.cpp GManifold.h GTransform.h GLearner.h GData.h \
  GMacros.h GHolders.h GRand.h GVec.h GActivation.h GBits.h GBitTable.h \
  GGraph.h GOptimizer.h GHillClimber.h GHeap.h GImage.h GKNN.h GMath.h \
  GNeighborFinder.h GNeuralNet.h GPlot.h GRect.h GSparseMatrix.h GTime.h \
  GTwt.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GManifold.cpp -o ../../obj/GClasses/dbg/GManifold.o
