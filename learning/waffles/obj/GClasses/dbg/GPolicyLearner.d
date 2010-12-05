../../obj/GClasses/dbg/GPolicyLearner.o: GPolicyLearner.cpp GPolicyLearner.h GData.h GMacros.h \
  GHolders.h GNeuralNet.h GLearner.h GKNN.h GDecisionTree.h \
  GNeighborFinder.h GOptimizer.h GVec.h GRand.h GHeap.h GHillClimber.h \
  GTwt.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GPolicyLearner.cpp -o ../../obj/GClasses/dbg/GPolicyLearner.o
