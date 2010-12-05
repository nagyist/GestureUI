../../obj/GClasses/dbg/GDecisionTree.o: GDecisionTree.cpp GDecisionTree.h GLearner.h GData.h \
  GMacros.h GHolders.h GVec.h GPolynomial.h GHillClimber.h GOptimizer.h \
  GRand.h GDistribution.h GTwt.h GHeap.h GTransform.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GDecisionTree.cpp -o ../../obj/GClasses/dbg/GDecisionTree.o
