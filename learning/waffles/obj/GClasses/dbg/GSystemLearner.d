../../obj/GClasses/dbg/GSystemLearner.o: GSystemLearner.cpp GSystemLearner.h GData.h GMacros.h \
  GHolders.h GActivation.h GPolicyLearner.h GNeuralNet.h GLearner.h \
  GNeighborFinder.h GHillClimber.h GOptimizer.h GVec.h GRand.h \
  GManifold.h GTransform.h GFile.h GHeap.h GImage.h GTwt.h GTime.h \
  GEvolutionary.h GApp.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GSystemLearner.cpp -o ../../obj/GClasses/dbg/GSystemLearner.o
