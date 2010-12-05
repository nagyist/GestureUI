../../obj/GClasses/opt/GNeuralNet.o: GNeuralNet.cpp GNeuralNet.h GLearner.h GData.h GMacros.h \
  GHolders.h GMath.h GActivation.h GDistribution.h GRand.h GVec.h GTwt.h \
  GHeap.h GHillClimber.h GOptimizer.h GTransform.h GSparseMatrix.h \
  GImage.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GNeuralNet.cpp -o ../../obj/GClasses/opt/GNeuralNet.o
