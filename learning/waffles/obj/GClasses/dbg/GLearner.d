../../obj/GClasses/dbg/GLearner.o: GLearner.cpp GLearner.h GData.h GMacros.h GHolders.h GVec.h \
  GHeap.h GTwt.h GImage.h GNeuralNet.h GKNN.h GDecisionTree.h \
  GNaiveInstance.h GLinear.h GNaiveBayes.h GEnsemble.h GPolynomial.h \
  GTransform.h GRand.h GPlot.h GRect.h GMath.h GDistribution.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GLearner.cpp -o ../../obj/GClasses/dbg/GLearner.o
