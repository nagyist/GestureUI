../../obj/GClasses/opt/GTransform.o: GTransform.cpp GTransform.h GLearner.h GData.h GMacros.h \
  GHolders.h GTwt.h GHeap.h GVec.h GRand.h GManifold.h GCluster.h \
  GString.h GNeuralNet.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GTransform.cpp -o ../../obj/GClasses/opt/GTransform.o
