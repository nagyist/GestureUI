../../obj/GClasses/opt/GEnsemble.o: GEnsemble.cpp GEnsemble.h GLearner.h GData.h GMacros.h \
  GHolders.h GVec.h GDistribution.h GTwt.h GHeap.h GRand.h \
  GDecisionTree.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GEnsemble.cpp -o ../../obj/GClasses/opt/GEnsemble.o
