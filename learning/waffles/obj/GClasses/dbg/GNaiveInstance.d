../../obj/GClasses/dbg/GNaiveInstance.o: GNaiveInstance.cpp GNaiveInstance.h GLearner.h GData.h \
  GMacros.h GHolders.h GVec.h GTwt.h GHeap.h GDistribution.h GRand.h \
  GTransform.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GNaiveInstance.cpp -o ../../obj/GClasses/dbg/GNaiveInstance.o
