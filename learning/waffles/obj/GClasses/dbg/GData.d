../../obj/GClasses/dbg/GData.o: GData.cpp GData.h GMacros.h GHolders.h GMath.h GDistribution.h \
  GVec.h GFile.h GHeap.h GTwt.h GHashTable.h GBits.h GLearner.h GRand.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GData.cpp -o ../../obj/GClasses/dbg/GData.o
