../../obj/GClasses/dbg/GDistribution.o: GDistribution.cpp GDistribution.h GMacros.h GTwt.h \
  GHeap.h GVec.h GRand.h GMath.h GData.h GHolders.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GDistribution.cpp -o ../../obj/GClasses/dbg/GDistribution.o
