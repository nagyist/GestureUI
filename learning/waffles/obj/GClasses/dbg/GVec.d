../../obj/GClasses/dbg/GVec.o: GVec.cpp GVec.h GRand.h GMacros.h GData.h GHolders.h GBits.h \
  GTwt.h GHeap.h GMath.h GImage.h GBitTable.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GVec.cpp -o ../../obj/GClasses/dbg/GVec.o
