../../obj/GClasses/opt/GVec.o: GVec.cpp GVec.h GRand.h GMacros.h GData.h GHolders.h GBits.h \
  GTwt.h GHeap.h GMath.h GImage.h GBitTable.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GVec.cpp -o ../../obj/GClasses/opt/GVec.o
