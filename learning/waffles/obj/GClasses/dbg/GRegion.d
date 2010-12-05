../../obj/GClasses/dbg/GRegion.o: GRegion.cpp GRegion.h GImage.h GMacros.h GHeap.h GFourier.h \
  GBits.h GRect.h GVec.h GHolders.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GRegion.cpp -o ../../obj/GClasses/dbg/GRegion.o
