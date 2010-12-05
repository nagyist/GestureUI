../../obj/GClasses/dbg/GFourier.o: GFourier.cpp GFourier.h GMacros.h GHolders.h GMath.h GImage.h \
  GBits.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GFourier.cpp -o ../../obj/GClasses/dbg/GFourier.o
