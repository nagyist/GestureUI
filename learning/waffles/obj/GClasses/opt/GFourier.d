../../obj/GClasses/opt/GFourier.o: GFourier.cpp GFourier.h GMacros.h GHolders.h GMath.h GImage.h \
  GBits.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GFourier.cpp -o ../../obj/GClasses/opt/GFourier.o
