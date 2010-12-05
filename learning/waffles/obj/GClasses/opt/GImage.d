../../obj/GClasses/opt/GImage.o: GImage.cpp GImage.h GMacros.h GBitTable.h GBits.h GFile.h \
  GRect.h GFourier.h GOptimizer.h GData.h GHolders.h GHillClimber.h \
  GVec.h GRand.h GMath.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GImage.cpp -o ../../obj/GClasses/opt/GImage.o
