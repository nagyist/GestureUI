../../obj/GClasses/dbg/GImage.o: GImage.cpp GImage.h GMacros.h GBitTable.h GBits.h GFile.h \
  GRect.h GFourier.h GOptimizer.h GData.h GHolders.h GHillClimber.h \
  GVec.h GRand.h GMath.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GImage.cpp -o ../../obj/GClasses/dbg/GImage.o
