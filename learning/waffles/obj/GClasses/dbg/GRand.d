../../obj/GClasses/dbg/GRand.o: GRand.cpp GRand.h GMacros.h sha1.h sha2.h uitypes.h GHistogram.h \
  GTime.h GMath.h GVec.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GRand.cpp -o ../../obj/GClasses/dbg/GRand.o
