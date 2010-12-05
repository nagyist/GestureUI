../../obj/GClasses/dbg/GHistogram.o: GHistogram.cpp GHistogram.h GMacros.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GHistogram.cpp -o ../../obj/GClasses/dbg/GHistogram.o
