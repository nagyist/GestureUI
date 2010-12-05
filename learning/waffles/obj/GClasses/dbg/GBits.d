../../obj/GClasses/dbg/GBits.o: GBits.cpp GBits.h GMacros.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GBits.cpp -o ../../obj/GClasses/dbg/GBits.o
