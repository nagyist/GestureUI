../../obj/GClasses/dbg/GTwt.o: GTwt.cpp GTwt.h GMacros.h GHeap.h GFile.h GHolders.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GTwt.cpp -o ../../obj/GClasses/dbg/GTwt.o
