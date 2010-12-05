../../obj/GClasses/dbg/GThread.o: GThread.cpp GThread.h GMacros.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GThread.cpp -o ../../obj/GClasses/dbg/GThread.o
