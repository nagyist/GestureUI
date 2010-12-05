../../obj/GClasses/dbg/GSpinLock.o: GSpinLock.cpp GSpinLock.h GMacros.h GThread.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GSpinLock.cpp -o ../../obj/GClasses/dbg/GSpinLock.o
