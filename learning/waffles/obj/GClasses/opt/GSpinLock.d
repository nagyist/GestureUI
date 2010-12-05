../../obj/GClasses/opt/GSpinLock.o: GSpinLock.cpp GSpinLock.h GMacros.h GThread.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GSpinLock.cpp -o ../../obj/GClasses/opt/GSpinLock.o
