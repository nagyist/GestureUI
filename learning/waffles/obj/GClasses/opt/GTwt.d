../../obj/GClasses/opt/GTwt.o: GTwt.cpp GTwt.h GMacros.h GHeap.h GFile.h GHolders.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GTwt.cpp -o ../../obj/GClasses/opt/GTwt.o
