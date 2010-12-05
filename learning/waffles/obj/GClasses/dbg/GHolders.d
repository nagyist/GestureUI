../../obj/GClasses/dbg/GHolders.o: GHolders.cpp GHolders.h GMacros.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GHolders.cpp -o ../../obj/GClasses/dbg/GHolders.o
