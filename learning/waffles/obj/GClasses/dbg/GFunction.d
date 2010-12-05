../../obj/GClasses/dbg/GFunction.o: GFunction.cpp GFunction.h GMacros.h GHolders.h GMath.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GFunction.cpp -o ../../obj/GClasses/dbg/GFunction.o
