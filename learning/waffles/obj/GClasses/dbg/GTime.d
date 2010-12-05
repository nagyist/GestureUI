../../obj/GClasses/dbg/GTime.o: GTime.cpp GTime.h GString.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GTime.cpp -o ../../obj/GClasses/dbg/GTime.o
