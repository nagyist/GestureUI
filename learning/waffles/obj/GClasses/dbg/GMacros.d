../../obj/GClasses/dbg/GMacros.o: GMacros.cpp GMacros.h GString.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GMacros.cpp -o ../../obj/GClasses/dbg/GMacros.o
