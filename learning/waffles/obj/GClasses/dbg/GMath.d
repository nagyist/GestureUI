../../obj/GClasses/dbg/GMath.o: GMath.cpp GMath.h GMacros.h GVec.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GMath.cpp -o ../../obj/GClasses/dbg/GMath.o
