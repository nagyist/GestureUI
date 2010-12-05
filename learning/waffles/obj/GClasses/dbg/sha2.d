../../obj/GClasses/dbg/sha2.o: sha2.cpp sha2.h uitypes.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c sha2.cpp -o ../../obj/GClasses/dbg/sha2.o
