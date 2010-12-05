../../obj/GClasses/dbg/sha1.o: sha1.cpp sha1.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c sha1.cpp -o ../../obj/GClasses/dbg/sha1.o
