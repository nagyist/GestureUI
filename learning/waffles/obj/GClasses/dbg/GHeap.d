../../obj/GClasses/dbg/GHeap.o: GHeap.cpp GHeap.h GMacros.h GHashTable.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GHeap.cpp -o ../../obj/GClasses/dbg/GHeap.o
