../../obj/GClasses/dbg/GHashTable.o: GHashTable.cpp GHashTable.h GMacros.h GHolders.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GHashTable.cpp -o ../../obj/GClasses/dbg/GHashTable.o
