../../obj/GClasses/opt/GHashTable.o: GHashTable.cpp GHashTable.h GMacros.h GHolders.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GHashTable.cpp -o ../../obj/GClasses/opt/GHashTable.o
