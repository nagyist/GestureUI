../../obj/GClasses/opt/GHeap.o: GHeap.cpp GHeap.h GMacros.h GHashTable.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GHeap.cpp -o ../../obj/GClasses/opt/GHeap.o
