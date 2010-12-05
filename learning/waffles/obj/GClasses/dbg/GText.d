../../obj/GClasses/dbg/GText.o: GText.cpp GText.h GHashTable.h GMacros.h GHeap.h GStemmer.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GText.cpp -o ../../obj/GClasses/dbg/GText.o
