../../obj/GClasses/opt/GText.o: GText.cpp GText.h GHashTable.h GMacros.h GHeap.h GStemmer.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GText.cpp -o ../../obj/GClasses/opt/GText.o
