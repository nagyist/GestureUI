../../obj/GClasses/opt/GHtml.o: GHtml.cpp GHtml.h GHashTable.h GMacros.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GHtml.cpp -o ../../obj/GClasses/opt/GHtml.o
