../../obj/GClasses/dbg/GFile.o: GFile.cpp GFile.h GMacros.h GHolders.h GString.h GApp.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GFile.cpp -o ../../obj/GClasses/dbg/GFile.o
