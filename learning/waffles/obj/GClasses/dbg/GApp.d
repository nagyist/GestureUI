../../obj/GClasses/dbg/GApp.o: GApp.cpp GApp.h GMacros.h GHolders.h GFile.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GApp.cpp -o ../../obj/GClasses/dbg/GApp.o
