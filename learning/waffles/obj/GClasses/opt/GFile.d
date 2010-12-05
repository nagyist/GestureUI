../../obj/GClasses/opt/GFile.o: GFile.cpp GFile.h GMacros.h GHolders.h GString.h GApp.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GFile.cpp -o ../../obj/GClasses/opt/GFile.o
