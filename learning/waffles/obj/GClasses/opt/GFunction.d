../../obj/GClasses/opt/GFunction.o: GFunction.cpp GFunction.h GMacros.h GHolders.h GMath.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GFunction.cpp -o ../../obj/GClasses/opt/GFunction.o
