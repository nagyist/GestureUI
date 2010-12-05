../../obj/GClasses/opt/GBits.o: GBits.cpp GBits.h GMacros.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GBits.cpp -o ../../obj/GClasses/opt/GBits.o
