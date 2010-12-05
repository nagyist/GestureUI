../../obj/GClasses/opt/GHistogram.o: GHistogram.cpp GHistogram.h GMacros.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GHistogram.cpp -o ../../obj/GClasses/opt/GHistogram.o
