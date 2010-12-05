../../obj/GClasses/opt/GRect.o: GRect.cpp GRect.h GMacros.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GRect.cpp -o ../../obj/GClasses/opt/GRect.o
