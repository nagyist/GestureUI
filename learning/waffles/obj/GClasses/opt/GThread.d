../../obj/GClasses/opt/GThread.o: GThread.cpp GThread.h GMacros.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GThread.cpp -o ../../obj/GClasses/opt/GThread.o
