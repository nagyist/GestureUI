../../obj/GClasses/opt/GTime.o: GTime.cpp GTime.h GString.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GTime.cpp -o ../../obj/GClasses/opt/GTime.o
