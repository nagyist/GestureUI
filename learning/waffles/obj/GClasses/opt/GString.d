../../obj/GClasses/opt/GString.o: GString.cpp GString.h GMacros.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GString.cpp -o ../../obj/GClasses/opt/GString.o
