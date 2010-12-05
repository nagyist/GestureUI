../../obj/GClasses/opt/GMacros.o: GMacros.cpp GMacros.h GString.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GMacros.cpp -o ../../obj/GClasses/opt/GMacros.o
