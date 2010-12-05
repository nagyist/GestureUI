../../obj/GClasses/opt/sha2.o: sha2.cpp sha2.h uitypes.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c sha2.cpp -o ../../obj/GClasses/opt/sha2.o
