../../obj/GClasses/opt/sha1.o: sha1.cpp sha1.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c sha1.cpp -o ../../obj/GClasses/opt/sha1.o
