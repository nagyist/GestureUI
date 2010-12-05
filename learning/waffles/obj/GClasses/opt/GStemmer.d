../../obj/GClasses/opt/GStemmer.o: GStemmer.cpp GStemmer.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GStemmer.cpp -o ../../obj/GClasses/opt/GStemmer.o
