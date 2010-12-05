../../obj/GClasses/opt/GMath.o: GMath.cpp GMath.h GMacros.h GVec.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GMath.cpp -o ../../obj/GClasses/opt/GMath.o
