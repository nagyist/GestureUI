../../obj/GClasses/dbg/GActivation.o: GActivation.cpp GActivation.h GMacros.h GMath.h GTwt.h \
  GHeap.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GActivation.cpp -o ../../obj/GClasses/dbg/GActivation.o
