../../obj/GClasses/dbg/GPriorityQueue.o: GPriorityQueue.cpp GPriorityQueue.h GMacros.h \
  GHolders.h GRand.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GPriorityQueue.cpp -o ../../obj/GClasses/dbg/GPriorityQueue.o
