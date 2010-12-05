../../obj/GClasses/opt/GPriorityQueue.o: GPriorityQueue.cpp GPriorityQueue.h GMacros.h \
  GHolders.h GRand.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GPriorityQueue.cpp -o ../../obj/GClasses/opt/GPriorityQueue.o
