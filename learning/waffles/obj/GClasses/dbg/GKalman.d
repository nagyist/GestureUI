../../obj/GClasses/dbg/GKalman.o: GKalman.cpp GKalman.h GData.h GMacros.h GHolders.h GVec.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GKalman.cpp -o ../../obj/GClasses/dbg/GKalman.o
