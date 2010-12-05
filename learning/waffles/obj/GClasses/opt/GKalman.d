../../obj/GClasses/opt/GKalman.o: GKalman.cpp GKalman.h GData.h GMacros.h GHolders.h GVec.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GKalman.cpp -o ../../obj/GClasses/opt/GKalman.o
