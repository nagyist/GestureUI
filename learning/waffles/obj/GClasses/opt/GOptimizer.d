../../obj/GClasses/opt/GOptimizer.o: GOptimizer.cpp GOptimizer.h GMacros.h GData.h GHolders.h \
  GVec.h GRand.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GOptimizer.cpp -o ../../obj/GClasses/opt/GOptimizer.o
