../../obj/GClasses/opt/GKernelTrick.o: GKernelTrick.cpp GKernelTrick.h GLearner.h GData.h \
  GMacros.h GHolders.h GVec.h GHillClimber.h GOptimizer.h GRand.h \
  GDistribution.h GMath.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GKernelTrick.cpp -o ../../obj/GClasses/opt/GKernelTrick.o
