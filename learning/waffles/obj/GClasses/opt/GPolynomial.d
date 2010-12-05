../../obj/GClasses/opt/GPolynomial.o: GPolynomial.cpp GMacros.h GPolynomial.h GLearner.h GData.h \
  GHolders.h GVec.h GDistribution.h GMath.h GHillClimber.h GOptimizer.h \
  GRand.h GTwt.h GHeap.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GPolynomial.cpp -o ../../obj/GClasses/opt/GPolynomial.o
