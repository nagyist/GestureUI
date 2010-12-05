../../obj/GClasses/dbg/GPolynomial.o: GPolynomial.cpp GMacros.h GPolynomial.h GLearner.h GData.h \
  GHolders.h GVec.h GDistribution.h GMath.h GHillClimber.h GOptimizer.h \
  GRand.h GTwt.h GHeap.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GPolynomial.cpp -o ../../obj/GClasses/dbg/GPolynomial.o
