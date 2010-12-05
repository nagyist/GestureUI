../../obj/GClasses/opt/GBayesianNetwork.o: GBayesianNetwork.cpp GBayesianNetwork.h GRand.h \
  GMath.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GBayesianNetwork.cpp -o ../../obj/GClasses/opt/GBayesianNetwork.o
