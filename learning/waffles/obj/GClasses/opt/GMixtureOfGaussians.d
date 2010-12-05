../../obj/GClasses/opt/GMixtureOfGaussians.o: GMixtureOfGaussians.cpp GMixtureOfGaussians.h \
  GDistribution.h GMacros.h GData.h GHolders.h GVec.h GRand.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GMixtureOfGaussians.cpp -o ../../obj/GClasses/opt/GMixtureOfGaussians.o
