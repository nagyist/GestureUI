../../obj/GClasses/opt/GParticleSwarm.o: GParticleSwarm.cpp GParticleSwarm.h GOptimizer.h \
  GMacros.h GData.h GHolders.h GVec.h GRand.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GParticleSwarm.cpp -o ../../obj/GClasses/opt/GParticleSwarm.o
