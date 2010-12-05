../../obj/GClasses/opt/GEvolutionary.o: GEvolutionary.cpp GMacros.h GEvolutionary.h GOptimizer.h \
  GData.h GHolders.h GBits.h GRand.h GVec.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GEvolutionary.cpp -o ../../obj/GClasses/opt/GEvolutionary.o
