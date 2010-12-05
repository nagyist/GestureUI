../../obj/GClasses/opt/GStabSearch.o: GStabSearch.cpp GStabSearch.h GOptimizer.h GMacros.h \
  GData.h GHolders.h GRand.h GVec.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GStabSearch.cpp -o ../../obj/GClasses/opt/GStabSearch.o
