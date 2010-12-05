../../obj/GClasses/opt/GGraph.o: GGraph.cpp GGraph.h GOptimizer.h GMacros.h GData.h GHolders.h \
  GBitTable.h GHashTable.h GRegion.h GHeap.h GNeighborFinder.h GRand.h \
  GVec.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GGraph.cpp -o ../../obj/GClasses/opt/GGraph.o
