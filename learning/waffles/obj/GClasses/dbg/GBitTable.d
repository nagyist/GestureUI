../../obj/GClasses/dbg/GBitTable.o: GBitTable.cpp GBitTable.h GMacros.h GRand.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GBitTable.cpp -o ../../obj/GClasses/dbg/GBitTable.o
