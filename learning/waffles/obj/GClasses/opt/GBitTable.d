../../obj/GClasses/opt/GBitTable.o: GBitTable.cpp GBitTable.h GMacros.h GRand.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GBitTable.cpp -o ../../obj/GClasses/opt/GBitTable.o
