../../obj/GClasses/dbg/GSparseMatrix.o: GSparseMatrix.cpp GMacros.h GSparseMatrix.h GData.h \
  GHolders.h GVec.h GFile.h GRand.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GSparseMatrix.cpp -o ../../obj/GClasses/dbg/GSparseMatrix.o
