../../obj/GClasses/opt/GSparseMatrix.o: GSparseMatrix.cpp GMacros.h GSparseMatrix.h GData.h \
  GHolders.h GVec.h GFile.h GRand.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GSparseMatrix.cpp -o ../../obj/GClasses/opt/GSparseMatrix.o
