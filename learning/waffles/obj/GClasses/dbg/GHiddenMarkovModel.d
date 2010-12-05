../../obj/GClasses/dbg/GHiddenMarkovModel.o: GHiddenMarkovModel.cpp GHiddenMarkovModel.h \
  GMacros.h GHolders.h GVec.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GHiddenMarkovModel.cpp -o ../../obj/GClasses/dbg/GHiddenMarkovModel.o
