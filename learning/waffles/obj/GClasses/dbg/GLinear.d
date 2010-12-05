../../obj/GClasses/dbg/GLinear.o: GLinear.cpp GLinear.h GLearner.h GData.h GMacros.h GHolders.h \
  GTwt.h GHeap.h GTransform.h GRand.h GVec.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GLinear.cpp -o ../../obj/GClasses/dbg/GLinear.o
