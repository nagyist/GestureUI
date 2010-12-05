../../obj/GClasses/opt/GLinear.o: GLinear.cpp GLinear.h GLearner.h GData.h GMacros.h GHolders.h \
  GTwt.h GHeap.h GTransform.h GRand.h GVec.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -O3 -c GLinear.cpp -o ../../obj/GClasses/opt/GLinear.o
