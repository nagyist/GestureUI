../../obj/GClasses/dbg/GReinforcement.o: GReinforcement.cpp GReinforcement.h GPolicyLearner.h \
  GData.h GMacros.h GHolders.h GLearner.h GVec.h GNeuralNet.h GKNN.h \
  GRand.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GReinforcement.cpp -o ../../obj/GClasses/dbg/GReinforcement.o
