../../obj/GClasses/dbg/GSelfOrganizingMap.o: GSelfOrganizingMap.cpp GSelfOrganizingMap.h \
  GTransform.h GLearner.h GData.h GMacros.h GHolders.h GVec.h GRand.h \
  GMath.h GImage.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GSelfOrganizingMap.cpp -o ../../obj/GClasses/dbg/GSelfOrganizingMap.o
