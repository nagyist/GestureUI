../obj/dbg/RayTrace.o: RayTrace.cpp RayTrace.h Gui.h ../../../src/GClasses/GRand.h \
  ../../../src/GClasses/GImage.h ../../../src/GClasses/GMacros.h \
  ../../../src/GClasses/GTime.h ../../../src/GClasses/GMacros.h \
  ../../../src/GClasses/GBits.h ../../../src/GClasses/GThread.h \
  GRibParser.h
	g++ -I/usr/local/include/SDL -D_THREAD_SAFE -DDARWIN -I/sw/include -I../../../src -no-cpp-precomp -g -D_DEBUG -c RayTrace.cpp -o ../obj/dbg/RayTrace.o
