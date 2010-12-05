../obj/dbg/Gui.o: Gui.cpp Gui.h ../../../src/GClasses/GImage.h \
  ../../../src/GClasses/GMacros.h ../../../src/GClasses/GHolders.h \
  ../../../src/GClasses/GTime.h ../../../src/GClasses/GFile.h \
  ../../../src/GClasses/GThread.h ../../../src/GClasses/GBits.h \
  ../../../src/GClasses/GString.h ../../../src/GClasses/GApp.h
	g++ -I/usr/local/include/SDL -D_THREAD_SAFE -DDARWIN -I/sw/include -I../../../src -no-cpp-precomp -g -D_DEBUG -c Gui.cpp -o ../obj/dbg/Gui.o
