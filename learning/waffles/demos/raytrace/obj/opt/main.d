../obj/opt/main.o: main.cpp ../../../src/GClasses/GApp.h \
  ../../../src/GClasses/GMacros.h RayTrace.h Gui.h \
  ../../../src/GClasses/GRand.h ../../../src/GClasses/GImage.h \
  ../../../src/GClasses/GMacros.h
	g++ -I/usr/local/include/SDL -D_THREAD_SAFE -DDARWIN -I/sw/include -I../../../src -no-cpp-precomp -O3 -c main.cpp -o ../obj/opt/main.o
