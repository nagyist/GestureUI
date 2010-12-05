../obj/opt/GRibParser.o: GRibParser.cpp GRibParser.h ../../../src/GClasses/GRand.h \
  ../../../src/GClasses/GImage.h ../../../src/GClasses/GMacros.h
	g++ -I/usr/local/include/SDL -D_THREAD_SAFE -DDARWIN -I/sw/include -I../../../src -no-cpp-precomp -O3 -c GRibParser.cpp -o ../obj/opt/GRibParser.o
