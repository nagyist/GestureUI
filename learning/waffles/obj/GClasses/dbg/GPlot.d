../../obj/GClasses/dbg/GPlot.o: GPlot.cpp GPlot.h GRect.h GMacros.h GMath.h GVec.h GImage.h \
  GRand.h
	umask 0;g++ -D_THREAD_SAFE -DDARWIN -I/sw/include -no-cpp-precomp -g -D_DEBUG -c GPlot.cpp -o ../../obj/GClasses/dbg/GPlot.o
