Build Instructions:
	Linux:
		1- Install the following dependency packages:
			g++
			make
			libpng3-dev
			libsdl1.2-dev
		2- cd waffles/src
		3- sudo make install

		If you also want to build the demo apps:
		4- cd ../src
		5- make opt

	Windows:
		1- Install Microsoft Visual C++ 2008 Express Edition.
		2- File->Open->Project/Solution
		3- Open waffles/src/waffles.sln
		4- Change the configuration from "Debug" to "Release".
		5- Build (F7)
		6- Set the startup app to the demo you want to see.
		7- Run it (F5)

	OSX
		1- Install Fink (a unix package manager).
		1- Install the following packages.
			g++
			make
			libpng3-dev
			libsdl1.2-dev
		2- cd waffles/src
		3- make opt

		If you also want to build the demo apps:
		4- cd ../src
		5- make opt

What to do after you build it:
	If you are already familiar with machine learning, and you know
	what you ultimately want to accomplish, then you might want to
	start by running waffles_wizard (or waffles_wizard.exe if you
	use Windows). It is a graphical tool that will help you to
	construct a command-line-interface (CLI) command to perform
	some machine learning task.

	If you would like to read some documentation, you can start by
	opening waffles/web/index.html in your favorite browser.

	You can run a collection of unit tests by executing
	waffles/bin/test (waffles\bin\test.exe on Windows). These tests
	also provide examples of how to use several of the classes in
	the class library.

	There is also a collection of demo apps in the waffles/demos
	folder that you can build and run. Some of the demos are related
	to machine learning. Others demonstrate classes that you may
	want to combine with your machine learning experiments.

Q and A:
	Linux:
		Q: I get a compiler error in GImage.cpp about png.h
		not found.
		A: You didn't install libpng3-dev as the instructions
		say to do.

		Q: I get a compiler error about sdl-config not found
		A: You didn't install libsdl1.2-dev, did you?

		Q: I'm trying to link GClasses.a with my app, and I
		get a linker error about png stuff not being found.
		A: You need to add "-lpng" after "GClasses.a" to your
		linking command. Yes, order matters on this line.

		Q: How do I uninstall Waffles?
		A: sudo make uninstall

		Q: How do I build optimized binaries?
		A: make opt

		Q: How do I build binaries with debug symbols?
		A: make dbg

		Q: Where can I get help?
		A: Go to http://waffles.sf.net, click on "Forums". Or,
		email me. My email address can be found at at
		http://waffles.sf.net

		Q: Why isn't Waffles in the apt-get/yum repositories?
		A: There's a lot of red-tape involved in getting an
		app into those repositories, and I'm too busy developing
		Waffles to bother with it. If you would like to become
		a package maintainer and do it, that would be a great
		contribution, and I would really appreciate it.

		Q: Why do you include the build dependencies for Windows,
		but not for Linux?
		A: It is easy for Linux developers to get the
		dependencies from the package repositories. This method
		is superior anyway because you get the latest versions.
		If you do not have permissions to install dependencies
		(perhaps because you are using a lab computer), then
		you have a scenario that we do not plan to support.
		You will either need to modify the LFLAGS line in each
		Makefile to look for the dependency libraries in a
		non-standard location, or bother your sysadmin to
		install the dependency libraries. 

		Q: When I run some of the graphical apps, I get an
		error message about a "bad adaptive filter type" or
		"extra compressed data".
		A: This appears to be caused by an old bug in the
		libpng3-dev package. It does not repro with newer
		versions. You should probably update your distro.

	Windows
		Q: I get errors while building GSocket.cpp.
		A: You have an outdated Windows SDK installed. The
		solution is to download and install the latest
		Windows SDK.

		Q: How do I install Waffles?
		A: There is no install process for Windows. Just manually
		copy the executables that you use into some folder in
		your environment's path.

		Q: Where can I get help?
		A: Go to http://waffles.sf.net, click on "Forums". Or,
		email me. My email address can be found at at
		http://waffles.sf.net

		Q: I'm trying to link the GClasses library with my
		code and I get linker errors.
		A: On Windows, the GClasses library is built with the
		static runtime libraries to reduce dll-hell. If you
		try to link it with your code, and your code uses the
		dynamic runtime libraries, you will get linker errors.
		The solution is to rebuild GClasses with dynamic
		runtime libraries, or rebuild your app with static
		runtime libraries. You can't mix static and dynamic
		libraries on Windows.

		Q: May we link the GClasses library with our proprietary
		product, and still keep our source code proprietary?
		A: Read the LGPL. The short answer is yes, but any
		changes you make to my modules must remain open source.
		So, you must link dynamically (.dll or .so), not
		statically (.lib or .a), so you can keep your
		proprietary code in your own modules. I do not have a
		build rule to make GClasses as a dynamic library, so
		you will need to add one, and since this is a change to
		one of my modules, you'll have to contribute it back to
		the project. This is an easy opportunity for you to
		tell your boss that he can keep his trade secrets, and
		benefit from the work of open source developers, as
		long as he lets you make some contributions to the open
		source community. If you use the phrase "win-win" with
		a straight face, this might even lead to a future
		promotion.
