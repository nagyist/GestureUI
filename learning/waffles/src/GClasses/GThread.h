/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#ifndef __GTHREAD_H__
#define __GTHREAD_H__

#include "GMacros.h"
#ifndef WIN32
#	include <unistd.h>
#	include <sched.h>
#endif // WIN32

namespace GClasses {

#ifdef WIN32
#	define BAD_HANDLE (void*)1
	typedef void* HANDLE;
#else // WIN32
#	ifdef DARWIN
#		define BAD_HANDLE (_opaque_pthread_t*)1
		typedef _opaque_pthread_t* HANDLE;
#	else
#		define HANDLE unsigned long
#		define BAD_HANDLE (size_t)-2
#	endif
#endif // else WIN32

/// A wrapper for PThreads on Linux and for some corresponding WIN32 api on Windows
class GThread
{
public:
	static HANDLE spawnThread(unsigned int (*pFunc)(void*), void* pData);

	/// it may be an error to sleep more than 976ms (1,000,000 / 1024) on Unix
	static void sleep(unsigned int nMiliseconds);
};

} // namespace GClasses

#endif // __GTHREAD_H__
