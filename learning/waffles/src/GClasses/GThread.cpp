/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#include "GThread.h"
#include "GMacros.h"
#ifdef WIN32
#	include <windows.h>
#else
#	include <pthread.h>
#endif // !WIN32

namespace GClasses {

// static
void GThread::sleep(unsigned int nMiliseconds)
{
#ifdef WIN32
	MSG aMsg;
	while(PeekMessage(&aMsg, NULL, WM_NULL, WM_NULL, PM_REMOVE))
	{
		TranslateMessage(&aMsg);
		DispatchMessage(&aMsg);
	}
	SleepEx(nMiliseconds, 1);
#else // WIN32
	nMiliseconds ? usleep(nMiliseconds*1024) : sched_yield();		// it is an error to sleep for more than 1,000,000
#endif // else WIN32
}

HANDLE GThread::spawnThread(unsigned int (*pFunc)(void*), void* pData)
{
#ifdef WIN32
	unsigned int nID;
	HANDLE hThread = (void*)CreateThread/*_beginthreadex*/(
							NULL,
							0,
							(LPTHREAD_START_ROUTINE)pFunc,
							pData,
							0,
							(unsigned long*)&nID
							);
	if(hThread == BAD_HANDLE)
		ThrowError("Failed to create thread");
	return hThread;
#else // WIN32
	pthread_t thread;
	if(pthread_create(&thread, NULL, (void*(*)(void*))pFunc, pData) != 0)
		ThrowError("Failed to create thread");
	pthread_detach(thread);
	return thread;
#endif // else WIN32
}

} // namespace GClasses

