/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#ifndef __GMACROS_H__
#define __GMACROS_H__


#include <stdio.h>
#ifdef WIN32
#	include <BaseTsd.h>
#endif

namespace GClasses {

#define BITS_PER_POINTER (sizeof(void*) * 8)
#define ALIGN_DOWN(p) (((p) / BITS_PER_POINTER) * BITS_PER_POINTER)
#define ALIGN_UP(p) ALIGN_DOWN((p) + BITS_PER_POINTER - 1)

#define INVALID_INDEX ((size_t)-1)

template<typename T>
inline const T& MIN(const T& a, const T& b)
{
	if(b < a)
		return b;
	return a;
}

template<typename T>
inline const T& MAX(const T& a, const T& b)
{
	if(a < b)
		return b;
	return a;
}

template<typename T>
inline T ABS(const T a)
{
	if(a < 0)
		return -a;
	return a;
}


#ifndef UCHAR
#define UCHAR(c) ((c) & (~32))
#endif // UCHAR



void GAssertFailed();
#ifdef _DEBUG
#ifdef WIN32
#define GAssert(x)			\
	{				\
		if(!(x))		\
		{			\
			GAssertFailed();\
			__asm int 3	\
		}			\
	}
#else // WIN32
#define GAssert(x)			\
	{				\
		if(!(x))                \
			GAssertFailed();\
	}
#endif // !WIN32
#else // _DEBUG
#define GAssert(x)	((void)0)
#endif // else _DEBUG






const char* gformat(char c);
const char* gformat(int n);
const char* gformat(unsigned int n);
#if _LP64
const char* gformat(size_t n);
#endif
#if DARWIN
const char* gformat(size_t n);
#endif
const char* gformat(double d, int precision = 14);

void ThrowError(const char* s1);
void ThrowError(const char* s1, const char* s2);
void ThrowError(const char* s1, const char* s2, const char* s3);
void ThrowError(const char* s1, const char* s2, const char* s3, const char* s4);
void ThrowError(const char* s1, const char* s2, const char* s3, const char* s4, const char* s5);
void ThrowError(const char* s1, const char* s2, const char* s3, const char* s4, const char* s5, const char* s6);
void ThrowError(const char* s1, const char* s2, const char* s3, const char* s4, const char* s5, const char* s6, const char* s7);



#define COMPILER_ASSERT(expr)  enum { CompilerAssertAtLine##__LINE__ = sizeof( char[(expr) ? +1 : -1] ) }

/*
/// Placing these on the stack can help catch buffer overruns
class GOverrunSentinel
{
protected:
	unsigned int m_sentinel;

public:
	GOverrunSentinel() : m_sentinel(0x5e47143a)
	{
	}

	~GOverrunSentinel()
	{
		Check();
	}

	void Check()
	{
		if(m_sentinel != 0x5e47143a)
		{
			GAssert(false); // buffer overrun!
			ThrowError("buffer overrun!");
		}
	}
};
*/



// ----------------------------
// Platform Compatability Stuff
// ----------------------------

#ifdef WIN32
typedef UINT_PTR uintptr_t;
// typedef INT_PTR ptrdiff_t;
#else // WIN32
int _stricmp(const char* szA, const char* szB);
int _strnicmp(const char* szA, const char* szB, int len);
long filelength(int filedes);
#endif // !WIN32



} // namespace GClasses

#endif // __GMACROS_H__
