/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#include "GMacros.h"
#include <stdarg.h>
#include <wchar.h>
#include <exception>
#ifdef WIN32
#else
#	include <unistd.h>
#endif
#include <signal.h>
#include <sys/stat.h>
#include <string.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include "GString.h"

using std::exception;
using std::string;
using std::cerr;

namespace GClasses {

string g_errorMessage;
class GException : public exception
{
public:
	virtual const char* what() const throw()
	{ return g_errorMessage.c_str(); }
};
GException g_exception;

void ThrowError(const char* s1)
{
	g_errorMessage = s1;

	// Behold! The central location from which all exceptions in this library are thrown!
	// (This might be a good place to put a breakpoint.)
	throw g_exception;
}


void ThrowError(const char* s1, const char* s2)
{
	string s = s1;
	s += s2;
	ThrowError(s.c_str());
}

void ThrowError(const char* s1, const char* s2, const char* s3)
{
	string s = s1;
	s += s2;
	s += s3;
	ThrowError(s.c_str());
}

void ThrowError(const char* s1, const char* s2, const char* s3, const char* s4)
{
	string s = s1;
	s += s2;
	s += s3;
	s += s4;
	ThrowError(s.c_str());
}

void ThrowError(const char* s1, const char* s2, const char* s3, const char* s4, const char* s5)
{
	string s = s1;
	s += s2;
	s += s3;
	s += s4;
	s += s5;
	ThrowError(s.c_str());
}

void ThrowError(const char* s1, const char* s2, const char* s3, const char* s4, const char* s5, const char* s6)
{
	string s = s1;
	s += s2;
	s += s3;
	s += s4;
	s += s5;
	s += s6;
	ThrowError(s.c_str());
}

void ThrowError(const char* s1, const char* s2, const char* s3, const char* s4, const char* s5, const char* s6, const char* s7)
{
	string s = s1;
	s += s2;
	s += s3;
	s += s4;
	s += s5;
	s += s6;
	s += s7;
	ThrowError(s.c_str());
}




#ifdef WIN32
void GAssertFailed()
{
	cerr << "Debug Assert Failed!\n";
	cerr.flush();
}
#else
void GAssertFailed()
{
	cerr << "Debug Assert Failed!\n";
	cerr.flush();
	kill(getpid(), SIGINT);
}

int _stricmp(const char* szA, const char* szB)
{
	while(*szA)
	{
		if((*szA | 32) < (*szB | 32))
			return -1;
		if((*szA | 32) > (*szB | 32))
			return 1;
		szA++;
		szB++;
	}
	if(*szB)
		return -1;
	return 0;
}

int _strnicmp(const char* szA, const char* szB, int len)
{
	int n;
	for(n = 0; n < len; n++)
	{
		if((*szA | 32) < (*szB | 32))
			return -1;
		if((*szA | 32) > (*szB | 32))
			return 1;
		szA++;
		szB++;
	}
	return 0;
}

long filelength(int filedes)
{
	struct stat s;
	if(fstat(filedes, &s) == -1)
		return 0;
	return s.st_size;
}
#endif

char g_convertToStringBuf[256];
int g_convertToStringPos = 0;

char* gformat_get_buf(int size)
{
	GAssert(size < 32);
	char* pBuf = g_convertToStringBuf + g_convertToStringPos;
	g_convertToStringPos += size;
	if(g_convertToStringPos >= 256 - 32)
		g_convertToStringPos = 0;
	return pBuf;
}

const char* gformat(char c)
{
	char* pBuf = gformat_get_buf(2);
	pBuf[0] = c;
	pBuf[1] = '\0';
	return pBuf;
}

const char* gformat(int n)
{
	char* pBuf = gformat_get_buf(14);
	std::ostringstream os;
	os << n;
	string tmp = os.str();
	safe_strcpy(pBuf, tmp.c_str(), 14);
	return pBuf;
}

const char* gformat(unsigned int n)
{
	char* pBuf = gformat_get_buf(14);
	std::ostringstream os;
	os << n;
	string tmp = os.str();
	safe_strcpy(pBuf, tmp.c_str(), 14);
	return pBuf;
}

#if _LP64
const char* gformat(size_t n)
{
	char* pBuf = gformat_get_buf(28);
	std::ostringstream os;
	os << n;
	string tmp = os.str();
	safe_strcpy(pBuf, tmp.c_str(), 28);
	return pBuf;
}
#endif

#if DARWIN
const char* gformat(size_t n)
{
	char* pBuf = gformat_get_buf(28);
	std::ostringstream os;
	os << n;
	string tmp = os.str();
	safe_strcpy(pBuf, tmp.c_str(), 28);
	return pBuf;
}
#endif

const char* gformat(double d, int precision)
{
	char* pBuf = gformat_get_buf(22);
	std::ostringstream os;
	os.precision(precision);
	os << d;
	string tmp = os.str();
	safe_strcpy(pBuf, tmp.c_str(), 22);
	return pBuf;
}

} // namespace GClasses

