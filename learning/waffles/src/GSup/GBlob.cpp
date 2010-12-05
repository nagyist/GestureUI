/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#include "GBlob.h"
#include "../GClasses/GBits.h"

using namespace GClasses;
using std::string;

GBlobIncoming::GBlobIncoming()
{
	m_pBuffer = NULL;
	m_nBufferSize = 0;
	m_nBufferPos = 0;
	m_bDeleteBuffer = false;
}

GBlobIncoming::GBlobIncoming(unsigned char* pBuffer, int nSize, bool bDeleteBuffer)
{
	m_bDeleteBuffer = false;
	setBlob(pBuffer, nSize, bDeleteBuffer);
}

GBlobIncoming::~GBlobIncoming()
{
	if(m_bDeleteBuffer)
		delete[] m_pBuffer;
}

void GBlobIncoming::setBlob(unsigned char* pBuffer, int nSize, bool bDeleteBuffer)
{
	if(m_bDeleteBuffer)
	{
		GAssert(pBuffer != m_pBuffer); // You gave me ownership of the same buffer twice
		delete(m_pBuffer);
	}
	m_pBuffer = pBuffer;
	m_nBufferSize = nSize;
	m_nBufferPos = 0;
	m_bDeleteBuffer = bDeleteBuffer;
}

void GBlobIncoming::get(unsigned char* pData, int nSize)
{
	if(m_nBufferSize - m_nBufferPos < nSize)
		ThrowError("GBlobIncoming blob is too small to contain the expected data");
	memcpy(pData, &m_pBuffer[m_nBufferPos], nSize);
	m_nBufferPos += nSize;
}

void GBlobIncoming::get(wchar_t* pwc)
{
	get((unsigned char*)pwc, sizeof(wchar_t));
#ifdef WIN32
	GAssert(sizeof(wchar_t) == 2);
	*pwc = GBits::littleEndianToN16((unsigned short)*pwc);
#else // WIN32
	GAssert(sizeof(wchar_t) == 4);
	*pwc = GBits::littleEndianToN32(*pwc);
#endif // !WIN32
}

void GBlobIncoming::get(char* pc)
{
	get((unsigned char*)pc, sizeof(char));
}

void GBlobIncoming::get(int* pn)
{
	get((unsigned char*)pn, sizeof(int));
	*pn = GBits::littleEndianToN32(*pn);
}

void GBlobIncoming::get(unsigned int* pui)
{
	get((unsigned char*)pui, sizeof(unsigned int));
	*pui = GBits::littleEndianToN32(*pui);
}

void GBlobIncoming::get(unsigned char* puc)
{
	get(puc, sizeof(unsigned char));
}

void GBlobIncoming::get(float* pf)
{
	get((unsigned char*)pf, sizeof(float));
	*pf = GBits::littleEndianToR32(*pf);
}

void GBlobIncoming::get(double* pd)
{
	get((unsigned char*)pd, sizeof(double));
	*pd = GBits::littleEndianToR64(*pd);
}

void GBlobIncoming::get(string* pString)
{
	int nLen;
	get(&nLen);
	if(m_nBufferSize - m_nBufferPos < nLen)
		ThrowError("GBlobIncoming blob is too small to contain the expected data");
	pString->assign((const char*)&m_pBuffer[m_nBufferPos], nLen);
	m_nBufferPos += nLen;
}

void GBlobIncoming::peek(int nIndex, unsigned char* pData, int nSize)
{
	if(nSize < 0 || nIndex < 0 || nIndex + nSize > m_nBufferSize)
		ThrowError("GBlobIncoming peek out of range");
	memcpy(pData, &m_pBuffer[nIndex], nSize);
}

// ----------------------------------------------------------------------

GBlobOutgoing::GBlobOutgoing(int nBufferSize, bool bOkToResizeBuffer)
{
	if(nBufferSize > 0)
		m_pBuffer = new unsigned char[nBufferSize];
	else
		m_pBuffer = NULL;
	m_nBufferSize = nBufferSize;
	m_nBufferPos = 0;
	m_bOkToResizeBuffer = bOkToResizeBuffer;
}

GBlobOutgoing::~GBlobOutgoing()
{
	delete[] m_pBuffer;
}

void GBlobOutgoing::resizeBuffer(int nRequiredSize)
{
	if(m_bOkToResizeBuffer)
	{
		int nNewSize = MAX(3 * m_nBufferSize, MAX(1024, nRequiredSize));
		unsigned char* pNewBuffer = new unsigned char[nNewSize];
		memcpy(pNewBuffer, m_pBuffer, m_nBufferPos);
		delete[] m_pBuffer;
		m_pBuffer = pNewBuffer;
		m_nBufferSize = nNewSize;
	}
	else
		ThrowError("GBlobOutgoing buffer too small to hold blob");
}

void GBlobOutgoing::add(const unsigned char* pData, int nSize)
{
	if(m_nBufferSize - m_nBufferPos < nSize)
		resizeBuffer(m_nBufferPos + nSize);
	memcpy(&m_pBuffer[m_nBufferPos], pData, nSize);
	m_nBufferPos += nSize;
}

void GBlobOutgoing::add(const wchar_t wc)
{
#ifdef WIN32
	GAssert(sizeof(wchar_t) == 2);
	wchar_t tmp = GBits::n16ToLittleEndian((unsigned short)wc);
	add((const unsigned char*)&tmp, sizeof(wchar_t));
#else // WIN32
	GAssert(sizeof(wchar_t) == 4);
	wchar_t tmp = GBits::n32ToLittleEndian(wc);
	add((const unsigned char*)&tmp, sizeof(wchar_t));
#endif // !WIN32
}

void GBlobOutgoing::add(const char c)
{
	add((const unsigned char*)&c, sizeof(char));
}

void GBlobOutgoing::add(const int n)
{
	int i = GBits::n32ToLittleEndian(n);
	add((const unsigned char*)&i, sizeof(int));
}

void GBlobOutgoing::add(const unsigned int n)
{
	unsigned int i = GBits::n32ToLittleEndian(n);
	add((const unsigned char*)&i, sizeof(unsigned int));
}

void GBlobOutgoing::add(const unsigned char uc)
{
	add((const unsigned char*)&uc, sizeof(unsigned char));
}

void GBlobOutgoing::add(const float f)
{
	float f2 = GBits::r32ToLittleEndian(f);
	add((const unsigned char*)&f2, sizeof(float));
}

void GBlobOutgoing::add(const double d)
{
	double d2 = GBits::r64ToLittleEndian(d);
	add((const unsigned char*)&d2, sizeof(double));
}

void GBlobOutgoing::add(const char* szString)
{
	int nLen = (int)strlen(szString);
	add(nLen);
	add((const unsigned char*)szString, nLen);
}

void GBlobOutgoing::poke(int nIndex, const unsigned char* pData, int nSize)
{
	if(nSize < 0 || nIndex < 0 || nIndex + nSize > m_nBufferPos)
		ThrowError("GBlobOutgoing poke out of range");
	memcpy(&m_pBuffer[nIndex], pData, nSize);
}

void GBlobOutgoing::poke(int nIndex, const int n)
{
	int tmp = GBits::n32ToLittleEndian(n);
	poke(nIndex, (const unsigned char*)&tmp, sizeof(int));
}
