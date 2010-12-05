/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#ifndef __GBLOB_H__
#define __GBLOB_H__

#include <string.h>
#include "../GClasses/GMacros.h"
#include <string>

namespace GClasses {

/// This class is for deserializing blobs. It takes care of
/// Endianness issues and protects against buffer overruns.
/// This class would be particularly useful for writing a network
/// protocol.
class GBlobIncoming
{
protected:
	int m_nBufferSize;
	int m_nBufferPos;
	unsigned char* m_pBuffer;
	bool m_bDeleteBuffer;

public:
	GBlobIncoming();
	GBlobIncoming(unsigned char* pBuffer, int nSize, bool bDeleteBuffer);
	~GBlobIncoming();

	unsigned char* getBlob() { return m_pBuffer; }
	int getPos() { return m_nBufferPos; }
	void setPos(int pos) { m_nBufferPos = pos; }
	int getBlobSize() { return m_nBufferSize; }
	void setBlob(unsigned char* pBuffer, int nSize, bool bDeleteBuffer);

	/// Pops a blob from the buffer (throws if buffer is too small)
	void get(unsigned char* pData, int nSize);

	/// Pops a single wide char from the buffer (throws if buffer is too small)
	void get(wchar_t* pwc);

	/// Pops a single char from the buffer (throws if buffer is too small)
	void get(char* pc);

	/// Pops an int from the buffer (throws if buffer is too small)
	void get(int* pn);

	/// Pops an unsigned int from the buffer (throws if buffer is too small)
	void get(unsigned int* pui);

	/// Pops an unsigned char from the buffer (throws if buffer is too small)
	void get(unsigned char* puc);

	/// Pops a float from the buffer (throws if buffer is too small)
	void get(float* pf);

	/// Pops a double from the buffer (throws if buffer is too small)
	void get(double* pd);

	/// Pops a string from the buffer
	void get(std::string* pOutString);

	/// Retrieves bytes from within the buffer
	void peek(int nIndex, unsigned char* pData, int nSize);
};

/// This class is for serializing objects
class GBlobOutgoing
{
protected:
	int m_nBufferSize;
	int m_nBufferPos;
	unsigned char* m_pBuffer;
	bool m_bOkToResizeBuffer;

public:
	GBlobOutgoing(int nBufferSize, bool bOkToResizeBuffer);
	~GBlobOutgoing();

	void setPos(int pos) { m_nBufferPos = pos; }
	unsigned char* getBlob() { return m_pBuffer; }
	int getBlobSize() { return m_nBufferPos; }

	/// Pushes a blob into the blob
	void add(const unsigned char* pData, int nSize);

	/// Pushes a wide char into the blob
	void add(const wchar_t wc);

	/// Pushes a char into the blob
	void add(const char c);

	/// Pushes an int into the blob
	void add(const int n);

	/// Pushes an unsigned int into the blob
	void add(const unsigned int n);

	/// Pushes an unsigned char into the blob
	void add(const unsigned char uc);

	/// Pushes a float into the blob
	void add(const float f);

	/// Pushes a double into the blob
	void add(const double d);

	/// Pushes a null-terminated string into the blob
	void add(const char* szString);


	/// Puts bytes into the buffer (overwriting existing data)
	void poke(int nIndex, const unsigned char* pData, int nSize);

	/// Puts an int into the buffer (overwriting existing data)
	void poke(int nIndex, const int n);

protected:
	void resizeBuffer(int nRequiredSize);
};


} // namespace GClasses

#endif // __GBLOB_H__
