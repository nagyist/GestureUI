/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#include "GFile.h"
#include "GMacros.h"
#include "GHolders.h"
//#include "GDirList.h"
#include "GString.h"
#include "GApp.h"
#ifdef WIN32
#	include <windows.h>
#	include <shlobj.h> // to get users' application data dir
#	include <sys/utime.h> // utime
#	include <direct.h>
#	include <io.h> // for "filelength"
#	include <process.h>
#else
#	include <unistd.h>
#	include <utime.h> // utime, which sets file times
#endif // WIN32
#include <stdio.h>
#include <sys/types.h>
#include <time.h>
#include <stdlib.h> // getenv
#include <string.h>
#include <sstream>
#include <fstream>


using namespace GClasses;
using std::string;

bool GFile::doesFileExist(const char *szFilename)
{
	return(access(szFilename, 0) == 0);
}

bool GFile::doesDirExist(const char *szDir)
{
	char szBuff[256];
	char* pCurDir = getcwd(szBuff, 255);
	int nVal = chdir(szDir);
	if(chdir(pCurDir) != 0)
		ThrowError("Failed to restore current dir");
	return(nVal == 0);
}

bool GFile::makeDir(const char* szPath)
{
	char szDir[256];
	safe_strcpy(szDir, szPath, 256);
	bool bOK = false;
	int n;
	for(n = 0; ; n++)
	{
		if(szDir[n] == '/' || szDir[n] == '\\' || szDir[n] == '\0')
		{
			char cTmp = szDir[n];
			szDir[n] = '\0';
#ifdef WIN32
			bOK = (mkdir(szDir) == 0);
#else // WIN32
			bOK = (mkdir(szDir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == 0); // read/write/search permissions for owner and group, and with read/search permissions for others
#endif // !WIN32
			szDir[n] = cTmp;
			if(cTmp == '\0')
				break;
		}
	}
	return bOK;
}

// Find the last slash and return what's past that
const char* GFile::clipPath(const char* szBuff)
{
	int n = 0;
	int i = 0;
	while(szBuff[i] != 0)
	{
		if(szBuff[i] == '\\' || szBuff[i] == '/')
			n = i + 1;
		i++;
	}
	return(szBuff + n);
}

// Find the last slash and set it to '\0'
char* GFile::clipFilename(char* szBuff)
{
	int n = -1;
	int i = 0;
	while(szBuff[i] != 0)
	{
		if(szBuff[i] == '\\' || szBuff[i] == '/')
			n = i;
		i++;
	}
   if(n > -1)
      szBuff[n + 1] = '\0';
	return szBuff;
}

bool GFile::copyFile(const char* szSrcPath, const char* szDestPath)
{
	FILE* pSrc = fopen(szSrcPath, "rb");
	if(!pSrc)
		return false;

	// create the subdirectory
	GTEMPBUF(char, szSubdirectoryString, (int)strlen(szDestPath)+ 1);
	strcpy(szSubdirectoryString, szDestPath);
	clipFilename(szSubdirectoryString);
	makeDir(szSubdirectoryString);
	
	FILE* pDest = fopen(szDestPath, "wb");
	if(!pDest)
	{
		fclose(pSrc);
		return false;
	}
	char szBuf[1024];
	int nFileSize = filelength(fileno(pSrc));
	int nChunkSize;
	while(nFileSize > 0)
	{
		nChunkSize = MIN(nFileSize, 1024);
		if(fread(szBuf, nChunkSize, 1, pSrc) != 1)
			ThrowError("Error reading");
		if(fwrite(szBuf, nChunkSize, 1, pDest) != 1)
			ThrowError("Error writing");
		nFileSize -= nChunkSize;
	}
	fclose(pSrc);
	fclose(pDest);
	return true;
}

bool GFile::localStorageDirectory(char *toHere)
{
	char *szReturnValue = NULL;
	toHere[0] = '\0';
#ifdef WIN32
	TCHAR szPath[MAX_PATH];
	if(SUCCEEDED(SHGetFolderPath(NULL, 
				CSIDL_LOCAL_APPDATA, 
				NULL, 
				0, 
				szPath))) 
	{
		strcpy(toHere, szPath);
		return true;
	}
	else
		return false;
#else
#	ifndef __linux__
	// DARWIN
	GApp::appPath(toHere, 360, true);
	return true;
#	else
	szReturnValue = getenv("HOME");
	if(!szReturnValue)
		szReturnValue = getenv("HOMEPATH");
	if(!szReturnValue)
		szReturnValue = getenv("USERPROFILE");

	strcpy(toHere, szReturnValue);

	return szReturnValue != NULL; // if we were successful this will have been set to something...
#	endif
#endif
}

/*static*/ char* GFile::loadFile(const char* szFilename, size_t* pnSize)
{
	std::ifstream s;
	s.exceptions(std::ios::failbit|std::ios::badbit);
	try
	{
		s.open(szFilename, std::ios::binary);
		s.seekg(0, std::ios::end);
		*pnSize = s.tellg();
		s.seekg(0, std::ios::beg);
	}
	catch(const std::exception&)
	{
		if(GFile::doesFileExist(szFilename))
			ThrowError("Error while trying to open the existing file: ", szFilename);
		else
			ThrowError("File not found: ", szFilename);
	}
	char* pBuf = new char[*pnSize + 1];
	ArrayHolder<char> hBuf(pBuf);
	s.read(pBuf, *pnSize);
	pBuf[*pnSize] = '\0'; // null terminate the file. todo: this is unnecessary. Make sure nothing relies on it, and remove it.
	return hBuf.release();
}

/*static*/ void GFile::saveFile(const char* pBuf, size_t size, const char* szFilename)
{
	std::ofstream s;
	s.exceptions(std::ios::failbit|std::ios::badbit);
	try
	{
		s.open(szFilename, std::ios::binary);
	}
	catch(const std::exception&)
	{
		ThrowError("Error creating file: ", szFilename);
	}
	s.write(pBuf, size);
}

void GFile::condensePath(char* szPath)
{
	int n;
#ifdef WIN32
	for(n = 0; szPath[n] != '\0'; n++)
	{
		if(szPath[n] == '\\')
			szPath[n] = '/';
	}
#endif // WIN32
	int nPrevSlash = -1;
	int nPrevPrevSlash = -1;
	for(n = 0; szPath[n] != '\0'; n++)
	{
		if(szPath[n] == '/')
		{
			nPrevPrevSlash = nPrevSlash;
			nPrevSlash = n;
			if(strncmp(szPath + n, "/./", 3) == 0)
			{
				int nDelSize = 2;
				int i;
				for(i = n + 1; ; i++)
				{
					szPath[i] = szPath[i + nDelSize];
					if(szPath[i] == '\0')
						break;
				}
				n--; // so we'll catch the current slash the next time around the loop -- handle it next time!
				
			}
			else if(nPrevPrevSlash >= 0 && strncmp(szPath + n, "/../", 4) == 0)
			{
				int nDelSize = n - nPrevPrevSlash + 3;
				int i;
				for(i = nPrevPrevSlash; ; i++)
				{
					szPath[i] = szPath[i + nDelSize];
					if(szPath[i] == '\0')
						break;
				}
				nPrevSlash = -1;
				nPrevPrevSlash = -1;
				n = -1;
			}
				
		}
	}
}

// static
time_t GFile::modifiedTime(const char* szFilename)
{
	struct stat buffer;
	/*int status = */stat(szFilename, &buffer);
	// todo: check return value
	return buffer.st_mtime; // number of seconds since 0000 UTC 1970... :-)
}

// static
void GFile::setModifiedTime(const char* szFilename, time_t t)
{
	struct stat bufferIn;
	/*int status = */stat(szFilename, &bufferIn);
	struct utimbuf bufferOut;
	bufferOut.actime = bufferIn.st_atime;
	bufferOut.modtime = t;
	utime(szFilename, &bufferOut);
}
/*
// static
bool GFile::shredFile(const char* szFilename)
{
#ifdef WIN32
	char szTmp[4096];
	int i;
	for(i = 0; i < 4096; i++)
		szTmp[i] = rand() % 256;
	FILE* pFile = fopen(szFilename, "r+");
	if(!pFile)
		return false;
	FileHolder hFile(pFile);
	int nFileSize = filelength(fileno(pFile));
	while(nFileSize > 0)
	{
		fwrite(szTmp, 1, 4096, pFile);
		nFileSize -= 4096;
	}
	fflush(pFile);
	hFile.release();
	return DeleteFile(szFilename) ? true : false;
#else // WIN32
	GTEMPBUF(char, szTmp, strlen(szFilename) + 32);
	strcpy(szTmp, "shred -fun1 ");
	strcat(szTmp, szFilename);
	bool bOK = (system(szTmp) == 0);
	return bOK;
#endif // !WIN32
}
*/
/*
// static
bool GFile::shredFolder(const char* szPath)
{
	char* szOldDir = new char[300]; // use heap so deep recursion won't overflow stack
	ArrayHolder<char> hOldDir(szOldDir);
	if(!getcwd(szOldDir, 300))
		ThrowError("Failed to read current dir");

	// Recurse subdirs
	bool bOK = true;
	{
		if(chdir(szPath) != 0)
			return false;
		GDirList dl(false, false, true, false);
		while(true)
		{
			const char* szDir = dl.GetNext();
			if(!szDir)
				break;
			if(!shredFolder(szDir))
				bOK = false;
		}
	}

	// Delete files
	{
		GDirList dl(false, true, false, false);
		while(true)
		{
			const char* szFile = dl.GetNext();
			if(!szFile)
				break;
			if(!GFile::shredFile(szFile))
				bOK = false;
		}
	}

	if(chdir(szOldDir) != 0)
		ThrowError("Failed to restore the old dir");
#ifdef WIN32
	RemoveDirectory(szPath);
#else // WIN32
	rmdir(szPath);
#endif // !WIN32
	return bOK;
}
*/
// static
void GFile::parsePath(const char* szPath, struct PathData* pData)
{
#ifdef WIN32
	if(szPath[0] != '\0' && szPath[1] == ':')
		pData->dirStart = 2;
	else
		pData->dirStart = 0;
#else
	pData->dirStart = 0;
#endif // !WIN32
	int n;
	int lastSlash = -1;
	int lastPeriod = -1;
	for(n = pData->dirStart; szPath[n] != '\0'; n++)
	{
		if(szPath[n] == '/'
#ifdef WIN32
			|| szPath[n] == '\\'
#endif // WIN32
			)
			lastSlash = n;
		else if(szPath[n] == '.')
			lastPeriod = n;
	}
	pData->fileStart = lastSlash + 1;
	if(lastPeriod <= lastSlash)
		pData->extStart = n;
	else
		pData->extStart = lastPeriod;
	pData->len = n;
}

int g_tempID = 0;

// static
void GFile::tempFilename(char* pBuf)
{
#ifdef WIN32
	char tmp[512];
	GetTempPath(512, tmp);
	tmp[256] = '\0';
	size_t len = strlen(tmp);
	if(len > 0 && tmp[len - 1] != '\\' && tmp[len - 1] != '/')
		strcpy(tmp + len, "/");
#endif
	time_t t;
	std::ostringstream os;
#ifdef WIN32
	os << tmp << "tmp";
#else
	os << "/tmp/tmp";
#endif
	os << rand() << time(&t) << getpid() << (g_tempID++) << ".tmp";
	string s = os.str();
	strcpy(pBuf, s.c_str());
}

class CompressPieceComparer
{
public:
	unsigned int m_size;

	CompressPieceComparer(unsigned int size) : m_size(size)
	{
	}

	bool operator() (const unsigned char* a, const unsigned char* b) const
	{
		for(unsigned int i = 0; i < m_size; i++)
		{
			if(a[i] < b[i])
				return true;
			else if(a[i] > b[i])
				return false;
		}
		return false;
	}
};

#include <set>
#include <map>

using std::set;
using std::map;
using std::make_pair;


unsigned char* compressWorker(unsigned char* pData, unsigned int len, unsigned int* pOutNewLen, unsigned int keySize, unsigned int pieceSize, unsigned int origLen)
{
	// Build a map that counts how often all the pieces occur
	GAssert(keySize > 0); // key size must be more than 0
	GAssert(pieceSize > keySize && pieceSize <= 255); // piece size must be more than key size and less than 256
	GAssert(len <= ((unsigned int)1 << (keySize * 8 - 1))); // len too big for this keySize
	CompressPieceComparer cpc1(pieceSize);
	map<unsigned char*,unsigned char,CompressPieceComparer> pieces(cpc1); // map from offset to number of occurrences
	CompressPieceComparer cpc2(keySize);
	set<unsigned char*,CompressPieceComparer> nonKeys(cpc2); // map from offset to number of occurrences
	unsigned int offset;
	int timeout = 0;
	for(offset = 0; offset + pieceSize <= len; offset++)
	{
		if(timeout > 0)
		{
			timeout--;
			set<unsigned char*,CompressPieceComparer>::iterator it2 = nonKeys.find(pData + offset);
			if(it2 == nonKeys.end())
				nonKeys.insert(pData + offset);
		}
		else
		{
			map<unsigned char*,unsigned char,CompressPieceComparer>::iterator it = pieces.find(pData + offset);
			if(it == pieces.end())
			{
				pieces.insert(make_pair(pData + offset, (unsigned char)1));
				set<unsigned char*,CompressPieceComparer>::iterator it2 = nonKeys.find(pData + offset);
				if(it2 == nonKeys.end())
					nonKeys.insert(pData + offset);
			}
			else
			{
				if(it->second < (unsigned char)255)
					it->second++;
				timeout = pieceSize - 1;
			}
		}
	}
	for(; offset + keySize <= len; offset++)
	{
		set<unsigned char*,CompressPieceComparer>::iterator it2 = nonKeys.find(pData + offset);
		if(it2 == nonKeys.end())
			nonKeys.insert(pData + offset);
	}

	// Find all the frequently occurring pieces
	unsigned char minOccurrences = (unsigned char)2;
	while(minOccurrences < 255 && keySize + pieceSize >= minOccurrences * (pieceSize - keySize)) // while the cost is more than the savings
		minOccurrences++;
	map<unsigned char*,unsigned char*,CompressPieceComparer> fops(cpc1); // map from frequently occurring pieces to the corresponding key index
	for(map<unsigned char*,unsigned char,CompressPieceComparer>::iterator it = pieces.begin(); it != pieces.end(); it++)
	{
		if(it->second >= minOccurrences)
		{
			fops.insert(make_pair(it->first, (unsigned char*)NULL));
			if(fops.size() >= 65535)
				break;
		}
	}
	if(fops.size() == 0)
		return NULL;

	// Write the original size, key size, and piece size
	if(len < sizeof(unsigned int) + sizeof(unsigned char) + sizeof(unsigned char))
		return NULL;
	unsigned char* pOut = new unsigned char[len];
	ArrayHolder<unsigned char> hOut(pOut);
	unsigned int pos = 0;
	*(unsigned int*)(pOut + pos) = origLen;
	pos += sizeof(unsigned int);
	pOut[pos++] = (unsigned char)keySize;
	pOut[pos++] = (unsigned char)pieceSize;

	// Write the number of entries in the key-piece table
	if(pos + sizeof(unsigned short) > len)
		return NULL;
	*(unsigned short*)(pOut + pos) = (unsigned short)fops.size();
	pos += sizeof(unsigned short);

	// Write the key-piece table
	unsigned char* pPrevKey = NULL;
	for(map<unsigned char*,unsigned char*,CompressPieceComparer>::iterator it = fops.begin(); it != fops.end(); it++)
	{
		if(pos + keySize + pieceSize > len)
			return NULL;
		unsigned char* pKey = pOut + pos;
		unsigned char* pPiece = pKey + keySize;
		pos += (keySize + pieceSize);
		if(pPrevKey)
		{
			memcpy(pKey, pPrevKey, keySize);
			pPrevKey = pKey;
		}
		else
		{
			memset(pKey, 2, keySize);
			*pKey = 1;
		}
		while(true)
		{
			// Increment the key
			for(unsigned int i = 0; i < keySize; i++)
				if(++pKey[i] != 0)
					break;
			if(nonKeys.find(pKey) == nonKeys.end())
				break;
		}
		nonKeys.insert(pKey);
		pPrevKey = pKey;
		memcpy(pPiece, it->first, pieceSize);
		it->second = pKey;
	}

	// Write the payload
	for(offset = 0; offset < len; offset++)
	{
		map<unsigned char*,unsigned char*,CompressPieceComparer>::iterator it;
		if(offset + pieceSize <= len)
			it = fops.find(pData + offset);
		else
			it = fops.end();
		if(it == fops.end())
		{
			if(pos >= len)
				return NULL;
			pOut[pos++] = pData[offset];
		}
		else
		{
			if(pos + keySize > len)
				return NULL;
			memcpy(pOut + pos, it->second, keySize);
			pos += keySize;
			offset += (pieceSize - 1);
		}
	}
	if(pos >= len)
		return NULL;
	*pOutNewLen = pos;
	return hOut.release();
}

const unsigned char g_pieceSizes[] = { 18, 12, 7, 6, 5, 4, 3, 2 };
#define PIECE_SIZE_COUNT 8

// static
unsigned char* GCompressor::compress(unsigned char* pIn, unsigned int len, unsigned int* pOutNewLen)
{
	unsigned int origLen = len;

	// todo: This isn't a very intelligent way to select keySize. For example, even large text files
	// could probably use very small keys since they don't use binary characters. I think this logic
	// could be better if it actually scanned the file to see of keySize could be 1 or 2, and then
	// used this logic if both of those failed to have unique keys.
	unsigned int keySize = 1;
	while(keySize < 5 && (unsigned int)1 << (keySize * 8 - 1) < origLen)
		keySize++;

	unsigned char* pOut = NULL;
	unsigned int newLen = 0;
	for(int i = 0; i < PIECE_SIZE_COUNT; i++)
	{
		unsigned char pieceSize = g_pieceSizes[i];
		if(pieceSize <= keySize)
			break;
		unsigned char* pOldOut = pOut;
		pOut = compressWorker(pIn, len, &newLen, keySize, (unsigned int)pieceSize, origLen);
		if(pOut)
		{
			delete[] pOldOut;
			pIn = pOut + sizeof(unsigned int); // don't include the final size again
			len = newLen - sizeof(unsigned int);
		}
		else
			pOut = pOldOut;
	}
	if(!pOut)
	{
		pOut = new unsigned char[sizeof(unsigned int) + sizeof(unsigned char) + len];
		*(unsigned int*)pOut = origLen;
		*(pOut + sizeof(unsigned int)) = '\0';
		memcpy(pOut + sizeof(unsigned int) + sizeof(unsigned char), pIn, len);
		newLen = len + sizeof(unsigned int) + sizeof(unsigned char);
	}
	*pOutNewLen = newLen;
	return pOut;
}

unsigned int uncompressWorker(unsigned char* pIn, unsigned int inLen, unsigned char* pOut, unsigned int outLen)
{
	// Read the key size
	if(inLen < sizeof(unsigned char))
		ThrowError("invalid data");
	unsigned int keySize = (unsigned int)*pIn;
	pIn++;
	inLen--;

	// Read the piece size
	if(inLen < sizeof(unsigned char))
		ThrowError("invalid data");
	unsigned int pieceSize = (unsigned int)*pIn;
	pIn++;
	inLen--;
	if(keySize >= pieceSize)
		ThrowError("invalid data");

	// Read the number of table entries
	if(inLen < sizeof(unsigned short))
		ThrowError("invalid data");
	unsigned int origLen = inLen;
	unsigned short keyCount = *(unsigned short*)pIn;
	pIn += sizeof(unsigned short);
	inLen -= sizeof(unsigned short);

	// Read the table
	CompressPieceComparer cpc(keySize);
	map<unsigned char*,unsigned char*,CompressPieceComparer> table(cpc);
	for(unsigned short i = 0; i < keyCount; i++)
	{
		if(inLen < keySize + pieceSize)
			ThrowError("invalid data");
		unsigned char* pKey = pIn;
		pIn += keySize;
		inLen -= keySize;
		unsigned char* pPiece = pIn;
		pIn += pieceSize;
		inLen -= pieceSize;
		table.insert(make_pair(pKey, pPiece));
	}

	// Parse the payload
	unsigned int newLen = 0;
	while(inLen > 0)
	{
		map<unsigned char*,unsigned char*,CompressPieceComparer>::iterator it;
		if(inLen >= keySize)
			it = table.find(pIn);
		else
			it = table.end();
		if(it == table.end())
		{
			if(outLen < 1)
				ThrowError("invalid data");
			*(pOut++) = *(pIn++);
			inLen--;
			outLen--;
			newLen++;
		}
		else
		{
			if(outLen < pieceSize)
				ThrowError("invalid data");
			memcpy(pOut, it->second, pieceSize);
			pIn += keySize;
			inLen -= keySize;
			pOut += pieceSize;
			outLen -= pieceSize;
			newLen += pieceSize;
		}
	}
	if(newLen <= origLen)
		ThrowError("invalid data");
	return newLen;
}

// static
unsigned char* GCompressor::uncompress(unsigned char* pIn, unsigned int len, unsigned int* pOutUncompressedLen)
{
	if(len < sizeof(unsigned int))
		ThrowError("invalid data");
	unsigned int origLen = *(unsigned int*)pIn;
	pIn += sizeof(unsigned int);
	len -= sizeof(unsigned int);
	if(len < sizeof(unsigned char))
		ThrowError("invalid data");
	if(*pIn == '\0')
	{
		pIn++;
		len--;
		if(len != origLen)
			ThrowError("invalid data");
		unsigned char* pOut = new unsigned char[len];
		memcpy(pOut, pIn, len);
		*pOutUncompressedLen = len;
		return pOut;
	}
	unsigned char* pOut = new unsigned char[origLen];
	ArrayHolder<unsigned char> hOut(pOut);
	ArrayHolder<unsigned char> hCur(NULL);
	unsigned char* pCur = pIn;
	while(true)
	{
		unsigned int newLen = uncompressWorker(pCur, len, pOut, origLen);
		if(newLen <= len)
			ThrowError("invalid data");
		if(newLen > origLen)
			ThrowError("invalid data");
		if(newLen == origLen)
			break;
		if(pCur == pIn)
		{
			pCur = new unsigned char[origLen];
			hCur.reset(pCur);
		}
		std::swap(pCur, pOut);
		len = newLen;
	}
	*pOutUncompressedLen = origLen;
	if(hOut.get() == pOut)
		return hOut.release();
	else
		return hCur.release();
}

#ifndef NO_TEST_CODE
// static
void GCompressor::test()
{
	const char* szTest = "a one and a two and a one two three four, one two three four xxx a one yyy a one zzz a one";
	unsigned int len = strlen(szTest);
	unsigned int compressedLen;
	unsigned char* pCompressed = GCompressor::compress((unsigned char*)szTest, len, &compressedLen);
	ArrayHolder<unsigned char> hCompressed(pCompressed);
	if(compressedLen >= len)
		ThrowError("failed to compress");
	unsigned int finalLen;
	unsigned char* pFinal = GCompressor::uncompress(pCompressed, compressedLen, &finalLen);
	ArrayHolder<unsigned char> hFinal(pFinal);
	if(finalLen != len)
		ThrowError("failed to uncompress");
	if(memcmp(szTest, pFinal, len) != 0)
		ThrowError("not the same");
}
#endif
