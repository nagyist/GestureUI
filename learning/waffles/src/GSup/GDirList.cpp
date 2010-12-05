/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#ifdef WIN32






// ----------------------
//  WINDOWS VERSION
// ----------------------

#include "GDirList.h"
#include <winsock2.h>

using std::ostringstream;
using std::string;


namespace GClasses {


class GFileFinder : public WIN32_FIND_DATA
{
protected:
	HANDLE m_hFind;

public:
	GFileFinder()                        { m_hFind = INVALID_HANDLE_VALUE; }
	GFileFinder(LPCTSTR pFile)           { m_hFind = INVALID_HANDLE_VALUE; GetFileInfo(pFile); }
	virtual ~GFileFinder()               { Close(); }
	
	void Close()                        { if(m_hFind != INVALID_HANDLE_VALUE){FindClose(m_hFind); m_hFind = INVALID_HANDLE_VALUE;} }
	
	BOOL FindFirst(LPCTSTR pFile)       { Close(); m_hFind = FindFirstFile(pFile, this); return m_hFind != INVALID_HANDLE_VALUE; }
	BOOL FindNext()                     { return m_hFind == INVALID_HANDLE_VALUE ? FALSE : FindNextFile(m_hFind, this); }
	BOOL Find(LPCTSTR pFile);
	BOOL Find()                         { return m_hFind == INVALID_HANDLE_VALUE ? FALSE : Find(NULL); }
	
	BOOL IsValid()                      { return m_hFind != INVALID_HANDLE_VALUE; }
	BOOL IsDots();
	BOOL IsNotDots()                    { return IsDots() ? FALSE : TRUE; }
	BOOL IsFile()                       { return IsNotDots(); }
	BOOL IsFile(BOOL excludeDirs)       { return (excludeDirs && IsDir()) ? FALSE : IsNotDots(); }
	BOOL IsDir()                        { return dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY ? TRUE : FALSE; }
	
	BOOL GetFileInfo(LPCTSTR pFile)     { BOOL rval = FindFirst(pFile); Close(); return rval; }
	
	DWORD    GetAttributes()            { return dwFileAttributes; }
	FILETIME GetCreationTime()          { return ftCreationTime; }
	FILETIME GetLastAccessTime()        { return ftLastAccessTime; }
	FILETIME GetLastWriteTime()         { return ftLastWriteTime; }
	DWORD    GetFileSizeHigh()          { return nFileSizeHigh; }
	DWORD    GetFileSizeLow()           { return nFileSizeLow; }
	LPSTR    GetFileName()              { return cFileName; }
	LPSTR    GetDOSName()               { return cAlternateFileName; }
};





BOOL GFileFinder::Find(LPCTSTR pFile)
{
   if(m_hFind == INVALID_HANDLE_VALUE)
   {
      m_hFind = FindFirstFile(pFile, this);
   }
   else
   {
      if(!FindNextFile(m_hFind, this))
         Close();
   }
   return m_hFind != INVALID_HANDLE_VALUE;
}

BOOL GFileFinder::IsDots()
{
   int i = 0;
   while(cFileName[i])
   {
      if(cFileName[i++] != '.')
         return FALSE;
   }
   return TRUE;
}

GDirList::GDirList(bool bRecurseSubDirs, bool bReportFiles, bool bReportDirs, bool bReportPaths)
{
   m_bReportFiles = bReportFiles;
   m_bReportDirs = bReportDirs;
   m_bRecurseSubDirs = bRecurseSubDirs;
   m_bReportPaths = bReportPaths;
   GetCurrentDirectory(255, m_szOldDir);
   m_pFinder[0] = new GFileFinder();
   m_nNests = 1;
}

GDirList::~GDirList()
{
   for( ; m_nNests > 0; m_nNests--)
      delete(m_pFinder[m_nNests - 1]);
   SetCurrentDirectory(m_szOldDir);
}

const char* GDirList::GetNext()
{
	GFileFinder* pFinder;
	if(m_nNests > 0)
		pFinder = m_pFinder[m_nNests - 1];
	else
		return NULL;
	if(pFinder->Find("*"))
	{
		if(pFinder->IsDir())
		{
			if(pFinder->IsDots())
				return(GetNext());
			if(m_bReportDirs)
			{
				if(m_bReportPaths)
				{
					char szBuff[256];
					GetCurrentDirectory(255, szBuff);
					m_buffer.str("");
					m_buffer.clear();
					m_buffer << szBuff;
					m_buffer << "\\";
				}
				else
				{
					m_buffer.str("");
					m_buffer.clear();
				}
			}
			if(m_bRecurseSubDirs)
			{
				SetCurrentDirectory(pFinder->GetFileName());
				m_pFinder[m_nNests] = new GFileFinder();
				m_nNests++;
			}
			if(m_bReportDirs)
			{
				m_buffer << pFinder->GetFileName();
				m_tempBuf = m_buffer.str();
				m_buffer.str("");
				m_buffer.clear();
				return m_tempBuf.c_str();
			}
			else
				return GetNext();
		}
		else
		{
			if(m_bReportFiles)
			{
				if(m_bReportPaths)
				{
					char szBuff[256];
					GetCurrentDirectory(255, szBuff);
					m_buffer.str("");
					m_buffer.clear();
					m_buffer << szBuff;
					m_buffer << "\\";
				}
				else
				{
					m_buffer.str("");
					m_buffer.clear();
				}
				m_buffer << pFinder->GetFileName();
				m_tempBuf = m_buffer.str();
				m_buffer.str("");
				m_buffer.clear();
				return m_tempBuf.c_str();
			}
			else
			{
				return GetNext();
			}
		}
	}
	else
	{
		if(m_nNests > 0)
		{
			SetCurrentDirectory("..");
			m_nNests--;
			delete(m_pFinder[m_nNests]);
			return(GetNext());
		}
		else
			return NULL;
	}
}


} // namespace GClasses












#else








// ----------------------
//  LINUX VERSION
// ----------------------


#include <unistd.h>
#include "GDirList.h"
#include "../GClasses/GMacros.h"
#include <string.h>

using namespace GClasses;
using std::ostringstream;
using std::string;

////////////////////////////////////////////////////////

GDirList::GDirList(bool bRecurseSubDirs, bool bReportFiles, bool bReportDirs, bool bReportPaths)
{
	m_bReportFiles = bReportFiles;
	m_bReportDirs = bReportDirs;
	m_bRecurseSubDirs = bRecurseSubDirs;
	m_bReportPaths = bReportPaths;
	if(!getcwd(m_szOldDir, 255))
		ThrowError("failed to read cur dir");
	m_nNests = 0;
	m_pCurDir = opendir( "." );
}

GDirList::~GDirList()
{
	if(m_pCurDir)
		closedir(m_pCurDir);
	for( ; m_nNests > 0; m_nNests--)
	{
		closedir(m_pDirs[m_nNests - 1]);
		if(chdir("..") != 0)
			ThrowError("Failed to move up one dir");
	}
}

const char* GDirList::GetNext()
{
	//The current directory isn't opening
	if(m_pCurDir == NULL)
		return NULL;

	struct dirent *pDirent;
	pDirent = readdir(m_pCurDir);
	if(pDirent != NULL)
	{
		if(pDirent->d_type == DT_DIR)
		{
			//skip the . and .. directories
			if(!strcmp(pDirent->d_name, ".") || !strcmp(pDirent->d_name, ".."))
				return(GetNext());

			//We need the full path if we want to open the next directory
			if(m_bReportPaths)
			{
				char szBuff[256];
				if(!getcwd(szBuff, 255))
					ThrowError("Failed to read cur dir");
				m_buffer.str("");
				m_buffer.clear();
				m_buffer << szBuff;
				m_buffer << "/";
			}
			else
			{
				m_buffer.str("");
				m_buffer.clear();
			}

			if(m_bRecurseSubDirs)
			{ 
				//Put the current Dir object on the recursion stack, 
				//change the current dir to the new one in preparation for next query
				m_pDirs[m_nNests] = m_pCurDir;
				m_nNests++;
				if(chdir(pDirent->d_name) != 0)
					ThrowError("Failed to change dir");
				m_pCurDir = opendir(".");
			}
			if(m_bReportDirs)
			{
				m_buffer << pDirent->d_name;
				if(m_bReportPaths)
					m_buffer << "/";
				m_tempBuf = m_buffer.str();
				m_buffer.str("");
				m_buffer.clear();
				return m_tempBuf.c_str();
			}
			else
				return GetNext();
		}
		else
		{
			if(m_bReportFiles)
			{
				if(m_bReportPaths)
				{
					char szBuff[256];
					if(!getcwd(szBuff, 255))
						ThrowError("Failed to read cur dir");
					m_buffer.str("");
					m_buffer.clear();
					m_buffer << szBuff;
					m_buffer << "/";
				}
				else
				{
					m_buffer.str("");
					m_buffer.clear();
				}

				m_buffer << pDirent->d_name;
				m_tempBuf = m_buffer.str();
				m_buffer.str("");
				m_buffer.clear();
				return m_tempBuf.c_str();
			}
			else
			{
				return GetNext();
			}
		}
	}
	else
	{
		//In here, there are no more files in the current directory
		//Step out of the current nest, recurse up
		if(m_nNests > 0)
		{
			if(chdir("..") != 0)
				ThrowError("Failed to move up one dir");
			closedir(m_pCurDir);
			m_pCurDir = m_pDirs[m_nNests - 1];
			m_nNests--;
			return(GetNext());
		}
		else //all done! No more files.
		{
			closedir(m_pCurDir);
			m_pCurDir = NULL;
			return NULL;
		}
	}
}




#endif // !WIN32
