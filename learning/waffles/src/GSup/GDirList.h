/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#ifndef __GDIRLIST_H__
#define __GDIRLIST_H__

#ifdef WIN32


/// ----------------------
///  WINDOWS VERSION
/// ----------------------
#include <sstream>
#include <string>

namespace GClasses {

class GFileFinder;
class GQueue;

/// Iterates through the files and/or folders in the current directory.
class GDirList
{
protected:
	GFileFinder* m_pFinder[256]; // it won't search more than 256 dir nests deep
	int m_nNests;
	std::ostringstream m_buffer;
	std::string m_tempBuf;
	char m_szOldDir[256];
	bool m_bReportFiles;
	bool m_bReportDirs;
	bool m_bRecurseSubDirs;
	bool m_bReportPaths;

public:
	GDirList(bool bRecurseSubDirs = true, bool bReportFiles = true, bool bReportDirs = false, bool bReportPaths = true);
	virtual ~GDirList();
	const char* GetNext();
};

} // namespace GClasses

#else // WIN32

// ----------------------
//  LINUX VERSION
// ----------------------
#include <sys/types.h>
#include <sys/dir.h>
#include <dirent.h>
#include <sstream>
#include <string>

namespace GClasses {

class GQueue;

class GDirList
{
protected:
	DIR* m_pDirs[256]; // it won't search more than 256 dir nests deep
	DIR* m_pCurDir;
	int m_nNests;
	std::ostringstream m_buffer;
	std::string m_tempBuf;
	char m_szOldDir[256];
	bool m_bReportFiles;
	bool m_bReportDirs;
	bool m_bRecurseSubDirs;
	bool m_bReportPaths;

public:
	GDirList(bool bRecurseSubDirs = true, bool bReportFiles = true, bool bReportDirs = false, bool bReportPaths = true);
	virtual ~GDirList();
	const char* GetNext();
};

} // namespace GClasses

#endif // !WIN32

#endif // __GDIRLIST_H__
