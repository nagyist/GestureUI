/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#ifndef __GTIME_H__
#define __GTIME_H__

#include <string>

namespace GClasses {

/// Provides some time-related functions
class GTime
{
public:
	/// Returns the number of seconds since "time X" with at least millisecond precision.
	/// On Linux, "time X" is the Epoch (midnight, Jan 1, 1970, GMT).
	/// On Windows, "time X" is the time when the system was started.
	static double seconds();

	/// Returns a string representation of the current time
	static const char* asciiTime(char* szBuf, int nSize, bool bGreenwichMeanTime = false);

	/// Returns a string suitable for use in a filename like "2009-12-31-23-59-59"
	static void appendTimeStampValue(std::string* pS, bool bGreenwichMeanTime);
};

} // namespace GClasses

#endif // __GTIME_H__
