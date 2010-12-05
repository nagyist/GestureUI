/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#ifndef __SELFORGANIZINGMAP_H__
#define __SELFORGANIZINGMAP_H__

#include "GTransform.h"

namespace GClasses {

class GData;
class GRand;

/// An implementation of a Kohonen map
class GSelfOrganizingMap : public GTransform
{
protected:
	int m_nMapDims;
	int m_nMapWidth;
	double m_dLearningRate;
	double m_dFocusFactor;
	GRand* m_pRand;

public:
	/// nMapDims specifies the number of dimensions for the map.
	/// nMapWidth specifies the size in one dimension of the map.
	/// (so if nMapDims is 3 and nMapWidth is 10, the map will contain 1000 (10^3) nodes.)
	GSelfOrganizingMap(int nMapDims, int nMapWidth, GRand* pRand);
	virtual ~GSelfOrganizingMap();

#ifndef NO_TEST_CODE
	/// Performs unit tests for this class. Throws an exception if there is a failure.
	static void test();
#endif

	/// Transforms pIn
	virtual GData* doit(GData* pIn);

	/// Makes the map.
	GData* makeMap(GData* pData, int nInOffset = 0);
};

} // namespace GClasses

#endif // __SELFORGANIZINGMAP_H__
