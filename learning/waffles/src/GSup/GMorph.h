/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#ifndef __GMORPH_H__
#define __GMORPH_H__

#include <vector>

namespace GClasses {

class GArffRelation;
class GData;
class GGraphCut;
class GRand;


/// This is an implementation of DP-warping, which interpolates between two functions
class GMorph
{
protected:
	int m_nInputDims;
	int m_nOutputDims;
	GData* m_pData1;
	GData* m_pData2;
	GData* m_pCorrespondence;
	double m_dSquaredRatio;
	double* m_pEdgeMask;
	double* m_pInputDelta;
	GRand* m_pRand;

	/// When dOutputToInputRatio is 0, it just does linear interpolation of the output values. When
	/// dOutputToInputRatio is very large, it does linear interpolation of the input values. When
	/// dOutputToInputRatio is 1, it weights inputs and outputs equally.
	GMorph(int nInputDims, int nOutputDims, GData* pData1, GData* pData2, double dOutputToInputRatio, GRand* pRand);
public:

	~GMorph();

	/// Produce a set of points that correspond between the two data sets.
	/// Each point p is encoded such that p = i * m_pData2->GetSize() + j,
	/// where i is the index in m_pData1 and j is the index in m_pData2.
	/// So i = p / m_pData2->GetSize(), and j = p % m_pData2->GetSize().
	static void correspondingPoints(int nInputDims, int nOutputDims, GData* pData1, GData* pData2, double dOutputToInputRatio, std::vector<int>* pPoints, GRand* pRand);

protected:
	/// Returns true iff the dot product of (pTarget - pOrigin) with
	/// m_pEdgeMask is > 0. (This effectively gives a direction to
	/// every edge. We can use this direction to make sure we only put
	/// the edge in the graph once, and to determine which node is
	/// responsible for the edge's capacity value.)
	bool checkEdge(double* pOrigin, double* pTarget);

	void makeGraph(GGraphCut* pGraph, int nSource);

	void gatherPathPoints(GGraphCut* pGraph, std::vector<int>* pPoints);

	void findCorrespondence(std::vector<int>* pPoints);
};

} // namespace GClasses

#endif // __GMORPH_H__
