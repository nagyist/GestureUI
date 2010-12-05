/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#include "GMorph.h"
#include "../GClasses/GRand.h"
#include "../GClasses/GData.h"
#include "../GClasses/GVec.h"
#include "../GClasses/GNeighborFinder.h"
#include "../GClasses/GGraph.h"

using namespace GClasses;

GMorph::GMorph(int nInputDims, int nOutputDims, GData* pData1, GData* pData2, double dOutputToInputRatio, GRand* pRand)
{
	m_pRand = pRand;
	m_nInputDims = nInputDims;
	m_nOutputDims = nOutputDims;
	m_pData1 = pData1;
	m_pData2 = pData2;
	m_dSquaredRatio = dOutputToInputRatio * dOutputToInputRatio;
	m_pEdgeMask = new double[2 * m_nInputDims];
	m_pInputDelta = &m_pEdgeMask[m_nInputDims];
	int i;
	for(i = 0; i < m_nInputDims; i++)
		m_pEdgeMask[i] = pRand->uniform();
}

GMorph::~GMorph()
{
	delete[] m_pEdgeMask;
}

// static
void GMorph::correspondingPoints(int nInputDims, int nOutputDims, GData* pData1, GData* pData2, double dOutputToInputRatio, std::vector<int>* pPoints, GRand* pRand)
{
	GMorph morph(nInputDims, nOutputDims, pData1, pData2, dOutputToInputRatio, pRand);
	morph.findCorrespondence(pPoints);
}

bool GMorph::checkEdge(double* pOrigin, double* pTarget)
{
	int i;
	for(i = 0; i < m_nInputDims; i++)
		m_pInputDelta[i] = pTarget[i] - pOrigin[i];
	return(GVec::dotProduct(m_pInputDelta, m_pEdgeMask, m_nInputDims) > 0);
}

void GMorph::makeGraph(GGraphCut* pGraph, int nSource)
{
	double dOpeningSize = .001; // .1
	double dWallThickness = .7; // .5

	// Find the K-nearest neighbors in both data sets
	GDissimilarityMetric* pMetric = new GRowDistance();
	int nNeighbors = 2 * m_nInputDims;
	GKdTree neighbors1(m_pData1, m_nOutputDims, nNeighbors, pMetric, true);
	GKdTree neighbors2(m_pData2, m_nOutputDims, nNeighbors, pMetric, false);
	GTEMPBUF(size_t, pNeighbors, nNeighbors);
	GTEMPBUF(double, pNeighborDistances, nNeighbors);

	// Get some initial values
	GTEMPBUF(double, pWorkinGVecs, m_nInputDims * 4);
	double* pMins1 = &pWorkinGVecs[0 * m_nInputDims];
	double* pRanges1 = &pWorkinGVecs[1 * m_nInputDims];
	double* pMins2 = &pWorkinGVecs[2 * m_nInputDims];
	double* pRanges2 = &pWorkinGVecs[3 * m_nInputDims];
	int nPos1, nPos2, i, index;
	for(i = 0; i < m_nInputDims; i++)
	{
		m_pData1->minAndRange(i, &pMins1[i], &pRanges1[i]);
		m_pData2->minAndRange(i, &pMins2[i], &pRanges2[i]);
	}

	// Make the graph
	int nNode = 0;
	double* pVec1;
	double* pVec2;
	double* pVec3;
	double cost, dClosestEdge1, dClosestEdge2, dInputGap, d;
	for(nPos1 = 0; nPos1 < (int)m_pData1->rows(); nPos1++)
	{
		pVec1 = m_pData1->row(nPos1);
		for(nPos2 = 0; nPos2 < (int)m_pData2->rows(); nPos2++)
		{
			pVec2 = m_pData2->row(nPos2);
			cost = GVec::squaredDistance(pVec1, pVec2, m_nInputDims) + GVec::squaredDistance(pVec1 + m_nInputDims, pVec2 + m_nInputDims, m_nOutputDims) * m_dSquaredRatio;

			// Make an edge to non-masked neighbors in m_pData1
			neighbors1.neighbors(pNeighbors, pNeighborDistances, nPos1);
			for(i = 0; i < nNeighbors; i++)
			{
				pVec3 = m_pData1->row(pNeighbors[i]);
				if(checkEdge(pVec1, pVec3))
				{
					index = (int)(pNeighbors[i] * m_pData2->rows() + nPos2);
					pGraph->addEdge(nNode, index, (float)cost);
				}
			}

			// Make an edge to non-masked neighbors in m_pData2
			neighbors2.neighbors(pNeighbors, pNeighborDistances, nPos2);
			for(i = 0; i < nNeighbors; i++)
			{
				pVec3 = m_pData2->row(pNeighbors[i]);
				if(checkEdge(pVec2, pVec3))
				{
					index = (int)(nPos1 * m_pData2->rows() + pNeighbors[i]);
					pGraph->addEdge(nNode, index, (float)cost);
				}
			}

			// Connect to the source or sink if it's near the edge
			dClosestEdge1 = 1e200;
			dClosestEdge2 = 1e200;
			dInputGap = 0;
			for(i = 0; i < m_nInputDims; i++)
			{
				d = (pVec1[i] - pMins1[i]) / pRanges1[i];
				if(d < dClosestEdge1)
					dClosestEdge1 = d;
				d = 1.0 - (pVec2[i] - pMins2[i]) / pRanges2[i];
				if(d < dClosestEdge1)
					dClosestEdge1 = d;
				d = (pVec2[i] - pMins2[i]) / pRanges2[i];
				if(d < dClosestEdge2)
					dClosestEdge2 = d;
				d = 1.0 - (pVec1[i] - pMins1[i]) / pRanges1[i];
				if(d < dClosestEdge2)
					dClosestEdge2 = d;
				dInputGap = MAX(dInputGap, ABS((pVec1[i] - pMins1[i]) / pRanges1[i] - (pVec2[i] - pMins2[i]) / pRanges2[i]));
			}
			dInputGap -= dOpeningSize;
			dInputGap *= dWallThickness;
			if(dInputGap > dClosestEdge1)
			{
				// Connect to the source
				GAssert(dInputGap <= dClosestEdge2); // oh no, the walls overlap. This means there's no path across
				pGraph->addEdge(nNode, nSource, (float)1e30);
			}
			else if(dInputGap > dClosestEdge2)
			{
				// Connect to the sink
				pGraph->addEdge(nNode, nSource + 1, (float)1e30);
			}
			nNode++;
		}
	}
}

void GMorph::gatherPathPoints(GGraphCut* pGraph, std::vector<int>* pPoints)
{
	GGraphEdgeIterator iter(pGraph, 0); 
	bool bSource;
	double* pVec1;
	double* pVec2;
	int nNode = 0;
	bool bOutgoing;
	float fWeight;
	int nPos1, nPos2, nNeighbor, i, j;
	for(nPos1 = 0; nPos1 < (int)m_pData1->rows(); nPos1++)
	{
		for(nPos2 = 0; nPos2 < (int)m_pData2->rows(); nPos2++)
		{
			if(pGraph->doesBorderTheCut(nNode))
			{
				bSource = pGraph->isSource(nNode);
				iter.reset(nNode);
				while(iter.next(&nNeighbor, &fWeight, &bOutgoing))
				{
					if(pGraph->isSource(nNeighbor) == bSource)
						continue;
					i = nNeighbor / (int)m_pData2->rows();
					j = nNeighbor % (int)m_pData2->rows();
					if(i == nPos1)
					{
						GAssert(j != nPos2); // expected a difference
						pVec1 = m_pData2->row(nPos2);
						pVec2 = m_pData2->row(j);
					}
					else
					{
						GAssert(j == nPos2); // expected one of the values to be the same
						pVec1 = m_pData1->row(nPos1);
						pVec2 = m_pData1->row(i);
					}
					if(!checkEdge(pVec1, pVec2))
						continue;
					pPoints->push_back(nNode);
					break;
				}
			}
			nNode++;
		}
	}
}

void GMorph::findCorrespondence(std::vector<int>* pPoints)
{
	int nSource = (int)m_pData1->rows() * (int)m_pData2->rows();
	GGraphCut graph(nSource + 2);
	makeGraph(&graph, nSource);
	graph.cut(nSource, nSource + 1);
	gatherPathPoints(&graph, pPoints);
}

