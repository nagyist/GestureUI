/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#include "GNeighborFinder.h"
#include "GVec.h"
#include "GRand.h"
#include "GPlot.h"
#include <stdlib.h>
#include <vector>
#include <queue>
#include <set>
#include "GOptimizer.h"
#include "GHillClimber.h"
#include <string.h>
#include "GGraph.h"
#include "GBitTable.h"
#include <deque>
#include "GTwt.h"
#include "GKNN.h"
#include "GTransform.h"
#include <sstream>
#include <string>
#include <iostream>
#include "GNeuralNet.h"

namespace GClasses {

using std::cout;
using std::vector;
using std::priority_queue;
using std::set;
using std::deque;
using std::make_pair;
using std::pair;
using std::string;

GDissimilarityMetric::GDissimilarityMetric(GTwtNode* pNode)
{
	m_pRelation = GRelation::fromTwt(pNode->field("relation"));
}

GTwtNode* GDissimilarityMetric::baseTwtNode(GTwtDoc* pDoc, const char* szClassName)
{
	GTwtNode* pNode = pDoc->newObj();
	pNode->addField(pDoc, "class", pDoc->newString(szClassName));
	pNode->addField(pDoc, "relation", m_pRelation->toTwt(pDoc));
	return pNode;
}

// static
GDissimilarityMetric* GDissimilarityMetric::fromTwt(GTwtNode* pNode)
{
	const char* szClass = pNode->field("class")->asString();
	if(strcmp(szClass, "GRowDistanceScaled") == 0)
		return new GRowDistanceScaled(pNode);
	if(strcmp(szClass, "GRowDistance") == 0)
		return new GRowDistance(pNode);
	if(strcmp(szClass, "GMinkowskiDistance") == 0)
		return new GMinkowskiDistance(pNode);
	ThrowError("Unrecognized class: ", szClass);
	return NULL;
}

// --------------------------------------------------------------------

GRowDistance::GRowDistance(GTwtNode* pNode)
: GDissimilarityMetric(pNode)
{
}

// virtual
GTwtNode* GRowDistance::toTwt(GTwtDoc* pDoc)
{
	GTwtNode* pNode = baseTwtNode(pDoc, "GRowDistance");
	return pNode;
}

// virtual
void GRowDistance::init(sp_relation& pRelation)
{
	m_pRelation = pRelation;
}

// virtual
double GRowDistance::dissimilarity(const double* pA, const double* pB)
{
	double sum = 0;
	int count = m_pRelation->size();
	double d;
	for(int i = 0; i < count; i++)
	{
		if(m_pRelation->valueCount(i) == 0)
			d = *pB - *pA;
		else
			d = ((int)*pB == (int)*pA ? 0 : 1);
		pA++;
		pB++;
		sum += (d * d);
	}
	return sum;
}

// --------------------------------------------------------------------

GRowDistanceScaled::GRowDistanceScaled(GTwtNode* pNode)
: GDissimilarityMetric(pNode)
{
	GTwtNode* pScaleFactors = pNode->field("scaleFactors");
	int dims = m_pRelation->size();
	if((int)pScaleFactors->itemCount() != dims)
		ThrowError("wrong number of scale factors");
	m_pScaleFactors = new double[dims];
	for(int i = 0; i < dims; i++)
		m_pScaleFactors[i] = pScaleFactors->item(i)->asDouble();
}

// virtual
GTwtNode* GRowDistanceScaled::toTwt(GTwtDoc* pDoc)
{
	GTwtNode* pNode = baseTwtNode(pDoc, "GRowDistance");
	int dims = m_pRelation->size();
	GTwtNode* pScaleFactors = pNode->addField(pDoc, "scaleFactors", pDoc->newList(dims));
	for(int i = 0; i < dims; i++)
		pScaleFactors->setItem(i, pDoc->newDouble(m_pScaleFactors[i]));
	return pNode;
}

// virtual
void GRowDistanceScaled::init(sp_relation& pRelation)
{
	m_pRelation = pRelation;
	delete[] m_pScaleFactors;
	m_pScaleFactors = new double[pRelation->size()];
	GVec::setAll(m_pScaleFactors, 1.0, pRelation->size());
}

// virtual
double GRowDistanceScaled::dissimilarity(const double* pA, const double* pB)
{
	double sum = 0;
	int count = m_pRelation->size();
	double d;
	const double* pSF = m_pScaleFactors;
	for(int i = 0; i < count; i++)
	{
		if(m_pRelation->valueCount(i) == 0)
			d = (*pB - *pA) * (*pSF);
		else
			d = ((int)*pB == (int)*pA ? 0 : *pSF);
		pA++;
		pB++;
		pSF++;
		sum += (d * d);
	}
	return sum;
}

// --------------------------------------------------------------------

GMinkowskiDistance::GMinkowskiDistance(GTwtNode* pNode)
: GDissimilarityMetric(pNode)
{
}

// virtual
GTwtNode* GMinkowskiDistance::toTwt(GTwtDoc* pDoc)
{
	GTwtNode* pNode = baseTwtNode(pDoc, "GMinkowskiDistance");
	return pNode;
}

// virtual
void GMinkowskiDistance::init(sp_relation& pRelation)
{
	if(!pRelation->areContinuous(0, pRelation->size()))
		ThrowError("Only continuous attributes are supported");
	m_pRelation = pRelation;
}

// virtual
double GMinkowskiDistance::dissimilarity(const double* pA, const double* pB)
{
	return GVec::minkowskiDistance(m_norm, pA, pB, m_pRelation->size());
}

// --------------------------------------------------------------------

void GNeighborFinder_InsertionSortNeighbors(int neighborCount, size_t* pNeighbors, double* pDistances)
{
	size_t tt;
	double t;
	for(int i = 1; i < neighborCount; i++)
	{
		for(int j = i; j > 0; j--)
		{
			if(pNeighbors[j] == INVALID_INDEX)
				break;
			if(pNeighbors[j - 1] != INVALID_INDEX && pDistances[j] >= pDistances[j - 1])
				break;

			// Swap
			tt = pNeighbors[j - 1];
			pNeighbors[j - 1] = pNeighbors[j];
			pNeighbors[j] = tt;
			t = pDistances[j - 1];
			pDistances[j - 1] = pDistances[j];
			pDistances[j] = t;
		}
	}
}

void GNeighborFinder::sortNeighbors(int neighborCount, size_t* pNeighbors, double* pDistances)
{
	// Use insertion sort if the list is small
	if(neighborCount < 7)
	{
		GNeighborFinder_InsertionSortNeighbors(neighborCount, pNeighbors, pDistances);
		return;
	}
	double t;
	size_t tt;
	int beg = 0;
	int end = neighborCount - 1;

	// Pick a pivot (using the median of 3 technique)
	double pivA = pDistances[0];
	double pivB = pDistances[neighborCount / 2];
	double pivC = pDistances[neighborCount - 1];
	double pivot;
	if(pivA < pivB)
	{
		if(pivB < pivC)
			pivot = pivB;
		else if(pivA < pivC)
			pivot = pivC;
		else
			pivot = pivA;
	}
	else
	{
		if(pivA < pivC)
			pivot = pivA;
		else if(pivB < pivC)
			pivot = pivC;
		else
			pivot = pivB;
	}

	// Do Quick Sort
	while(true)
	{
		while(beg < end && pNeighbors[beg] != INVALID_INDEX && pDistances[beg] < pivot)
			beg++;
		while(end > beg && (pNeighbors[end] == INVALID_INDEX || pDistances[end] > pivot))
			end--;
		if(beg >= end)
			break;
		t = pDistances[beg];
		pDistances[beg] = pDistances[end];
		pDistances[end] = t;
		tt = pNeighbors[beg];
		pNeighbors[beg] = pNeighbors[end];
		pNeighbors[end] = tt;
		beg++;
		end--;
	}

	// Recurse
	if(pNeighbors[beg] != INVALID_INDEX && pDistances[beg] < pivot)
		beg++;
	else if(beg == 0) // This could happen if they're all -1 (bad neighbors)
	{
		GNeighborFinder_InsertionSortNeighbors(neighborCount, pNeighbors, pDistances);
		return;
	}
	GNeighborFinder::sortNeighbors(beg, pNeighbors, pDistances);
	GNeighborFinder::sortNeighbors(neighborCount - beg, pNeighbors + beg, pDistances + beg);
}

void GNeighborFinder::sortNeighbors(size_t* pNeighbors, double* pDistances)
{
	GNeighborFinder::sortNeighbors(m_neighborCount, pNeighbors, pDistances);
}









GNeighborFinderCacheWrapper::GNeighborFinderCacheWrapper(GNeighborFinder* pNF, bool own)
: GNeighborFinder(pNF->data(), pNF->neighborCount()), m_pNF(pNF), m_own(own)
{
	m_pCache = new size_t[m_pData->rows() * m_neighborCount];
	m_pDissims = new double[m_pData->rows() * m_neighborCount];
	for(size_t i = 0; i < m_pData->rows(); i++)
		m_pCache[i * m_neighborCount] = m_pData->rows();
}

// virtual
GNeighborFinderCacheWrapper::~GNeighborFinderCacheWrapper()
{
	delete[] m_pCache;
	delete[] m_pDissims;
	if(m_own)
		delete(m_pNF);
}

// virtual
void GNeighborFinderCacheWrapper::neighbors(size_t* pOutNeighbors, size_t index)
{
	size_t* pCache = m_pCache + m_neighborCount * index;
	if(*pCache == m_pData->rows())
	{
		double* pDissims = m_pDissims + m_neighborCount * index;
		((GNeighborFinder*)m_pNF)->neighbors(pCache, pDissims, index);
	}
	memcpy(pOutNeighbors, pCache, sizeof(size_t) * m_neighborCount);
}

// virtual
void GNeighborFinderCacheWrapper::neighbors(size_t* pOutNeighbors, double* pOutDistances, size_t index)
{
	size_t* pCache = m_pCache + m_neighborCount * index;
	double* pDissims = m_pDissims + m_neighborCount * index;
	if(*pCache == m_pData->rows())
		((GNeighborFinder*)m_pNF)->neighbors(pCache, pDissims, index);
	memcpy(pOutNeighbors, pCache, sizeof(size_t) * m_neighborCount);
	memcpy(pOutDistances, pDissims, sizeof(double) * m_neighborCount);
}

void GNeighborFinderCacheWrapper::fillCache()
{
	size_t rowCount = m_pData->rows();
	size_t* pCache = m_pCache;
	double* pDissims = m_pDissims;
	for(size_t i = 0; i < rowCount; i++)
	{
		if(*pCache == m_pData->rows())
			((GNeighborFinder*)m_pNF)->neighbors(pCache, pDissims, i);
		pCache += m_neighborCount;
		pDissims += m_neighborCount;
	}
}

void GNeighborFinderCacheWrapper::fillDistances(GDissimilarityMetric* pMetric)
{
	pMetric->init(m_pData->relation());
	double* pDissim = m_pDissims;
	size_t* pHood = m_pCache;
	for(size_t i = 0; i < m_pData->rows(); i++)
	{
		double* pA = m_pData->row(i);
		for(int j = 0; j < m_neighborCount; j++)
		{
			double* pB = m_pData->row(pHood[j]);
			*pDissim = pMetric->dissimilarity(pA, pB);
			pDissim++;
		}
		pHood += m_neighborCount;
	}
}

size_t GNeighborFinderCacheWrapper::cutShortcuts(int cycleLen)
{
	GCycleCut cc(m_pCache, m_pData, m_neighborCount);
	cc.setCycleThreshold(cycleLen);
	return cc.cut();
}

void GNeighborFinderCacheWrapper::patchMissingSpots(GRand* pRand)
{
	size_t rowCount = m_pData->rows();
	size_t* pCache = m_pCache;
	double* pDissims = m_pDissims;
	for(size_t i = 0; i < rowCount; i++)
	{
		if(*pCache == m_pData->rows())
			ThrowError("cache not filled out");
		for(int j = 0; j < m_neighborCount; j++)
		{
			if(pCache[j] >= m_pData->rows())
			{
				int k = (int)pRand->next(m_neighborCount);
				int l;
				for(l = k; l < m_neighborCount; l++)
				{
					if(pCache[l] < m_pData->rows())
						break;
				}
				if(l >= m_neighborCount)
				{
					for(l = 0; l < k; l++)
					{
						if(pCache[l] < m_pData->rows())
							break;
					}
				}
				if(pCache[l] >= m_pData->rows())
					ThrowError("row has zero valid neighbors");
				if(pDissims)
					pDissims[j] = pDissims[l];
				pCache[j] = pCache[l];
			}
		}
		pCache += m_neighborCount;
		pDissims += m_neighborCount;
	}
}

void GNeighborFinderCacheWrapper::normalizeDistances()
{
	size_t rowCount = m_pData->rows();
	size_t* pCache = m_pCache;
	double* pDissims = m_pDissims;
	double total = 0.0;
	for(size_t i = 0; i < rowCount; i++)
	{
		if(*pCache == m_pData->rows())
			ThrowError("cache not filled out");
		for(int j = 0; j < m_neighborCount; j++)
		{
			pDissims[j] = sqrt(pDissims[j]);
			total += pDissims[j];
		}
		pCache += m_neighborCount;
		pDissims += m_neighborCount;
	}
	pDissims = m_pDissims;
	for(size_t i = 0; i < rowCount; i++)
	{
		double s = 0;
		for(int j = 0; j < m_neighborCount; j++)
			s += pDissims[j];
		s = 1.0 / s;
		for(int j = 0; j < m_neighborCount; j++)
			pDissims[j] *= s;
		pDissims += m_neighborCount;
	}
	total /= rowCount;
	pDissims = m_pDissims;
	for(size_t i = 0; i < rowCount; i++)
	{
		for(int j = 0; j < m_neighborCount; j++)
		{
			double d = pDissims[j] * total;
			pDissims[j] = (d * d);
		}
		pDissims += m_neighborCount;
	}
}

bool GNeighborFinderCacheWrapper::isConnected()
{
	GBitTable bt(m_pData->rows());
	deque<size_t> q;
	bt.set(0);
	q.push_back(0);
	size_t count = 1;
	while(q.size() > 0)
	{
		size_t n = q.front();
		q.pop_front();
		size_t* pHood = m_pCache + m_neighborCount * n;
		for(size_t i = 0; i < (size_t)m_neighborCount; i++)
		{
			size_t neigh = *(pHood++);
			if(neigh < m_pData->rows())
			{
				if(!bt.bit(neigh))
				{
					bt.set(neigh);
					count++;
					if(count >= m_pData->rows())
						return true;
					q.push_back(neigh);
				}
			}
		}
	}
	return false;
}

// --------------------------------------------------------------------

// This helper class keeps neighbors sorted as a binary heap, such that the most dissimilar
// of the k-current-neighbors is always at the front of the heap.
class GClosestNeighborFindingHelper
{
protected:
	int m_found;
	int m_neighbors;
	size_t* m_pNeighbors;
	double* m_pDistances;

public:
	GClosestNeighborFindingHelper(int neighbors, size_t* pNeighbors, double* pDistances)
	 : m_found(0), m_neighbors(neighbors), m_pNeighbors(pNeighbors), m_pDistances(pDistances)
	{
		int i;
		for(i = 0; i < m_neighbors; i++)
		{
			m_pNeighbors[i] = -1;
			m_pDistances[i] = 1e308;
		}
	}

	~GClosestNeighborFindingHelper()
	{
	}

	// Adds a point to the set of current neighbors if it is closer than the
	// most dissimilar of the k-current-neighbors
	void TryPoint(size_t index, double distance)
	{
		double* pHeapDist = m_pDistances - 1;
		size_t* pHeapNeigh = m_pNeighbors - 1;
		int heapPos;
		if(m_found < m_neighbors)
			heapPos = ++m_found;
		else
		{
			// Compare with the front of the heap, which holds the most dissimilar of the k-current-neighbors
			if(distance >= m_pDistances[0])
				return;

			// Release the most dissimilar of the k-current neighbors
			heapPos = 1;
			while(2 * heapPos <= m_neighbors)
			{
				if(2 * heapPos == m_neighbors || pHeapDist[2 * heapPos] > pHeapDist[2 * heapPos + 1])
				{
					pHeapDist[heapPos] = pHeapDist[2 * heapPos];
					pHeapNeigh[heapPos] = pHeapNeigh[2 * heapPos];
					heapPos = 2 * heapPos;
				}
				else
				{
					pHeapDist[heapPos] = pHeapDist[2 * heapPos + 1];
					pHeapNeigh[heapPos] = pHeapNeigh[2 * heapPos + 1];
					heapPos = 2 * heapPos + 1;
				}
			}
		}

		// Insert into heap
		pHeapDist[heapPos] = distance;
		pHeapNeigh[heapPos] = index;
		while(heapPos > 1 && pHeapDist[heapPos / 2] < pHeapDist[heapPos])
		{
			std::swap(pHeapDist[heapPos / 2], pHeapDist[heapPos]);
			std::swap(pHeapNeigh[heapPos / 2], pHeapNeigh[heapPos]);
			heapPos /= 2;
		}
	}

	double GetWorstDist()
	{
		return m_found >= m_neighbors ? m_pDistances[0] : 1e308;
	}

#ifndef NO_TEST_CODE
#	define TEST_NEIGHBOR_COUNT 33
	static void test()
	{
		size_t neighbors[TEST_NEIGHBOR_COUNT];
		double distances[TEST_NEIGHBOR_COUNT];
		GData values1(1);
		GData values2(1);
		GClosestNeighborFindingHelper ob(TEST_NEIGHBOR_COUNT, neighbors, distances);
		GRand prng(0);
		for(size_t i = 0; i < 300; i++)
		{
			double d = prng.uniform();
			ob.TryPoint(i, d);
			values1.newRow()[0] = d;
			values1.sort(0);
			values2.flush();
			for(size_t j = 0; j < MIN((size_t)TEST_NEIGHBOR_COUNT, values1.rows()); j++)
				values2.newRow()[0] = distances[j];
			values2.sort(0);
			for(size_t j = 0; j < MIN((size_t)TEST_NEIGHBOR_COUNT, values1.rows()); j++)
			{
				if(ABS(values1[j][0] - values2[j][0]) > 1e-12)
					ThrowError("something is wrong");
			}
		}
	}
#endif
};

// --------------------------------------------------------------------------------

GNeighborFinderGeneralizing::GNeighborFinderGeneralizing(GData* pData, int labelDims, int neighborCount, GDissimilarityMetric* pMetric, bool ownMetric)
: GNeighborFinder(pData, neighborCount), m_pMetric(pMetric), m_ownMetric(ownMetric)
{
	if(!m_pMetric)
	{
		m_pMetric = new GRowDistance();
		m_ownMetric = true;
	}
	sp_relation pRel;
	if(labelDims > 0)
	{
		GMixedRelation* pMixedRel = new GMixedRelation();
		pMixedRel->addAttrs(pData->relation().get(), 0, pData->cols() - labelDims);
		pRel = pMixedRel;
	}
	else
		pRel = pData->relation();
	m_pMetric->init(pRel);
}

// virtual
GNeighborFinderGeneralizing::~GNeighborFinderGeneralizing()
{
	if(m_ownMetric)
		delete(m_pMetric);
}

// --------------------------------------------------------------------------------

GBruteForceNeighborFinder::GBruteForceNeighborFinder(GData* pData, int labelDims, int neighborCount, GDissimilarityMetric* pMetric, bool ownMetric)
: GNeighborFinderGeneralizing(pData, labelDims, neighborCount, pMetric, ownMetric)
{
}

GBruteForceNeighborFinder::~GBruteForceNeighborFinder()
{
}

size_t GBruteForceNeighborFinder::addVector(double* pVector)
{
	size_t index = m_pData->rows();
	m_pData->takeRow(pVector);
	return index;
}

double* GBruteForceNeighborFinder::releaseVector(size_t nIndex)
{
	return m_pData->releaseRow(nIndex);
}

// virtual
void GBruteForceNeighborFinder::reoptimize()
{
}

// virtual
void GBruteForceNeighborFinder::neighbors(size_t* pOutNeighbors, size_t index)
{
	GTEMPBUF(double, distances, m_neighborCount);
	neighbors(pOutNeighbors, distances, index);
}

// virtual
void GBruteForceNeighborFinder::neighbors(size_t* pOutNeighbors, double* pOutDistances, size_t index)
{
	GClosestNeighborFindingHelper helper(m_neighborCount, pOutNeighbors, pOutDistances);
	double* pCand;
	double dist;
	double* pInputVector = m_pData->row(index);
	for(size_t i = 0; i < m_pData->rows(); i++)
	{
		if(i == index)
			continue;
		pCand = m_pData->row(i);
		dist = m_pMetric->dissimilarity(pInputVector, pCand);
		helper.TryPoint(i, dist);
	}
}

// virtual
void GBruteForceNeighborFinder::neighbors(size_t* pOutNeighbors, double* pOutDistances, const double* pInputVector)
{
	GClosestNeighborFindingHelper helper(m_neighborCount, pOutNeighbors, pOutDistances);
	double* pCand;
	double dist;
	for(size_t i = 0; i < m_pData->rows(); i++)
	{
		pCand = m_pData->row(i);
		dist = m_pMetric->dissimilarity(pInputVector, pCand);
		helper.TryPoint(i, dist);
	}
}

// --------------------------------------------------------------------------------

class GKdNode
{
protected:
	double m_minDist;
	double* m_pOffset;
	int m_dims;

public:
	GKdNode(int dims)
	{
		m_dims = dims;
		m_pOffset = new double[dims];
		GVec::setAll(m_pOffset, 0.0, dims);
		m_minDist = 0;
	}

	virtual ~GKdNode()
	{
		delete[] m_pOffset;
	}

	virtual bool IsLeaf() = 0;

	// Builds an array of all the indexes in all of the leaf nodes that descend from me
	virtual int Gather(size_t* pOutIndexes) = 0;

	virtual void Insert(GKdTree* pTree, size_t index, double* pRow) = 0;

	virtual void Remove(GKdTree* pTree, size_t index, double* pRow) = 0;

	virtual void Rename(GKdTree* pTree, size_t oldIndex, size_t newIndex, double* pRow) = 0;

	double GetMinDist()
	{
		return m_minDist;
	}

	int GetDims()
	{
		return m_dims;
	}

	void CopyOffset(GKdNode* pParent)
	{
		GVec::copy(m_pOffset, pParent->m_pOffset, m_dims);
		m_minDist = pParent->m_minDist;
	}

	void AdjustOffset(int attr, double offset, const double* m_pScaleFactors)
	{
		if(offset > m_pOffset[attr])
		{
			if(m_pScaleFactors)
			{
				m_minDist -= (m_pOffset[attr] * m_pOffset[attr] * m_pScaleFactors[attr] * m_pScaleFactors[attr]);
				m_pOffset[attr] = offset;
				m_minDist += (m_pOffset[attr] * m_pOffset[attr] * m_pScaleFactors[attr] * m_pScaleFactors[attr]);
			}
			else
			{
				m_minDist -= (m_pOffset[attr] * m_pOffset[attr]);
				m_pOffset[attr] = offset;
				m_minDist += (m_pOffset[attr] * m_pOffset[attr]);
			}
		}
	}
};


class GKdInteriorNode : public GKdNode
{
protected:
	GKdNode* m_pLess;
	GKdNode* m_pGreaterOrEqual;
	int m_size;
	int m_attr;
	double m_pivot;

public:
	int m_timeLeft;

	GKdInteriorNode(int dims, GKdNode* pLess, GKdNode* pGreaterOrEqual, int size, int attr, double pivot)
	 : GKdNode(dims), m_pLess(pLess), m_pGreaterOrEqual(pGreaterOrEqual), m_size(size), m_attr(attr), m_pivot(pivot)
	{
		m_timeLeft = (int)MIN((double)0x7fffffff, ((double)size * size) / 36 + 6);
	}

	virtual ~GKdInteriorNode()
	{
		delete(m_pLess);
		delete(m_pGreaterOrEqual);
	}

	virtual bool IsLeaf() { return false; }

	GKdNode* Rebuild(GKdTree* pTree)
	{
		size_t* pIndexes = new size_t[m_size];
		ArrayHolder<size_t> hIndexes(pIndexes);
		int used = Gather(pIndexes);
		GAssert(used == m_size); // m_size is wrong. This may corrupt memory.
		return pTree->buildTree(used, pIndexes);
	}

	virtual void Insert(GKdTree* pTree, size_t index, double* pRow)
	{
		m_timeLeft--;
		if(pTree->isGreaterOrEqual(pRow, m_attr, m_pivot))
		{
			m_pGreaterOrEqual->Insert(pTree, index, pRow);
			m_size++;
			if(!m_pGreaterOrEqual->IsLeaf() && ((GKdInteriorNode*)m_pGreaterOrEqual)->m_timeLeft <= 0 && m_timeLeft >= m_size / 4)
			{
				GKdNode* pNewNode = ((GKdInteriorNode*)m_pGreaterOrEqual)->Rebuild(pTree);
				delete(m_pGreaterOrEqual);
				m_pGreaterOrEqual = pNewNode;
			}
		}
		else
		{
			m_pLess->Insert(pTree, index, pRow);
			m_size++;
			if(!m_pLess->IsLeaf() && ((GKdInteriorNode*)m_pLess)->m_timeLeft <= 0 && m_timeLeft >= m_size / 4)
			{
				GKdNode* pNewNode = ((GKdInteriorNode*)m_pLess)->Rebuild(pTree);
				delete(m_pLess);
				m_pLess = pNewNode;
			}
		}
	}

	virtual void Remove(GKdTree* pTree, size_t index, double* pRow)
	{
		m_timeLeft--;
		if(pTree->isGreaterOrEqual(pRow, m_attr, m_pivot))
			m_pGreaterOrEqual->Remove(pTree, index, pRow);
		else
			m_pLess->Remove(pTree, index, pRow);
		m_size--;
	}

	virtual void Rename(GKdTree* pTree, size_t oldIndex, size_t newIndex, double* pRow)
	{
		if(pTree->isGreaterOrEqual(pRow, m_attr, m_pivot))
			m_pGreaterOrEqual->Rename(pTree, oldIndex, newIndex, pRow);
		else
			m_pLess->Rename(pTree, oldIndex, newIndex, pRow);
	}

	GKdNode* GetLess() { return m_pLess; }
	GKdNode* GetGreaterOrEqual() { return m_pGreaterOrEqual; }
	int GetSize() { return m_size; }

	void GetDivision(int* pAttr, double* pPivot)
	{
		*pAttr = m_attr;
		*pPivot = m_pivot;
	}

	int Gather(size_t* pOutIndexes)
	{
		int n = m_pLess->Gather(pOutIndexes);
		return m_pGreaterOrEqual->Gather(pOutIndexes + n) + n;
	}
};


class GKdLeafNode : public GKdNode
{
protected:
	vector<size_t> m_indexes;

public:
	GKdLeafNode(int count, size_t* pIndexes, int dims, int maxLeafSize)
	 : GKdNode(dims)
	{
		m_indexes.reserve(MAX(count, maxLeafSize));
		int i;
		for(i = 0; i < count; i++)
			m_indexes.push_back(pIndexes[i]);
	}

	virtual ~GKdLeafNode()
	{
	}

	virtual bool IsLeaf() { return true; }

	int GetSize()
	{
		return (int)m_indexes.size();
	}

	virtual void Insert(GKdTree* pTree, size_t index, double* pRow)
	{
		m_indexes.push_back(index);
	}

	virtual void Remove(GKdTree* pTree, size_t index, double* pRow)
	{
		size_t count = m_indexes.size();
		for(size_t i = 0; i < count; i++)
		{
			if(m_indexes[i] == index)
			{
				m_indexes[i] = m_indexes[count - 1];
				m_indexes.pop_back();
				return;
			}
		}
		GAssert(false); // failed to find index. Did the row change?
	}

	virtual void Rename(GKdTree* pTree, size_t oldIndex, size_t newIndex, double* pRow)
	{
		size_t count = m_indexes.size();
		for(size_t i = 0; i < count; i++)
		{
			if(m_indexes[i] == oldIndex)
			{
				m_indexes[i] = newIndex;
				return;
			}
		}
		GAssert(false); // failed to find index. Did the row change?
	}

	vector<size_t>* GetIndexes() { return &m_indexes; }

	int Gather(size_t* pOutIndexes)
	{
		for(vector<size_t>::iterator i = m_indexes.begin(); i < m_indexes.end(); i++)
		{
			*pOutIndexes = *i;
			pOutIndexes++;
		}
		return (int)m_indexes.size();
	}
};


// --------------------------------------------------------------------------------------------------------

GKdTree::GKdTree(GData* pData, int labelDims, int neighborCount, GDissimilarityMetric* pMetric, bool ownMetric)
: GNeighborFinderGeneralizing(pData, labelDims, neighborCount, pMetric, ownMetric)
{
	m_maxLeafSize = 6;
	size_t count = pData->rows();
	GTEMPBUF(size_t, tmp, count);
	for(size_t i = 0; i < count; i++)
		tmp[i] = i;
	m_pRoot = buildTree(count, tmp);
}

// virtual
GKdTree::~GKdTree()
{
	delete(m_pRoot);
}

void GKdTree::computePivotAndGoodness(size_t count, size_t* pIndexes, int attr, double* pOutPivot, double* pOutGoodness)
{
	int valueCount = m_pMetric->relation()->valueCount(attr);
	if(valueCount > 0)
	{
		// Count the ocurrences of each value
		double* pPat;
		GTEMPBUF(int, counts, valueCount);
		memset(counts, '\0', sizeof(int) * valueCount);
		for(size_t i = 0; i < count; i++)
		{
			pPat = m_pData->row(pIndexes[i]);
			if((int)pPat[attr] >= 0)
			{
				GAssert((unsigned int)pPat[attr] < (unsigned int)valueCount); // out of range
				if((unsigned int)pPat[attr] < (unsigned int)valueCount)
					counts[(int)pPat[attr]]++;
			}
		}

		// Total up the entropy
		int max = 0;
		int maxcount = -1;
		double entropy = 0;
		double ratio;
		for(int i = 0; i < valueCount; i++)
		{
			if(counts[i] > maxcount)
			{
				maxcount = counts[i];
				max = i;
			}
			if(counts[i] > 0)
			{
				ratio = (double)counts[i] / count;
				entropy -= ratio * log(ratio);
			}
		}
		const double* pScaleFactors = m_pMetric->scaleFactors();
		if(pScaleFactors)
			entropy *= (pScaleFactors[attr] * pScaleFactors[attr]);

		*pOutPivot = max;
		*pOutGoodness = entropy;
	}
	else
	{
		// Compute the mean
		double mean = 0;
		double* pPat;
		for(size_t i = 0; i < count; i++)
		{
			pPat = m_pData->row(pIndexes[i]);
			mean += pPat[attr];
		}
		mean /= count;

		// Compute the scaled variance
		double var = 0;
		double d;
		const double* pScaleFactors = m_pMetric->scaleFactors();
		if(pScaleFactors)
		{
			for(size_t i = 0; i < count; i++)
			{
				pPat = m_pData->row(pIndexes[i]);
				d = (pPat[attr] - mean) * pScaleFactors[attr];
				var += (d * d);
			}
		}
		else
		{
			for(size_t i = 0; i < count; i++)
			{
				pPat = m_pData->row(pIndexes[i]);
				d = (pPat[attr] - mean);
				var += (d * d);
			}
		}
		var /= count; // (the biassed estimator of variance is better for this purpose)

		*pOutPivot = mean;
		*pOutGoodness = var;
	}
}

size_t GKdTree::splitIndexes(size_t count, size_t* pIndexes, int attr, double pivot)
{
	double* pPat;
	size_t t;
	size_t beg = 0;
	size_t end = count - 1;
	if(m_pMetric->relation()->valueCount(attr) == 0)
	{
		while(end >= beg && end < count)
		{
			pPat = m_pData->row(pIndexes[beg]);
			if(pPat[attr] >= pivot)
			{
				t = pIndexes[beg];
				pIndexes[beg] = pIndexes[end];
				pIndexes[end] = t;
				end--;
			}
			else
				beg++;
		}
	}
	else
	{
		while(end >= beg && end < count)
		{
			pPat = m_pData->row(pIndexes[beg]);
			if((int)pPat[attr] == (int)pivot)
			{
				t = pIndexes[beg];
				pIndexes[beg] = pIndexes[end];
				pIndexes[end] = t;
				end--;
			}
			else
				beg++;
		}
	}
	return beg;
}

// static
bool GKdTree::isGreaterOrEqual(const double* pPat, int attr, double pivot)
{
	if(m_pMetric->relation()->valueCount(attr) == 0)
		return (pPat[attr] >= pivot);
	else
		return ((int)pPat[attr] == (int)pivot);
}

GKdNode* GKdTree::buildTree(size_t count, size_t* pIndexes)
{
	int dims = m_pMetric->relation()->size();
	if(count <= (size_t)m_maxLeafSize)
		return new GKdLeafNode((int)count, pIndexes, dims, m_maxLeafSize);

	// Find a good place to split
	double pivot, goodness, p, g;
	int attr = 0;
	computePivotAndGoodness(count, pIndexes, 0, &pivot, &goodness);
	int i;
	for(i = 1; i < dims; i++)
	{
		computePivotAndGoodness(count, pIndexes, i, &p, &g);
		if(g > goodness)
		{
			pivot = p;
			goodness = g;
			attr = i;
		}
	}

	// Split the data
	size_t lessCount = splitIndexes(count, pIndexes, attr, pivot);
	size_t greaterOrEqualCount = count - lessCount;
	if(lessCount == 0 || greaterOrEqualCount == 0)
		return new GKdLeafNode((int)count, pIndexes, dims, m_maxLeafSize);

	// Make an interior node
	GKdNode* pLess = buildTree(lessCount, pIndexes);
	GKdNode* greaterOrEqual = buildTree(greaterOrEqualCount, pIndexes + lessCount);
	return new GKdInteriorNode(dims, pLess, greaterOrEqual, (int)count, attr, pivot);
}

// virtual
size_t GKdTree::addVector(double* pVector)
{
	size_t index = m_pData->rows();
	m_pData->takeRow(pVector);
	m_pRoot->Insert(this, index, pVector);
	if(m_pRoot->IsLeaf())
	{
		if(((GKdLeafNode*)m_pRoot)->GetSize() > m_maxLeafSize)
		{
			size_t* pIndexes = new size_t[((GKdLeafNode*)m_pRoot)->GetSize()];
			ArrayHolder<size_t> hIndexes(pIndexes);
			int used = m_pRoot->Gather(pIndexes);
			GAssert(used == ((GKdLeafNode*)m_pRoot)->GetSize()); // m_size is wrong. This may corrupt memory.
			GKdNode* pNewNode = buildTree(used, pIndexes);
			delete(m_pRoot);
			m_pRoot = pNewNode;
		}
	}
	else
	{
		if(((GKdInteriorNode*)m_pRoot)->m_timeLeft <= 0)
		{
			GKdNode* pNewNode = ((GKdInteriorNode*)m_pRoot)->Rebuild(this);
			delete(m_pRoot);
			m_pRoot = pNewNode;
		}
	}
	return index;
}

// virtual
double* GKdTree::releaseVector(size_t index)
{
	double* pPat = m_pData->row(index);
	m_pRoot->Remove(this, index, pPat);
	size_t last = m_pData->rows() - 1;
	if(index != last)
	{
		double* pPatLast = m_pData->row(last);
		m_pRoot->Rename(this, last, index, pPatLast);
	}
	return m_pData->releaseRow(index); // (releaseRow moves the last row to the index position)
}

class KdTree_Compare_Nodes_Functor
{
public:
	bool operator() (GKdNode* pA, GKdNode* pB) const
	{
		double a = pA->GetMinDist();
		double b = pB->GetMinDist();
		return (a > b);
	}
};

void GKdTree::findNeighbors(size_t* pOutNeighbors, double* pOutSquaredDistances, const double* pInputVector, size_t nExclude)
{
	GClosestNeighborFindingHelper helper(m_neighborCount, pOutNeighbors, pOutSquaredDistances);
	KdTree_Compare_Nodes_Functor comparator;
	priority_queue< GKdNode*, vector<GKdNode*>, KdTree_Compare_Nodes_Functor > q(comparator);
	q.push(m_pRoot);
	while(q.size() > 0)
	{
		GKdNode* pNode = q.top();
		q.pop();
		if(pNode->GetMinDist() >= helper.GetWorstDist())
			break;
		if(pNode->IsLeaf())
		{
			double squaredDist;
			double* pCand;
			vector<size_t>* pIndexes = ((GKdLeafNode*)pNode)->GetIndexes();
			size_t count = pIndexes->size();
			for(size_t i = 0; i < count; i++)
			{
				size_t index = (*pIndexes)[i];
				if(index == nExclude)
					continue;
				pCand = m_pData->row(index);
				squaredDist = m_pMetric->dissimilarity(pInputVector, pCand);
				helper.TryPoint(index, squaredDist);
			}
		}
		else
		{
			int attr;
			double pivot;
			GKdInteriorNode* pParent = (GKdInteriorNode*)pNode;
			pParent->GetDivision(&attr, &pivot);
			GKdNode* pLess = pParent->GetLess();
			pLess->CopyOffset(pParent);
			GKdNode* pGreaterOrEqual = pParent->GetGreaterOrEqual();
			pGreaterOrEqual->CopyOffset(pParent);
			if(isGreaterOrEqual(pInputVector, attr, pivot))
				pLess->AdjustOffset(attr, pInputVector[attr] - pivot, m_pMetric->scaleFactors());
			else
				pGreaterOrEqual->AdjustOffset(attr, pivot - pInputVector[attr], m_pMetric->scaleFactors());
			q.push(pLess);
			q.push(pGreaterOrEqual);
		}
	}
}

// virtual
void GKdTree::neighbors(size_t* pOutNeighbors, size_t index)
{
	GTEMPBUF(double, distances, m_neighborCount);
	neighbors(pOutNeighbors, distances, index);
}

// virtual
void GKdTree::neighbors(size_t* pOutNeighbors, double* pOutDistances, size_t index)
{
	findNeighbors(pOutNeighbors, pOutDistances, m_pData->row(index), index);
}

// virtual
void GKdTree::neighbors(size_t* pOutNeighbors, double* pOutDistances, const double* pInputVector)
{
	findNeighbors(pOutNeighbors, pOutDistances, pInputVector, INVALID_INDEX);
}

// virtual
void GKdTree::reoptimize()
{
	if(!m_pRoot->IsLeaf())
	{
		GKdNode* pNewNode = ((GKdInteriorNode*)m_pRoot)->Rebuild(this);
		delete(m_pRoot);
		m_pRoot = pNewNode;
	}
}

#ifndef NO_TEST_CODE
#	include "GImage.h"
#	include "GHeap.h"

void MeasureBounds(GData* pData, GKdNode* pNode, int attr, double* pMin, double* pMax)
{
	if(pNode->IsLeaf())
	{
		double min = 1e200;
		double max = -1e200;
		vector<size_t>* pIndexes = ((GKdLeafNode*)pNode)->GetIndexes();
		double* pPat;
		for(size_t i = 0; i < pIndexes->size(); i++)
		{
			pPat = pData->row((*pIndexes)[i]);
			min = MIN(pPat[attr], min);
			max = MAX(pPat[attr], max);
		}
		*pMin = min;
		*pMax = max;
	}
	else
	{
		double min1, min2, max1, max2;
		GKdNode* pChild = ((GKdInteriorNode*)pNode)->GetLess();
		MeasureBounds(pData, pChild, attr, &min1, &max1);
		pChild = ((GKdInteriorNode*)pNode)->GetGreaterOrEqual();
		MeasureBounds(pData, pChild, attr, &min2, &max2);
		*pMin = MIN(min1, min2);
		*pMax = MAX(max1, max2);
	}
}

void DrawKdNode(GPlotWindow* pw, GKdNode* pNode, GData* pData)
{
	if(pNode->IsLeaf())
	{
		vector<size_t>* pIndexes = ((GKdLeafNode*)pNode)->GetIndexes();
		double* pPat;
		for(size_t i = 0; i < pIndexes->size(); i++)
		{
			pPat = pData->row((*pIndexes)[i]);
			pw->dot(pPat[0], pPat[1], 5, 0xff00ff00, 0xff000000);
			std::ostringstream os;
			os << (int)(*pIndexes)[i];
			string tmp = os.str();
			pw->label(pPat[0], pPat[1], tmp.c_str(), 1.0f, 0xffffffff);
		}
	}
	else
	{
		int attr;
		double pivot, min, max;
		((GKdInteriorNode*)pNode)->GetDivision(&attr, &pivot);
		if(attr == 0)
		{
			MeasureBounds(pData, pNode, 1, &min, &max);
			pw->line(pivot, min, pivot, max, 0xffff0000);
		}
		else
		{
			GAssert(attr == 1); // unsupported value
			MeasureBounds(pData, pNode, 0, &min, &max);
			pw->line(min, pivot, max, pivot, 0xffff0000);
		}
		GKdNode* pChild = ((GKdInteriorNode*)pNode)->GetLess();
		DrawKdNode(pw, pChild, pData);
		pChild = ((GKdInteriorNode*)pNode)->GetGreaterOrEqual();
		DrawKdNode(pw, pChild, pData);
	}
}

class GDontGoFarMetric : public GDissimilarityMetric
{
public:
	double m_squaredMaxDist;

	GDontGoFarMetric(double maxDist)
	: GDissimilarityMetric(), m_squaredMaxDist(maxDist * maxDist)
	{
	}

	virtual ~GDontGoFarMetric()
	{
	}

	virtual GTwtNode* toTwt(GTwtDoc* pDoc)
	{
		ThrowError("not implemented");
		return NULL;
	}

	virtual void init(sp_relation& pRelation)
	{
		m_pRelation = pRelation;
	}

	virtual double dissimilarity(const double* pA, const double* pB)
	{
		double squaredDist = GVec::squaredDistance(pA, pB, m_pRelation->size());
		if(squaredDist > m_squaredMaxDist)
			ThrowError("a kd-tree shouldn't have to look this far away");
		return squaredDist;
	}
};

void GKdTree_testThatItDoesntLookFar()
{
	GRand prng(0);
	GData tmp(2);
	for(size_t i = 0; i < 100000; i++)
	{
		double* pRow = tmp.newRow();
		pRow[0] = prng.uniform();
		pRow[1] = prng.uniform();
	}
	GDontGoFarMetric metric(0.05);
	GKdTree kdTree(&tmp, 0, 5, &metric, false);
	double row[2];
	size_t neighs[5];
	double dists[5];
	for(size_t i = 0; i < 100; i++)
	{
		row[0] = prng.uniform();
		row[1] = prng.uniform();
		kdTree.neighbors(neighs, dists, row);
	}
}

#	define TEST_DIMS 4
#	define TEST_PATTERNS 1000
#	define TEST_NEIGHBORS 24
// static
void GKdTree::test()
{
	GClosestNeighborFindingHelper::test();
	GKdTree_testThatItDoesntLookFar();

	sp_relation rel;
	rel = new GUniformRelation(TEST_DIMS, 0);
	GHeap heap(2048);
	GData data(rel, &heap);
	GRand prng(0);
	int i, j;
	for(i = 0; i < TEST_PATTERNS; i++)
	{
		double* pPat = data.newRow();
		prng.spherical(pPat, TEST_DIMS);
	}
	GBruteForceNeighborFinder bf(&data, 0, TEST_NEIGHBORS, NULL, true);
	GKdTree kd(&data, 0, TEST_NEIGHBORS, NULL, true);
/*
	GAssert(TEST_DIMS == 2); // You must change TEST_DIMS to 2 if you're going to plot the tree
	GImage image;
	image.SetSize(1000, 1000);
	image.Clear(0xff000000);
	GPlotWindow pw(&image, -1.1, -1.1, 1.1, 1.1);
	DrawKdNode(&pw, kd.GetRoot(), &data);
	image.SavePNGFile("kdtree.png");
*/
	size_t bfNeighbors[TEST_NEIGHBORS];
	size_t kdNeighbors[TEST_NEIGHBORS];
	double bfDistances[TEST_NEIGHBORS];
	double kdDistances[TEST_NEIGHBORS];
	for(i = 0; i < TEST_PATTERNS; i++)
	{
		bf.neighbors(bfNeighbors, bfDistances, i);
		bf.sortNeighbors(bfNeighbors, bfDistances);
		kd.neighbors(kdNeighbors, kdDistances, i);
		kd.sortNeighbors(kdNeighbors, kdDistances);
		for(j = 0; j < TEST_DIMS; j++)
		{
			if(bfNeighbors[j] != kdNeighbors[j])
				ThrowError("wrong answer!");
		}
	}
}
#endif // !NO_TEST_CODE

// --------------------------------------------------------------------------------------------------------





















class GShortcutPrunerAtomicCycleDetector : public GAtomicCycleFinder
{
protected:
	GShortcutPruner* m_pThis;
	int m_thresh;

public:
	GShortcutPrunerAtomicCycleDetector(size_t nodes, GShortcutPruner* pThis, int thresh) : GAtomicCycleFinder(nodes), m_pThis(pThis), m_thresh(thresh)
	{
	}

	virtual ~GShortcutPrunerAtomicCycleDetector()
	{
	}

	virtual bool onDetectAtomicCycle(vector<size_t>& cycle)
	{
		if(cycle.size() >= (size_t)m_thresh)
		{
			m_pThis->onDetectBigAtomicCycle(cycle);
			return false;
		}
		else
			return true;
	}
};

GShortcutPruner::GShortcutPruner(size_t* pNeighborhoods, size_t n, int k)
: m_pNeighborhoods(pNeighborhoods), m_n(n), m_k(k), m_cycleThresh(10), m_subGraphRange(6), m_cuts(0)
{
}

GShortcutPruner::~GShortcutPruner()
{
}

bool GShortcutPruner::isEveryNodeReachable()
{
	GBitTable visited(m_n);
	deque<size_t> q;
	visited.set(0);
	q.push_back(0);
	while(q.size() > 0)
	{
		size_t cur = q.front();
		q.pop_front();
		for(int j = 0; j < m_k; j++)
		{
			size_t neigh = m_pNeighborhoods[m_k * cur + j];
			if(neigh < m_n && !visited.bit(neigh))
			{
				visited.set(neigh);
				q.push_back(neigh);
			}
		}
	}
	for(size_t i = 0; i < m_n; i++)
	{
		if(!visited.bit(i))
			return false;
	}
	return true;
}

size_t GShortcutPruner::prune()
{
	while(true)
	{
		bool everyNodeReachable = isEveryNodeReachable();
		GShortcutPrunerAtomicCycleDetector g(m_n, this, m_cycleThresh);
		size_t* pHood = m_pNeighborhoods;
		for(size_t i = 0; i < m_n; i++)
		{
			for(int j = 0; j < m_k; j++)
			{
				if(pHood[j] < m_n)
					g.addEdgeIfNotDupe(i, pHood[j]);
			}
			pHood += m_k;
		}
		size_t oldCuts = m_cuts;
		g.compute();
		if(everyNodeReachable && !isEveryNodeReachable())
			ThrowError("Cutting shortcuts should not segment the graph");
		if(m_cuts == oldCuts)
			break;
	}
	return m_cuts;
}

void GShortcutPruner::onDetectBigAtomicCycle(vector<size_t>& cycle)
{
	// Make a subgraph containing only nodes close to the cycle
	size_t* mapIn = new size_t[m_n];
	ArrayHolder<size_t> hMapIn(mapIn);
	vector<size_t> mapOut;
	GBitTable visited(m_n);
	deque<size_t> q;
	for(vector<size_t>::iterator it = cycle.begin(); it != cycle.end(); it++)
	{
		q.push_back(*it);
		q.push_back(1);
	}
	while(q.size() > 0)
	{
		size_t cur = q.front();
		q.pop_front();
		size_t depth = q.front();
		q.pop_front();
		mapIn[cur] = mapOut.size();
		mapOut.push_back(cur);
		if(depth <= (size_t)m_subGraphRange)
		{
			for(int j = 0; j < m_k; j++)
			{
				size_t neigh = m_pNeighborhoods[cur * m_k + j];
				if(neigh < m_n && !visited.bit(neigh))
				{
					visited.set(neigh);
					q.push_back(neigh);
					q.push_back(depth + 1);
				}
			}
		}
	}

	// Compute betweenness of all edges
	GBrandesBetweennessCentrality g(mapOut.size());
	for(size_t i = 0; i < mapOut.size(); i++)
	{
		size_t* pHood = m_pNeighborhoods + mapOut[i] * m_k;
		for(int j = 0; j < m_k; j++)
		{
			size_t neigh = pHood[j];
			if(neigh < m_n && visited.bit(neigh))
			{
				g.addDirectedEdgeIfNotDupe(i, mapIn[neigh]);
				g.addDirectedEdgeIfNotDupe(mapIn[neigh], i);
			}
		}
	}
	g.compute();

	// Find the edge on the cycle with the largest betweenness
	size_t shortcutFrom = 0;
	size_t shortcutTo = 0;
	double shortcutBetweenness = 0;
	for(size_t i = 0; i < cycle.size(); i++)
	{
		size_t from = cycle[i];
		size_t to = cycle[(i + 1) % cycle.size()];
		size_t forwIndex = g.neighborIndex(mapIn[from], mapIn[to]);
		size_t revIndex = g.neighborIndex(mapIn[to], mapIn[from]);
		double d = g.edgeBetweenness(mapIn[from], forwIndex) + g.edgeBetweenness(mapIn[to], revIndex);
		if(i == 0 || d > shortcutBetweenness)
		{
			shortcutBetweenness = d;
			shortcutFrom = from;
			shortcutTo = to;
		}
	}

	// Cut the shortcut
	bool cutForward = false;
	for(int j = 0; j < m_k; j++)
	{
		if(m_pNeighborhoods[shortcutFrom * m_k + j] == shortcutTo)
		{
			m_pNeighborhoods[shortcutFrom * m_k + j] = INVALID_INDEX;
			cutForward = true;
			m_cuts++;
			break;
		}
	}
	bool cutReverse = false;
	for(int j = 0; j < m_k; j++)
	{
		if(m_pNeighborhoods[shortcutTo * m_k + j] == shortcutFrom)
		{
			m_pNeighborhoods[shortcutTo * m_k + j] = INVALID_INDEX;
			cutReverse = true;
			m_cuts++;
			break;
		}
	}
	if(!cutForward && !cutReverse)
		ThrowError("Failed to find the offending edge");
}

#ifndef NO_TEST_CODE
// static
void GShortcutPruner::test()
{
	// Make a fully-connected grid
	int w = 6;
	int h = 6;
	int n = w * h;
	int k = 4;
	size_t* pNeighbors = new size_t[n * k];
	ArrayHolder<size_t> hNeighbors(pNeighbors);
	size_t i = 0;
	size_t* pHood = pNeighbors;
	for(int y = 0; y < h; y++)
	{
		for(int x = 0; x < w; x++)
		{
			int j = 0;
			pHood[j++] = (x > 0 ? i - 1 : INVALID_INDEX);
			pHood[j++] = (x < w - 1 ? i + 1 : INVALID_INDEX);
			pHood[j++] = (y > 0 ? i - w : INVALID_INDEX);
			pHood[j++] = (y < h - 1 ? i + w : INVALID_INDEX);
			pHood += k;
			i++;
		}
	}

	// Add 3 shortcuts
	pNeighbors[(0 * w + 0) * k + 0] = n - 1; // connect (0,0) to (w-1, h-1)
	pNeighbors[(0 * w + (w - 1)) * k + 1] = n - 1; // connect (w-1,0) to (w-1,h-1)
	pNeighbors[((h - 1) * w + (w - 1)) * k + 0] = w - 1; // connect (w-1,h-1) to (w-1,0)

	// Cut the shortcuts
	GShortcutPruner pruner(pNeighbors, n, k);
	pruner.setCycleThreshold(h);
	pruner.setSubGraphRange(3);
	size_t cuts = pruner.prune();
	if(pNeighbors[(0 * w + 0) * k + 0] != INVALID_INDEX)
		ThrowError("missed a shortcut");
	if(pNeighbors[(0 * w + (w - 1)) * k + 1] != INVALID_INDEX)
		ThrowError("missed a shortcut");
	if(pNeighbors[((h - 1) * w + (w - 1)) * k + 0] != INVALID_INDEX)
		ThrowError("missed a shortcut");
	if(cuts != 3)
		ThrowError("wrong number of cuts");
}
#endif // NO_TEST_CODE











class GCycleCutAtomicCycleDetector : public GAtomicCycleFinder
{
protected:
	GCycleCut* m_pThis;
	int m_thresh;
	bool m_restore;
	bool m_gotOne;

public:
	GCycleCutAtomicCycleDetector(size_t nodes, GCycleCut* pThis, int thresh, bool restore) : GAtomicCycleFinder(nodes), m_pThis(pThis), m_thresh(thresh), m_restore(restore), m_gotOne(false)
	{
	}

	virtual ~GCycleCutAtomicCycleDetector()
	{
	}

	bool gotOne() { return m_gotOne; }

	virtual bool onDetectAtomicCycle(vector<size_t>& cycle)
	{
		if(cycle.size() >= (size_t)m_thresh)
		{
			if(m_restore)
				m_gotOne = true;
			else
				m_pThis->onDetectBigAtomicCycle(cycle);
			return false;
		}
		else
			return true;
	}
};

GCycleCut::GCycleCut(size_t* pNeighborhoods, GData* pPoints, int k)
: m_pNeighborhoods(pNeighborhoods), m_pPoints(pPoints), m_k(k), m_cycleThresh(10), m_cutCount(0)
{
	// Compute the mean neighbor distance
/*	size_t* pNeigh = m_pNeighborhoods;
	int colCount = m_pPoints->cols();
	size_t count = 0;
	double sum = 0;
	for(size_t i = 0; i < m_pPoints->rows(); i++)
	{
		for(int j = 0; j < k; j++)
		{
			if(*pNeigh < m_pPoints->rows())
			{
				sum += sqrt(GVec::squaredDistance(m_pPoints->row(i), m_pPoints->row(*pNeigh), colCount));
				count++;
			}
			pNeigh++;
		}
	}
	m_aveDist = sum / count;
*/
	// Compute the capacities
	size_t* pNeigh = m_pNeighborhoods;
	for(size_t i = 0; i < m_pPoints->rows(); i++)
	{
		for(int j = 0; j < k; j++)
		{
			if(*pNeigh < m_pPoints->rows())
			{
/*
				double cap = 1.0 / (m_aveDist + sqrt(GVec::squaredDistance(m_pPoints->row(i), m_pPoints->row(*pNeigh), colCount)));
				m_capacities[make_pair(i, *pNeigh)] = cap;
				m_capacities[make_pair(*pNeigh, i)] = cap;
*/
				m_capacities[make_pair(i, *pNeigh)] = 1.0;
				m_capacities[make_pair(*pNeigh, i)] = 1.0;
			}
			pNeigh++;
		}
	}
}

GCycleCut::~GCycleCut()
{
}

bool GCycleCut::doAnyBigAtomicCyclesExist()
{
	// Make the graph
	GCycleCutAtomicCycleDetector g(m_pPoints->rows(), this, m_cycleThresh, true);
	size_t* pHood = m_pNeighborhoods;
	for(size_t i = 0; i < m_pPoints->rows(); i++)
	{
		for(int j = 0; j < m_k; j++)
		{
			if(pHood[j] < m_pPoints->rows())
				g.addEdgeIfNotDupe(i, pHood[j]);
		}
		pHood += m_k;
	}

	// Find a large atomic cycle (calls onDetectBigAtomicCycle when found)
	g.compute();
	return g.gotOne();
}

size_t GCycleCut::cut()
{
	m_cuts.clear();

	// Cut the graph
	while(true)
	{
		// Make the graph
		GCycleCutAtomicCycleDetector g(m_pPoints->rows(), this, m_cycleThresh, false);
		size_t* pHood = m_pNeighborhoods;
		for(size_t i = 0; i < m_pPoints->rows(); i++)
		{
			for(int j = 0; j < m_k; j++)
			{
				if(pHood[j] < m_pPoints->rows())
					g.addEdgeIfNotDupe(i, pHood[j]);
			}
			pHood += m_k;
		}

		// Find a large atomic cycle (calls onDetectBigAtomicCycle when found)
		size_t oldCuts = m_cutCount;
		g.compute();
		if(m_cutCount == oldCuts)
			break;
	}

	// Restore superfluous cuts
	for(vector<size_t>::iterator it = m_cuts.begin(); it != m_cuts.end(); )
	{
		size_t point = *it;
		it++;
		GAssert(it != m_cuts.end()); // broken cuts list
		int neigh = (int)*it;
		it++;
		GAssert(it != m_cuts.end()); // broken cuts list
		size_t other = *it;
		it++;

		// Restore the edge if it doesn't create a big atomic cycle
		m_pNeighborhoods[point * m_k + neigh] = other;
		if(!doAnyBigAtomicCyclesExist())
			m_cutCount--;
		else
			m_pNeighborhoods[point * m_k + neigh] = INVALID_INDEX;
	}

	return m_cutCount;
}

void GCycleCut::onDetectBigAtomicCycle(vector<size_t>& cycle)
{
	// Find the bottleneck
	double bottleneck = 1e308;
	for(size_t i = 0; i < cycle.size(); i++)
	{
		size_t from = cycle[i];
		size_t to = cycle[(i + 1) % cycle.size()];
		pair<size_t, size_t> p = make_pair(from, to);
		double d = m_capacities[p];
		if(i == 0 || d < bottleneck)
			bottleneck = d;
	}
	GAssert(bottleneck > 0); // all capacities should be greater than zero

	// Reduce every edge in the cycle by the bottleneck's capacity
	for(size_t i = 0; i < cycle.size(); i++)
	{
		size_t from = cycle[i];
		size_t to = cycle[(i + 1) % cycle.size()];
		pair<size_t, size_t> p1 = make_pair(from, to);
		pair<size_t, size_t> p2 = make_pair(to, from);
		double d = m_capacities[p1];
		if(d - bottleneck > 1e-12)
		{
			// Reduce the capacity
			m_capacities[p1] = d - bottleneck;
			m_capacities[p2] = d - bottleneck;
		}
		else
		{
			// Remove the edge
			m_capacities.erase(p1);
			m_capacities.erase(p2);
			int forw = -1;
			size_t* pHood = m_pNeighborhoods + from * m_k;
			for(int j = 0; j < m_k; j++)
			{
				if(pHood[j] == to)
				{
					forw = j;
					break;
				}
			}
			int rev = -1;
			pHood = m_pNeighborhoods + to * m_k;
			for(int j = 0; j < m_k; j++)
			{
				if(pHood[j] == from)
				{
					rev = j;
					break;
				}
			}
			GAssert(rev >= 0 || forw >= 0); // couldn't find the edge
			if(forw >= 0)
			{
				m_pNeighborhoods[from * m_k + forw] = INVALID_INDEX;
				m_cuts.push_back(from);
				m_cuts.push_back(forw);
				m_cuts.push_back(to);
				m_cutCount++;
			}
			if(rev >= 0)
			{
				m_pNeighborhoods[to * m_k + rev] = INVALID_INDEX;
				m_cuts.push_back(to);
				m_cuts.push_back(rev);
				m_cuts.push_back(from);
				m_cutCount++;
			}
		}
	}
}

#ifndef NO_TEST_CODE
// static
void GCycleCut::test()
{
	// Make a fully-connected grid
	int w = 6;
	int h = 6;
	int n = w * h;
	int k = 4;
	size_t* pNeighbors = new size_t[n * k];
	ArrayHolder<size_t> hNeighbors(pNeighbors);
	size_t i = 0;
	size_t* pHood = pNeighbors;
	for(int y = 0; y < h; y++)
	{
		for(int x = 0; x < w; x++)
		{
			int j = 0;
			pHood[j++] = (x > 0 ? i - 1 : INVALID_INDEX);
			pHood[j++] = (x < w - 1 ? i + 1 : INVALID_INDEX);
			pHood[j++] = (y > 0 ? i - w : INVALID_INDEX);
			pHood[j++] = (y < h - 1 ? i + w : INVALID_INDEX);
			pHood += k;
			i++;
		}
	}

	// Add 3 shortcuts
	pNeighbors[(0 * w + 0) * k + 0] = n - 1; // connect (0,0) to (w-1, h-1)
	pNeighbors[(0 * w + (w - 1)) * k + 1] = n - 1; // connect (w-1,0) to (w-1,h-1)
	pNeighbors[((h - 1) * w + (w - 1)) * k + 0] = w - 1; // connect (w-1,h-1) to (w-1,0)

	// Make some random data
	GData data(5);
	GRand prng(0);
	for(int i = 0; i < n; i++)
	{
		double* pRow = data.newRow();
		prng.spherical(pRow, 5);
	}

	// Cut the shortcuts
	GCycleCut pruner(pNeighbors, &data, k);
	pruner.setCycleThreshold(h);
	size_t cuts = pruner.cut();
	if(pNeighbors[(0 * w + 0) * k + 0] != INVALID_INDEX)
		ThrowError("missed a shortcut");
	if(pNeighbors[(0 * w + (w - 1)) * k + 1] != INVALID_INDEX)
		ThrowError("missed a shortcut");
	if(pNeighbors[((h - 1) * w + (w - 1)) * k + 0] != INVALID_INDEX)
		ThrowError("missed a shortcut");
	if(cuts != 3)
		ThrowError("wrong number of cuts");
}
#endif // NO_TEST_CODE










//#define DEBUG_SPEW

GManifoldNeighborFinder::GManifoldNeighborFinder(GData* pData, int littleK, int bigK, int intrinsicDims, double alpha, double beta, bool prune, GRand* pRand)
: GNeighborFinder(pData, littleK), m_littleK(littleK)
{
	m_a = alpha;
	m_b = beta;
	m_learningRate = 0.2;
	m_windowSize = 30;
	m_minImprovement = 0.02;
	m_prune = prune;

	// Make a table of the bigK nearest neighbors to each point
	int dims = pData->relation()->size();
#ifdef DEBUG_SPEW
	cout << "k=" << littleK << ", m=" << bigK << ", dims=" << dims << ", t=" << intrinsicDims << ", alpha=" << m_a << ",beta=" << m_b << "\n";
	cout << "Building kd-tree...\n";
#endif
	GKdTree neighborFinder(pData, 0, bigK, NULL, true);
#ifdef DEBUG_SPEW
	cout << "Building neighbor table...\n";
#endif
//	GBruteForceNeighborFinder neighborFinder(pData, bigK, NULL, true, NULL);
	size_t* pBigNeighborTable = new size_t[bigK * pData->rows()];
	ArrayHolder<size_t> hBigNeighborTable(pBigNeighborTable);
	double* pDistances = new double[bigK];
	ArrayHolder<double> hDistances(pDistances);
	size_t* pNeighbors = pBigNeighborTable;
	m_meanSquaredDist = 0.0;
	for(size_t i = 0; i < pData->rows(); i++)
	{
		neighborFinder.neighbors(pNeighbors, pDistances, i);
		neighborFinder.sortNeighbors(pNeighbors, pDistances);
		if(pNeighbors[bigK - 1] >= pData->rows())
			ThrowError("neighbor finder failed to find enough neighbors");
		m_meanSquaredDist += pDistances[bigK / 2];
		pNeighbors += bigK;
	}
	m_meanSquaredDist /= pData->rows();

	// Initialize the weights
	double* pWeights = new double[bigK * pData->rows()];
	ArrayHolder<double> hWeights(pWeights);
	GVec::setAll(pWeights, (double)littleK / bigK, bigK * pData->rows());

	// Find better neighbors
	GDataArray tanSpaces(dims);
	tanSpaces.newSets(pData->rows(), intrinsicDims);
	GTEMPBUF(double, pBuf, dims + bigK);
	double* pDissim = pBuf + dims;
	double besterr = 1e200;
	int window = 0;
	int iters = 0;
	while(true)
	{
#ifdef DEBUG_SPEW
		cout << "--------------- Iter = " << iters << "\n";
#endif
		// Compute the tangeant hyperplane at each point
		double* pNeighborWeights = pWeights;
		GData neighborhood(dims);
		for(size_t i = 0; i < pData->rows(); i++)
		{
			// Make the neighborhood
			neighborhood.flush();
			for(int j = 0; j < bigK; j++)
				neighborhood.copyRow(pData->row(pBigNeighborTable[bigK * i + j]));

			// Compute the tangeant hyperplane
#ifdef DEBUG_SPEW
//			cout << "Tan space " << (int)i << "\n";
#endif
			GData* pTanSpace = tanSpaces.sets()[i];
			for(int j = 0; j < intrinsicDims; j++)
			{
				//neighborhood.principalComponent(pBasis, dims, pData->row(i), pRand);
				double* pBasis = pTanSpace->row(j);
				neighborhood.weightedPrincipalComponent(pBasis, dims, pData->row(i), pNeighborWeights, pRand);
				neighborhood.removeComponent(pData->row(i), pBasis, dims);
#ifdef DEBUG_SPEW
//				cout << " <";
//				GVec::print(stdout, 3, pBasis, dims);
//				cout << ">";
#endif
			}
			pNeighborWeights += bigK;
#ifdef DEBUG_SPEW
//			cout << "\n";
#endif
		}

		// Refine the weights
		pNeighborWeights = pWeights;
		size_t* pHood = pBigNeighborTable;
		double err = 0;
		for(size_t i = 0; i < pData->rows(); i++)
		{
#ifdef DEBUG_SPEW
//			cout << "Dissim " << i << ":";
#endif
			// Compute the dissimilarity with each neighbor
			double* pMe = pData->row(i);
			GData* pMyTan = tanSpaces.sets()[i];
			for(int j = 0; j < bigK; j++)
			{
				size_t neighIndex = pHood[j];
#ifdef DEBUG_SPEW
//				if(j == littleK)
//					cout << " ### ";
//				cout << " " << neighIndex << "(";
#endif
				double* pNeigh = pData->row(neighIndex);
				GData* pNeighTan = tanSpaces.sets()[neighIndex];

				// Compute dissimilarity between the two neighbors
				double sqdist = GVec::squaredDistance(pMe, pNeigh, dims);
				double sqsin = 0;
				if(sqdist > 0)
				{
					pNeighTan->project(pBuf, pMe, pNeigh);
					sqsin = GVec::squaredDistance(pBuf, pMe, dims) / sqdist;
					pMyTan->project(pBuf, pNeigh, pMe);
					sqsin = MAX(sqsin, GVec::squaredDistance(pBuf, pNeigh, dims) / sqdist);
				}
				double d = pMyTan->dihedralCorrelation(pNeighTan, pRand);
#ifdef DEBUG_SPEW
//				cout << "mono=" << sqsin << " , di=" << (1.0 - (d * d));
#endif
				sqsin = MAX(sqsin, 1.0 - (d * d));
				double normDist = sqdist / m_meanSquaredDist;
				pDissim[j] = m_a * sqsin + m_b * sqsin * normDist + normDist;
#ifdef DEBUG_SPEW
//				cout << ", del=" << normDist << ", d=" << pDissim[j];
#endif
			}
#ifdef DEBUG_SPEW
//			cout << "\n";
#endif
			GVec::smallestToFront(pDissim, littleK, bigK, pNeighborWeights, pHood);
			err += GVec::dotProduct(pNeighborWeights, pDissim, bigK);

			// Adjust weights, and move biggest littleK to front
			for(int j = 0; j < littleK; j++)
				pNeighborWeights[j] = (1.0 - m_learningRate) * pNeighborWeights[j] + m_learningRate;
			GVec::sumToOne(pNeighborWeights, bigK);
			GVec::multiply(pNeighborWeights, littleK, bigK);
#ifdef DEBUG_SPEW
//			cout << "Weights " << i << ":";
//			for(int j = 0; j < bigK; j++)
//			{
//				if(j == littleK)
//					cout << " ### ";
//				cout << " " << pHood[j] << "(" << pNeighborWeights[j] << ")";
//			}
//			cout << "\n";
#endif

			// Advance
			pNeighborWeights += bigK;
			pHood += bigK;
		}

		// Detect convergence
#ifdef DEBUG_SPEW
		cout << "g=" << err << " (improvement=" << ((1.0 - err / besterr) * 100.0) << ", window=" << window << "\n";
#endif
		if(1.0 - err / besterr >= m_minImprovement)
		{
			window = 0;
			besterr = err;
		}
		else if(window >= m_windowSize)
			break;
		window++;
		iters++;

	}
#ifdef DEBUG_SPEW
	cout << "Iters: " << iters << "\n";
	cout << "Storing neighbors...\n";
#endif

	// Move biggest weights to the front
	double* pNeighborWeights = pWeights;
	size_t* pHood = pBigNeighborTable;
	for(size_t i = 0; i < pData->rows(); i++)
	{
		for(int j = 0; j < bigK; j++)
			pNeighborWeights[j] = 1.0 - pNeighborWeights[j];
		GVec::smallestToFront(pNeighborWeights, littleK, bigK, NULL, pHood);
		for(int j = 0; j < bigK; j++)
			pNeighborWeights[j] = 1.0 - pNeighborWeights[j];
		pNeighborWeights += bigK;
		pHood += bigK;
	}

	// Store the neighborhood
	m_pNeighborhoods = new size_t[littleK * pData->rows()];
	size_t* pBigNeigh = pBigNeighborTable;
	size_t* pLittleNeigh = m_pNeighborhoods;
	for(size_t i = 0; i < pData->rows(); i++)
	{
		memcpy(pLittleNeigh, pBigNeigh, sizeof(size_t) * littleK);
		pLittleNeigh += littleK;
		pBigNeigh += bigK;
	}

	// Prune the shortcuts
	if(m_prune)
	{
		GCycleCut pruner(m_pNeighborhoods, pData, littleK);
		pruner.cut();
	}
#ifdef DEBUG_SPEW
	cout << "Done finding neighbors.\n";
#endif
}

// virtual
GManifoldNeighborFinder::~GManifoldNeighborFinder()
{
	delete[] m_pNeighborhoods;
}

// virtual
void GManifoldNeighborFinder::neighbors(size_t* pOutNeighbors, size_t index)
{
	memcpy(pOutNeighbors, m_pNeighborhoods + index * m_neighborCount, sizeof(size_t) * m_neighborCount);
}

// virtual
void GManifoldNeighborFinder::neighbors(size_t* pOutNeighbors, double* pOutDistances, size_t index)
{
	neighbors(pOutNeighbors, index);
	for(int i = 0; i < m_neighborCount; i++)
		pOutDistances[i] = 1.0;
}

#ifndef NO_TEST_CODE
#	include "GPlot.h"
#	define NEIGHBORS 4
//#	define SAVE_VISUALIZATIONS

void GManifoldNeighborFinder_Test_Helper(GData* pData, GNeighborFinder* pNF, const char* szFilename)
{
#	ifdef SAVE_VISUALIZATIONS
	GImage image;
	image.setSize(1024, 1024);
	image.clear(0xffffffff);
	GPlotWindow pw(&image, -4, -4, 4, 4);
//	GPlotWindow pw(&image, -.3, -.3, .3, .3);

	// Plot the neighbors
	size_t neighbors[NEIGHBORS];
	for(size_t i = 0; i < pData->rows(); i++)
	{
		double* pPat = pData->row(i);
		pNF->neighbors(neighbors, i);
		unsigned int col = gAHSV(0xff, (float)i / pData->rows(), 1.0f, 0.5f);
		for(int j = 0; j < NEIGHBORS; j++)
		{
			if(neighbors[j] >= pData->rows())
				continue;
			double* pPat2 = pData->row(neighbors[j]);
			pw.line(pPat[0], pPat[1], pPat2[0], pPat2[1], col);
		}
	}

	// Plot the points
	for(size_t i = 0; i < pData->rows(); i++)
	{
		double* pPat = pData->row(i);
		unsigned int col = gAHSV(0xff, (float)i / pData->rows(), 1.0f, 0.5f);
		pw.dot(pPat[0], pPat[1], 3, col, 0xffffffff);
//		pw.dot(pPat[0], pPat[1], 12, col, 0xffffffff);

		// Label the point
		pw.label(pPat[0], pPat[1], gformat(i), 1.0f, 0xff000000);

	}
	image.savePng(szFilename);
#	endif
}

#	ifdef SAVE_VISUALIZATIONS
#		include "../GClasses/G3D.h"
#	endif
void GManifoldNeighborFinder_Test_Helper3(GData* pData, GNeighborFinder* pNF, const char* szFilename)
{
#	ifdef SAVE_VISUALIZATIONS
	GImage image;
	image.setSize(1024, 1024);
	image.clear(0xffffffff);
	GCamera camera;
	camera.setImageSize(image.width(), image.height());
	camera.setViewAngle(PI / 3);
//	camera.setViewAngle(PI / 5);
	G3DVector* pCameraPos = camera.lookFromPoint();
	pCameraPos->set(0, 1, 3);
//	pCameraPos->set(0, 0.5, 1.5);
	G3DVector dir;
	dir.set(0, -1, -3);
	camera.setDirection(&dir, 0.0);

	// Plot the neighbors
	size_t neighbors[NEIGHBORS];
	G3DVector in1, out1, in2, out2;
	for(size_t i = 0; i < pData->rows(); i++)
	{
		double* pPat = pData->row(i);
		in1.set(pPat[0], pPat[1], pPat[2]);
		camera.project(&in1, &out1);
		pNF->neighbors(neighbors, i);
		unsigned int col = gAHSV(0xff, (float)i / pData->rows(), 1.0f, 0.5f);
		for(int j = 0; j < NEIGHBORS; j++)
		{
			if(neighbors[j] >= pData->rows())
				continue;
			double* pPat2 = pData->row(neighbors[j]);
			in2.set(pPat2[0], pPat2[1], pPat2[2]);
			camera.project(&in2, &out2);
			image.line((int)floor(out1.m_vals[0] + 0.5), (int)floor(out1.m_vals[1] + 0.5), (int)floor(out2.m_vals[0] + 0.5), (int)floor(out2.m_vals[1] + 0.5), col);
		}
	}

	// Plot the points
	for(size_t i = 0; i < pData->rows(); i++)
	{
		double* pPat = pData->row(i);
		in1.m_vals[0] = pPat[0];
		in1.m_vals[1] = pPat[1];
		in1.m_vals[2] = pPat[2];
		camera.project(&in1, &out1);
		unsigned int col = gAHSV(0xff, (float)i / pData->rows(), 1.0f, 0.5f);
		float radius = 10.0 / (float)out1.m_vals[2];
		image.dot((float)out1.m_vals[0], (float)out1.m_vals[1], radius, col, 0xffffffff);
/*
		// Label the point
		GRect r((int)out1.m_vals[0], (int)out1.m_vals[1], 50, 16);
		image.text(&r, gformat(i), 0xff000000, 1.0f);
*/
	}
	image.savePng(szFilename);
#	endif
}

size_t GManifoldNeighborFinder_test_countShortcuts(GNeighborFinder* pNF, GData* pData, int k)
{
	// Build the table of neighbors
	size_t n = pData->rows();
	size_t* pNeighborhoods = new size_t[n * k];
	ArrayHolder<size_t> hNeighborhoods(pNeighborhoods);
	size_t* pHood = pNeighborhoods;
	for(size_t i = 0; i < n; i++)
	{
		pNF->neighbors(pHood, i);
		pHood += k;
	}

	// Prune the shortcuts
	GCycleCut pruner(pNeighborhoods, pData, k);
	pruner.setCycleThreshold(8);
	size_t count = pruner.cut();
	return count;
}

// static
void GManifoldNeighborFinder::test()
{
/*
	{
		// SLINKEY
		GData data(3);
		for(double t = 0; t < 7 * PI; t += 0.05)
		{
			double* pPat = data.newRow();
			pPat[0] = sin(t);
			pPat[1] = 0.005 * t;
			pPat[2] = cos(t);
		}
		GRand prng(0);
		GKdTree nf1(&data, NEIGHBORS, NULL, true);
		GManifoldNeighborFinder_Test_Helper3(&data, &nf1, "slinkey1.png");
		if(GManifoldNeighborFinder_test_countShortcuts(&nf1, &data, NEIGHBORS) == 0)
			ThrowError("perfect baseline?");
		GManifoldNeighborFinder nf2(&data, NEIGHBORS, NEIGHBORS * 10,
				1, // intrinsicDims
				24, // a
				0, // b
				false, // prune?
				&prng);
		GManifoldNeighborFinder_Test_Helper3(&data, &nf2, "slinkey2.png");
		if(GManifoldNeighborFinder_test_countShortcuts(&nf2, &data, NEIGHBORS) != 0)
			ThrowError("got a shortcut");
	}
*/
/*
	{
		// HOUR-GLASS
		GData data(2);
		for(double x = -0.5 * PI; x < 0.5 * PI; x += 0.1)
		{
			for(double y = 0; y < PI; y += 0.1)
			{
				double* pPat = data.newRow();
				pPat[0] = 2.0 * x * cos(y);
				pPat[1] = 2.0 * (y - 0.5 * PI);
			}
		}
		GRand prng(0);

		GKdTree nf1(&data, NEIGHBORS, NULL, true);
		GManifoldNeighborFinder_Test_Helper(&data, &nf1, "hourglass1.png");
		if(GManifoldNeighborFinder_test_countShortcuts(&nf1, &data, NEIGHBORS) == 0)
			ThrowError("perfect baseline?");
		GManifoldNeighborFinder nf2(&data, NEIGHBORS, NEIGHBORS * 10,
				1, // intrinsicDims
				2, // a
				1, // b
				false, // prune?
				&prng);
		GManifoldNeighborFinder_Test_Helper(&data, &nf2, "hourglass2.png");
		if(GManifoldNeighborFinder_test_countShortcuts(&nf2, &data, NEIGHBORS) != 0)
			ThrowError("got a shortcut");
	}
*/
	{
		// TWIST
		GData data(2);
		for(double t = M_PI / 4; t < 2 * M_PI - M_PI / 4; t += 0.1)
		{
			double* pPat = data.newRow();
			pPat[0] = sin(t * 2);
			pPat[1] = -2 * cos(t);
		}
		GRand prng(0);
		GKdTree nf1(&data, 0, NEIGHBORS, NULL, true);
		GManifoldNeighborFinder_Test_Helper(&data, &nf1, "twist1.png");
		if(GManifoldNeighborFinder_test_countShortcuts(&nf1, &data, NEIGHBORS) == 0)
			ThrowError("perfect baseline?");
		GManifoldNeighborFinder nf2(&data, NEIGHBORS, NEIGHBORS * 5,
				1, // intrinsicDims
				0.5, // a
				4.0, // b
				false, // prune?
				&prng);
		GManifoldNeighborFinder_Test_Helper(&data, &nf2, "twist2.png");
		if(GManifoldNeighborFinder_test_countShortcuts(&nf2, &data, NEIGHBORS) != 0)
			ThrowError("got a shortcut");
	}
/*
	{
		// TWIST3D
		GRand prng(0);
		GData data(3);
		for(double t = PI / 4; t < 2 * PI - PI / 4; t += 0.001)
		{
			double* pPat = data.newRow();
			pPat[0] = sin(t * 2);
			pPat[1] = -2 * cos(t);
			pPat[2] = 0.5 * prng.uniform() - 0.25;
		}
		GKdTree nf1(&data, NEIGHBORS, NULL, true);
		GManifoldNeighborFinder_Test_Helper3(&data, &nf1, "twist3d1.png");
		if(GManifoldNeighborFinder_test_countShortcuts(&nf1, &data, NEIGHBORS) == 0)
			ThrowError("perfect baseline?");
		GManifoldNeighborFinder nf2(&data, NEIGHBORS, NEIGHBORS * 4,
				2, // intrinsicDims
				1.0, // a
				8.0, // b
				false, // prune?
				&prng);
		GManifoldNeighborFinder_Test_Helper3(&data, &nf2, "twist3d2.png");
		if(GManifoldNeighborFinder_test_countShortcuts(&nf2, &data, NEIGHBORS) != 0)
			ThrowError("got a shortcut");
	}
*/
	{
		// FLASK
		GData data(2);
		for(double t = -M_PI; t < M_PI; t += 0.1)
		{
			double* pPat = data.newRow();
			pPat[0] = 0.5 * (t + 2.1 * sin(t * 2));
			pPat[1] = -2 * cos(t);
		}
		GRand prng(0);
		GKdTree nf1(&data, 0, NEIGHBORS, NULL, true);
		GManifoldNeighborFinder_Test_Helper(&data, &nf1, "flask1.png");
		if(GManifoldNeighborFinder_test_countShortcuts(&nf1, &data, NEIGHBORS) == 0)
			ThrowError("perfect baseline?");
		GManifoldNeighborFinder nf2(&data, NEIGHBORS, NEIGHBORS * 5,
				1, // intrinsicDims
				0.25, // a
				0.5, // b
				false, // prune?
				&prng);
		GManifoldNeighborFinder_Test_Helper(&data, &nf2, "flask2.png");
		if(GManifoldNeighborFinder_test_countShortcuts(&nf2, &data, NEIGHBORS) != 0)
			ThrowError("got a shortcut");
	}
/*
	{
		// ELBOW
		GData data(2);
		for(double t = 0; t < 1; t += 0.1)
		{
			double* pPat = data.newRow();
			pPat[0] = 0.0;
			pPat[1] = 3.0 - 3.0 * t;
		}
		for(double t = 0; t < 1; t += 0.1)
		{
			double* pPat = data.newRow();
			pPat[0] = 3.0 * t;
			pPat[1] = 0;
		}
		GRand prng(0);
		GKdTree nf1(&data, NEIGHBORS, NULL, true);
		GManifoldNeighborFinder_Test_Helper(&data, &nf1, "elbow1.png");
		if(GManifoldNeighborFinder_test_countShortcuts(&nf1, &data, NEIGHBORS) == 0)
			ThrowError("perfect baseline?");
		GManifoldNeighborFinder nf2(&data, NEIGHBORS, NEIGHBORS * 5,
				1, // intrinsicDims
				0.25, // a
				0.0, // b
				false, // prune?
				&prng);
		GManifoldNeighborFinder_Test_Helper(&data, &nf2, "elbow2.png");
		if(GManifoldNeighborFinder_test_countShortcuts(&nf2, &data, NEIGHBORS) != 0)
			ThrowError("got a shortcut");
	}
*/
/*
	{
		// CLOVER
		GData data(2);
		for(double t = 0; t < 4 * PI; t += 0.05)
		{
			double* pPat = data.newRow();
			double scale = cos(t * 3 / 2) + 2;
			pPat[0] = cos(t) * scale;
			pPat[1] = sin(t) * scale;
		}
		GRand prng(0);
		GKdTree nf1(&data, NEIGHBORS, NULL, true);
		GManifoldNeighborFinder_Test_Helper(&data, &nf1, "clover1.png");
		if(GManifoldNeighborFinder_test_countShortcuts(&nf1, &data, NEIGHBORS) == 0)
			ThrowError("perfect baseline?");
		GManifoldNeighborFinder nf2(&data,
				NEIGHBORS, // littleK
				NEIGHBORS * 3, // bigK
				1, // intrinsicDims
				2, // a
				2, // b
				false, // prune?
				&prng);
		GManifoldNeighborFinder_Test_Helper(&data, &nf2, "clover2.png");
		if(GManifoldNeighborFinder_test_countShortcuts(&nf2, &data, NEIGHBORS) != 0)
			ThrowError("got a shortcut");
	}
*/
}
#endif // NO_TEST_CODE








GDynamicSystemNeighborFinder::GDynamicSystemNeighborFinder(GData* pObservations, GData* pActions, bool ownActionsData, int neighborCount, GRand* pRand)
: GNeighborFinder(pObservations, neighborCount), m_ownActionsData(ownActionsData), m_pActions(pActions), m_pRand(pRand)
{
	if(pObservations->rows() != pActions->rows())
		ThrowError("Expected the same number of observations as control vectors");
	if(pActions->cols() != 1)
		ThrowError("Sorry, only one action dim is currently supported");
	int actionValues = m_pActions->relation()->valueCount(0);
	if(actionValues < 2)
		ThrowError("Sorry, only nominal actions are currently supported");

	// Train the consequence maps
	int obsDims = pObservations->cols();
	sp_relation spRel = new GUniformRelation(obsDims * 2, 0);
	for(int j = 0; j < actionValues; j++)
	{
		GData consequenceData(spRel);
		for(size_t i = 0; i < pObservations->rows() - 1; i++)
		{
			if((int)pActions->row(i)[0] == j)
			{
				double* pObs = consequenceData.newRow();
				GVec::copy(pObs, pObservations->row(i), obsDims);
				double* pDelta = pObs + obsDims;
				GVec::copy(pDelta, pObservations->row(i + 1), obsDims);
				GVec::subtract(pDelta, pObs, obsDims);
			}
		}
		GAssert(consequenceData.rows() > 20); // not much data
/*
		GNeuralNet* pNN = new GNeuralNet(pRand);
		pNN->addLayer(40);
		GFilter* pMap = new GFilter(pNN, true);
		pMap->setFeatureTransform(new GNormalize(-2.0, 2.0), true);
		pMap->setLabelTransform(new GNormalize(0.0, 1.0), true);
*/
		GKNN* pMap = new GKNN(1, pRand);
		//pMap->setInterpolationMethod(GKNN::Mean);
//		GDecisionTree* pMap = new GDecisionTree(pRand);
		m_consequenceMaps.push_back(pMap);
		pMap->train(&consequenceData, obsDims);
#ifdef _DEBUG
/*
		string s = "h";
		s += gformat(j);
		s += ".arff";
		consequenceData.saveArff(s.c_str());
*/
/*
		// Check accuracy
		double* pBuf = new double[pObservations->cols()];
		ArrayHolder<double> hBuf(pBuf);
		for(size_t i = 0; i < pObservations->rows() - 1; i++)
		{
			if((int)pActions->row(i)[0] == j)
			{
				m_consequenceMaps[j]->predict(pObservations->row(i), pBuf);
				GVec::add(pBuf, pObservations->row(i), pObservations->cols());
				GVec::subtract(pBuf, pObservations->row(i + 1), pObservations->cols());
				double meanErr = sqrt(GVec::squaredMagnitude(pBuf, pObservations->cols()) / pObservations->cols());
				GAssert(meanErr < 2.0); // bad estimate
			}
		}
*/
#endif
	}
}

// virtual
GDynamicSystemNeighborFinder::~GDynamicSystemNeighborFinder()
{
	if(m_ownActionsData)
		delete(m_pActions);
	for(vector<GSupervisedLearner*>::iterator it = m_consequenceMaps.begin(); it != m_consequenceMaps.end(); it++)
		delete(*it);
}

bool GDynamicSystemNeighborFinder::findPath(size_t from, size_t to, double* path, double maxDist)
{
	// Find the path
	int actionValues = m_pActions->relation()->valueCount(0);
	double* pStart = m_pData->row(from);
	double* pGoal = m_pData->row(to);
	int dims = m_pData->cols();
	double origSquaredDist = GVec::squaredDistance(pStart, pGoal, dims);
	GTEMPBUF(double, pObs, dims + dims + dims);
	double* pDelta = pObs + dims;
	double* pRemaining = pDelta + dims;
	GVec::copy(pObs, pStart, dims);
	GBitTable usedActions(actionValues);
	GVec::setAll(path, 0.0, actionValues);
	while(true)
	{
		GVec::copy(pRemaining, pGoal, dims);
		GVec::subtract(pRemaining, pObs, dims);
		if(GVec::squaredMagnitude(pRemaining, dims) < 1e-9)
			break; // We have arrived at the destination
		double biggestCorr = 1e-6;
		int bestAction = -1;
		double stepSize = 0.0;
		int lastPredicted = -1;
		for(int i = 0; i < actionValues; i++)
		{
			if(usedActions.bit(i))
				continue;
			m_consequenceMaps[i]->predict(pObs, pDelta);
			lastPredicted = i;
			double d = GVec::correlation(pDelta, pRemaining, dims);
			if(d <= 0)
				usedActions.set(i);
			else if(d > biggestCorr)
			{
				biggestCorr = d;
				bestAction = i;
				stepSize = MIN(1.0, GVec::dotProduct(pDelta, pRemaining, dims) / GVec::squaredMagnitude(pDelta, dims));
			}
		}
		if(bestAction < 0)
			break; // There are no good actions, so we're done
		if(stepSize < 1.0)
			usedActions.set(bestAction); // let's not do microscopic zig-zagging

		// Advance the current observation
		if(bestAction != lastPredicted)
			m_consequenceMaps[bestAction]->predict(pObs, pDelta);
		GVec::addScaled(pObs, stepSize, pDelta, dims);
		path[bestAction] += stepSize;
		if(GVec::squaredMagnitude(path, actionValues) > maxDist * maxDist)
			return false;
	}
	if(GVec::squaredMagnitude(pRemaining, dims) >= 0.2 * 0.2 * origSquaredDist)
		return false; // Too imprecise. Throw this one out.
	return true;
}

// virtual
void GDynamicSystemNeighborFinder::neighbors(size_t* pOutNeighbors, size_t index)
{
	GTEMPBUF(double, dissims, m_neighborCount);
	neighbors(pOutNeighbors, dissims, index);
}

// virtual
void GDynamicSystemNeighborFinder::neighbors(size_t* pOutNeighbors, double* pOutDistances, size_t index)
{
	int valueCount = m_pActions->relation()->valueCount(0);
	if(m_pActions->cols() > 1 || valueCount == 0)
		ThrowError("continuous and multi-dim actions not supported yet");
	int actionValues = m_pActions->relation()->valueCount(0);
	int pos = 0;
	GTEMPBUF(double, path, actionValues);
	for(size_t i = 0; pos < m_neighborCount && i < m_pData->rows(); i++)
	{
		if(index == i)
			continue;
		if(!findPath(index, i, path, 2.0 //distCap
			))
		{
			if(index + 1 == i)
			{
				pOutNeighbors[pos] = i;
				pOutDistances[pos] = 1.0;
			}
			else if(index == i + 1)
			{
				pOutNeighbors[pos] = i;
				pOutDistances[pos] = 1.0;
			}
			else
				continue;
		}
		pOutNeighbors[pos] = i;
		pOutDistances[pos] = GVec::squaredMagnitude(path, actionValues);
		//GAssert(ABS(pOutDistances[pos]) < 0.001 || ABS(pOutDistances[pos] - 1.0) < 0.001 || ABS(pOutDistances[pos] - 1.4142) < 0.001 || ABS(pOutDistances[pos] - 2.0) < 0.001); // Noisy result. Does the transition function have noise? If so, then this is expected, so comment me out.
		pos++;
	}

	// Fill the remaining slots with nothing
	while(pos < m_neighborCount)
	{
		pOutNeighbors[pos] = INVALID_INDEX;
		pOutDistances[pos] = 0.0;
		pos++;
	}
}
















/*
GTemporalNeighborFinder::GTemporalNeighborFinder(GData* pObservations, int neighborCount, GRand* pRand)
: GNeighborFinder(pObservations, neighborCount), m_pRand(pRand)
{
	// Make the data
	GMixedRelation* pR = new GMixedRelation();
	sp_relation pRel = pR;
	pR->addAttrs(pObservations->relation().get());
	pR->addAttrs(pObservations->cols(), 0);
	pR->addAttr(0);
	int dims = pObservations->cols();
	GData consequenceData(pRel);
	consequenceData.newRows(2 * pObservations->rows() - 2);
	for(size_t i = 0; i < pObservations->rows(); i++)
	{
		if(i > 0)
		{
			double* pRow = consequenceData.row(2 * i - 1);
			GVec::copy(pRow, pObservations->row(i), dims);
			GVec::copy(pRow + dims, pObservations->row(i - 1), dims);
			GVec::subtract(pRow + dims, pObservations->row(i), dims);
			double mag = sqrt(GVec::squaredMagnitude(pRow + dims, dims));
			GVec::safeNormalize(pRow + dims, dims, m_pRand);
			pRow[dims + dims] = mag;
		}
		if(i < pObservations->rows() - 1)
		{
			double* pRow = consequenceData.row(2 * i);
			GVec::copy(pRow, pObservations->row(i), dims);
			GVec::copy(pRow + dims, pObservations->row(i + 1), dims);
			GVec::subtract(pRow + dims, pObservations->row(i), dims);
			double mag = sqrt(GVec::squaredMagnitude(pRow + dims, dims));
			GVec::safeNormalize(pRow + dims, dims, m_pRand);
			pRow[dims + dims] = mag;
		}
	}

	// Make the model
	GNeuralNet* pNN = new GNeuralNet(m_pRand);
	pNN->addLayer(60); // todo: don't hard-code this
	GFilter* pFilter1 = new GFilter(pNN, true);
	pFilter1->setFeatureTransform(new GNormalize(), true);
	pFilter1->setLabelTransform(new GNormalize(), true);
	GFilter* pFilter2 = new GFilter(pFilter1, true);
	pFilter2->setFeatureTransform(new GNominalToCat(12), true);
	m_pMap = pFilter2;

	// Train the model
	pFilter2->train(&consequenceData, 1);
	m_pBuf = new double[dims + dims];
}

// virtual
GTemporalNeighborFinder::~GTemporalNeighborFinder()
{
	delete(m_pMap);
	delete[] m_pBuf;
}

// virtual
void GTemporalNeighborFinder::neighbors(size_t* pOutNeighbors, size_t index)
{
	neighbors(pOutNeighbors, NULL, index);
}

// virtual
void GTemporalNeighborFinder::neighbors(size_t* pOutNeighbors, double* pOutDistances, size_t index)
{
	int dims = m_pData->cols();
	int pos = 0;
	GVec::copy(m_pBuf, m_pData->row(index), dims);
	double mag;
	for(size_t i = 0; pos < m_neighborCount && i < m_pData->rows(); i++)
	{
		if(index == i)
			continue;
		GVec::copy(m_pBuf + dims, m_pData->row(i), dims);
		GVec::subtract(m_pBuf + dims, m_pBuf, dims);
		double d = sqrt(GVec::squaredMagnitude(m_pBuf + dims, dims));
		GVec::safeNormalize(m_pBuf + dims, dims, m_pRand);
		m_pMap->predict(m_pBuf, &mag);
		d /= mag;
		if(d < 1.5)
		{
			pOutNeighbors[pos] = i;
			if(pOutDistances)
				pOutDistances[pos] = d * d;
			pos++;
		}
	}

	// Fill the remaining slots with nothing
	while(pos < m_neighborCount)
	{
		pOutNeighbors[pos] = INVALID_INDEX;
		if(pOutDistances)
			pOutDistances[pos] = 0.0;
		pos++;
	}
}
*/








GSequenceNeighborFinder::GSequenceNeighborFinder(GData* pData, int neighborCount)
: GNeighborFinder(pData, neighborCount)
{
}

// virtual
GSequenceNeighborFinder::~GSequenceNeighborFinder()
{
}

// virtual
void GSequenceNeighborFinder::neighbors(size_t* pOutNeighbors, size_t index)
{
	return neighbors(pOutNeighbors, NULL, index);
}

// virtual
void GSequenceNeighborFinder::neighbors(size_t* pOutNeighbors, double* pOutDistances, size_t index)
{
	int prevPos = -1;
	int pos = 0;
	int i = 1;
	while(true)
	{
		if(pos == prevPos)
		{
			while(pos < m_neighborCount)
				pOutNeighbors[pos++] = INVALID_INDEX;
			break;
		}
		prevPos = pos;
		if(index - i < m_pData->rows())
		{
			pOutNeighbors[pos] = index - i;
			if(++pos >= m_neighborCount)
				break;
		}
		if(index + i < m_pData->rows())
		{
			pOutNeighbors[pos] = index + i;
			if(++pos >= m_neighborCount)
				break;
		}
		i++;
	}
	if(pOutDistances)
	{
		for(int i = 0; i < m_neighborCount; i++)
			pOutDistances[i] = (double)((i + 2) / 2);
	}
}



} // namespace GClasses


