/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#include "GCluster.h"
#include "GNeighborFinder.h"
#include "GBitTable.h"
#include "GHeap.h"
#include "GMath.h"
#include "GVec.h"
#include <math.h>
#include <stdlib.h>
#include "GNeuralNet.h"
#include "GHillClimber.h"
#include "GBitTable.h"
#include "GKNN.h"
#include "GTime.h"
#include "GGraph.h"
#include "GTwt.h"
#include <iostream>

using namespace GClasses;
using std::cout;
using std::vector;


GAgglomerativeClusterer::GAgglomerativeClusterer(int nTargetClusterCount)
: GClusterer(nTargetClusterCount)
{
	m_pNextCluster = NULL;
	m_pPrevCluster = NULL;
	m_pNextNeighbor = NULL;
	m_pClusters = NULL;
	m_pData = NULL;
	m_nVectorCount = 0;
	m_nCurrentClusterCount = 0;
	m_nTargetClusterCount = nTargetClusterCount;
	m_nFirstCluster = -1;
	m_dBalancingFactor = 5;
}

GAgglomerativeClusterer::~GAgglomerativeClusterer()
{
	delete[] m_pNextCluster;
}

double GAgglomerativeClusterer::computeClusterDistance(int a, int b, bool* pAIsBigger)
{
	GAssert(a != b);
	double* pPatA = m_pData->row(a);
	double* pPatB = m_pData->row(b);
	int best;
	double dClosestDist = GVec::estimateSquaredDistanceWithUnknowns(pPatA, pPatB, m_nDims);
	double d;
	int i, j, aSize, bSize;
	for(i = 0; i < 3; i++)
	{
		// Find the point in cluster A that is closest to the current point in cluster B
		best = a;
		aSize = 1;
		for(j = m_pNextNeighbor[a]; j != a; j = m_pNextNeighbor[j])
		{
			d = GVec::estimateSquaredDistanceWithUnknowns(m_pData->row(j), m_pData->row(b), m_nDims);
			if(d < dClosestDist)
			{
				best = j;
				dClosestDist = d;
			}
			aSize++;
		}
		a = best;

		// Find the point in cluster B that is closest to the current point in cluster A
		best = b;
		bSize = 1;
		for(j = m_pNextNeighbor[b]; j != b; j = m_pNextNeighbor[j])
		{
			d = GVec::estimateSquaredDistanceWithUnknowns(m_pData->row(a), m_pData->row(j), m_nDims);
			if(d < dClosestDist)
			{
				best = j;
				dClosestDist = d;
			}
			bSize++;
		}
		b = best;
	}
	if(aSize >= bSize)
		*pAIsBigger = true;
	else
		*pAIsBigger = false;
	return dClosestDist * pow((double)(aSize + bSize), m_dBalancingFactor);
}

void GAgglomerativeClusterer::findMerges()
{
	// Find the nearest neighbor of every cluster
	GHeap heap(1000);
	GData data(3, &heap);
	data.reserve(m_nCurrentClusterCount);
	int i, j, posi, posj;
	double* pVec;
	for(i = m_nFirstCluster; i >= 0; i = m_pNextCluster[i])
	{
		pVec = data.newRow();
		pVec[0] = i;
		pVec[1] = -1;
		pVec[2] = 1e200;
	}
	posi = 0;
	double d;
	bool firstIsBigger;
	for(i = m_nFirstCluster; i >= 0; i = m_pNextCluster[i])
	{
		posj = posi + 1;
		for(j = m_pNextCluster[i]; j >= 0; j = m_pNextCluster[j])
		{
			d = computeClusterDistance(i, j, &firstIsBigger);
			if(firstIsBigger)
			{
				pVec = data.row(posi);
				GAssert(pVec[0] == i);
				if(d < pVec[2])
				{
					pVec[1] = j;
					pVec[2] = d;
				}
			}
			else
			{
				pVec = data.row(posj);
				GAssert(pVec[0] == j);
				if(d < pVec[2])
				{
					pVec[1] = i;
					pVec[2] = d;
				}
			}
			posj++;
		}
		posi++;
	}

	// Add the best merges to the queue
	data.sort(2);
	int nMax = MAX(1, m_nCurrentClusterCount / 4);
	for(i = 0; i < m_nCurrentClusterCount && (int)m_q.size() < nMax; i++)
	{
		pVec = data.row(i);
		if(pVec[2] >= 1e200)
		{
			GAssert(i > 0);
			break;
		}
		m_q.push_back((int)pVec[0]);
		m_q.push_back((int)pVec[1]);
	}
}

void GAgglomerativeClusterer::merge()
{
	// Determine which two clusters should be merged
	if(m_nCurrentClusterCount < 2)
		ThrowError("There's only one cluster");
	if(m_q.size() == 0)
	{
		findMerges();
		GAssert(m_q.size() > 0);
	}
	m_pClusters[0] = -1;
	int a1 = m_q.front();
	m_q.pop_front();
	int b1 = m_q.front();
	m_q.pop_front();
	GAssert(a1 != b1);
	GAssert(a1 >= 0 && b1 >= 0);

	// Merge them
	if(m_pNextCluster[a1] == -2 || m_pNextCluster[b1] == -2)
		return; // one of them isn't a cluster head anymore
	int a2 = m_pNextNeighbor[a1];
	int b2 = m_pNextNeighbor[b1];
	m_pNextNeighbor[a1] = b2;
	m_pNextNeighbor[b1] = a2;

	// Drop one of the clusters from the list
	m_nCurrentClusterCount--;
	if(m_pPrevCluster[b1] >= 0)
	{
		GAssert(m_pNextCluster[b1] >= -1);
		m_pNextCluster[m_pPrevCluster[b1]] = m_pNextCluster[b1];
	}
	else
	{
		GAssert(m_pNextCluster[b1] >= -1);
		m_nFirstCluster = m_pNextCluster[b1];
	}
	if(m_pNextCluster[b1] >= 0)
	{
		GAssert(m_pPrevCluster[b1] >= -1);
		m_pPrevCluster[m_pNextCluster[b1]] = m_pPrevCluster[b1];
	}
	m_pPrevCluster[b1] = -2;
	m_pNextCluster[b1] = -2;
	GAssert(m_nFirstCluster >= 0);
}

// virtual
void GAgglomerativeClusterer::cluster(GData* pData)
{
	if(!pData->relation()->areContinuous(0, pData->cols()))
		ThrowError("GAgglomerativeClusterer doesn't support nominal attributes. You should filter with the NominalToCat transform to convert nominals to real values.");
	m_nDims = pData->relation()->size();
	m_q.clear();
	m_pData = pData;
	m_nVectorCount = (int)m_pData->rows();
	m_nCurrentClusterCount = (int)m_pData->rows();
	m_nFirstCluster = 0;
	delete[] m_pNextCluster;
	m_pNextCluster = new int[4 * m_nVectorCount];
	m_pPrevCluster = m_pNextCluster + m_nVectorCount;
	m_pNextNeighbor = m_pPrevCluster + m_nVectorCount;
	m_pClusters = m_pNextNeighbor + m_nVectorCount;
	m_pClusters[0] = -1;
	int i;
	for(i = 0; i < m_nVectorCount; i++)
	{
		m_pNextCluster[i] = i + 1;
		m_pPrevCluster[i] = i - 1;
		m_pNextNeighbor[i] = i;
	}
	m_pNextCluster[m_nVectorCount - 1] = -1;
	while(m_nCurrentClusterCount > m_nTargetClusterCount)
		merge();
}

// virtual
int GAgglomerativeClusterer::whichCluster(size_t nVector)
{
	GAssert(nVector >= 0 && nVector < (size_t)m_nVectorCount);
	if(m_pClusters[0] < 0)
	{
		int j, head;
		int nCluster = 0;
		for(head = m_nFirstCluster; head >= 0; head = m_pNextCluster[head])
		{
			j = head;
			while(true)
			{
				m_pClusters[j] = nCluster;
				j = m_pNextNeighbor[j];
				if(j == head)
					break;
			}
			nCluster++;
		}
		GAssert(m_pClusters[0] >= 0);
	}
	GAssert(m_pClusters[nVector] >= 0 && m_pClusters[nVector] < m_nCurrentClusterCount);
	return m_pClusters[nVector];
}

#ifndef NO_TEST_CODE

#define SPRIAL_POINTS 250
#define SPIRAL_HEIGHT 3

#include "GImage.h"
//#include "G3D.h"

// static
void GAgglomerativeClusterer::test()
{
	// Make a 3D data set with 3 entwined spirals
	GHeap heap(1000);
	GData data(3, &heap);
	double dThirdCircle = M_PI * 2 / 3;
	double* pVector;
	double rads;
	int i;
	for(i = 0; i < SPRIAL_POINTS; i += 3)
	{
		rads = (double)i * 2 * M_PI / SPRIAL_POINTS;

		pVector = data.newRow();
		pVector[0] = cos(rads);
		pVector[2] = sin(rads);
		pVector[1] = (double)i * SPIRAL_HEIGHT / SPRIAL_POINTS;

		pVector = data.newRow();
		pVector[0] = cos(rads + dThirdCircle);
		pVector[2] = sin(rads + dThirdCircle);
		pVector[1] = (double)i * SPIRAL_HEIGHT / SPRIAL_POINTS;

		pVector = data.newRow();
		pVector[0] = cos(rads + dThirdCircle + dThirdCircle);
		pVector[2] = sin(rads + dThirdCircle + dThirdCircle);
		pVector[1] = (double)i * SPIRAL_HEIGHT / SPRIAL_POINTS;
	}

	// Cluster the points
	GAgglomerativeClusterer clust(3);
	clust.cluster(&data);

	// Test for correctness
	if(clust.whichCluster(0) == clust.whichCluster(1))
		ThrowError("failed");
	if(clust.whichCluster(1) == clust.whichCluster(2))
		ThrowError("failed");
	if(clust.whichCluster(2) == clust.whichCluster(0))
		ThrowError("failed");
	for(i = 3; i < SPRIAL_POINTS; i += 3)
	{
		if(clust.whichCluster(i) != clust.whichCluster(0))
			ThrowError("Wrong cluster");
		if(clust.whichCluster(i + 1) != clust.whichCluster(1))
			ThrowError("Wrong cluster");
		if(clust.whichCluster(i + 2) != clust.whichCluster(2))
			ThrowError("Wrong cluster");
	}

/*  // Uncomment this to make a spiffy visualization of the entwined spirals

	// Draw the classifications
	GImage image;
	image.SetSize(600, 600);
	image.Clear(0xff000000);
	GCamera camera;
	camera.SetViewAngle(PI / 2);
	camera.GetLookFromPoint()->Set(2, 1.5, 3);
	camera.GetLookDirection()->Set(-2, 0, -3);
	camera.ComputeSideVector();
	camera.SetImageSize(600, 600);
	double* pVec;
	G3DVector point;
	GColor col = 0;
	for(i = 0; i < SPRIAL_POINTS; i++)
	{
		pVec = data.row(i);
		point.Set(pVec[0], pVec[1], pVec[2]);
		switch(clust.whichCluster(i))
		{
			case 0: col = 0xffff0000; break;
			case 1: col = 0xff00ff00; break;
			case 2: col = 0xff0000ff; break;
		}
		image.Draw3DLine(&point, &point, &camera, col);
	}

	// Draw the bounding box
	G3DVector point2;
	int x, y, z;
	for(z = 0; z < 2; z++)
	{
		for(y = 0; y < 2; y++)
		{
			for(x = 0; x < 2; x++)
			{
				if(x == 0)
				{
					point.Set(-1, 3 * y, 2 * z - 1);
					point2.Set(1, 3 * y, 2 * z - 1);
					image.Draw3DLine(&point, &point2, &camera, 0xff808080);
				}
				if(y == 0)
				{
					point.Set(2 * x - 1, 0, 2 * z - 1);
					point2.Set(2 * x - 1, 3, 2 * z - 1);
					image.Draw3DLine(&point, &point2, &camera, 0xff808080);
				}
				if(z == 0)
				{
					point.Set(2 * x - 1, 3 * y, -1);
					point2.Set(2 * x - 1, 3 * y, 1);
					image.Draw3DLine(&point, &point2, &camera, 0xff808080);
				}
			}
		}
	}
	image.SavePNGFile("spirals.png");
*/
}
#endif // !NO_TEST_CODE

// -----------------------------------------------------------------------------------------

GAgglomerativeTransducer::GAgglomerativeTransducer()
: GTransducer()
{
	m_pNextCluster = NULL;
	m_pPrevCluster = NULL;
	m_pNextNeighbor = NULL;
	m_pData = NULL;
	m_nVectorCount = 0;
	m_nCurrentClusterCount = 0;
	m_nFirstCluster = -1;
	m_dBalancingFactor = 1;
}

GAgglomerativeTransducer::~GAgglomerativeTransducer()
{
	delete[] m_pNextCluster;
}

double GAgglomerativeTransducer::computeClusterDistance(int a, int b, bool* pAIsBigger)
{
	GAssert(a != b); // same point
	double* pPatA = m_pData->row(a);
	double* pPatB = m_pData->row(b);
	if(pPatA[m_featureDims] != UNKNOWN_DISCRETE_VALUE && pPatB[m_featureDims] != UNKNOWN_DISCRETE_VALUE)
	{
		if(pPatA[m_featureDims] == pPatB[m_featureDims])
			return 0;
		else
			return 1e200;
	}
	int best;
	double dClosestDist = GVec::estimateSquaredDistanceWithUnknowns(pPatA, pPatB, m_featureDims);
	double d;
	int i, j, aSize, bSize;
	for(i = 0; i < 3; i++)
	{
		// Find the point in cluster A that is closest to the current point in cluster B
		best = a;
		aSize = 1;
		for(j = m_pNextNeighbor[a]; j != a; j = m_pNextNeighbor[j])
		{
			d = GVec::estimateSquaredDistanceWithUnknowns(m_pData->row(j), m_pData->row(b), m_featureDims);
			if(d < dClosestDist)
			{
				best = j;
				dClosestDist = d;
			}
			aSize++;
		}
		a = best;

		// Find the point in cluster B that is closest to the current point in cluster A
		best = b;
		bSize = 1;
		for(j = m_pNextNeighbor[b]; j != b; j = m_pNextNeighbor[j])
		{
			d = GVec::estimateSquaredDistanceWithUnknowns(m_pData->row(a), m_pData->row(j), m_featureDims);
			if(d < dClosestDist)
			{
				best = j;
				dClosestDist = d;
			}
			bSize++;
		}
		b = best;
	}
	if(aSize >= bSize)
		*pAIsBigger = true;
	else
		*pAIsBigger = false;
	return dClosestDist * pow((double)(aSize + bSize), m_dBalancingFactor);
}

void GAgglomerativeTransducer::findMerges()
{
	// Find the nearest neighbor of every cluster
	GHeap heap(1000);
	GData data(3, &heap);
	data.reserve(m_nCurrentClusterCount);
	int i, j, posi, posj;
	double* pVec;
	for(i = m_nFirstCluster; i >= 0; i = m_pNextCluster[i])
	{
		pVec = data.newRow();
		pVec[0] = i;
		pVec[1] = -1;
		pVec[2] = 1e200;
	}
	posi = 0;
	double d;
	bool firstIsBigger;
	for(i = m_nFirstCluster; i >= 0; i = m_pNextCluster[i])
	{
/*
if(posi % 20 == 0)
{
cout << ((float)posi / m_nCurrentClusterCount) << "     \r";
cout.flush();
}
*/
		posj = posi + 1;
		for(j = m_pNextCluster[i]; j >= 0; j = m_pNextCluster[j])
		{
			d = computeClusterDistance(i, j, &firstIsBigger);
			if(firstIsBigger)
			{
				pVec = data.row(posi);
				GAssert(pVec[0] == i);
				if(d < pVec[2])
				{
					pVec[1] = j;
					pVec[2] = d;
				}
			}
			else
			{
				pVec = data.row(posj);
				GAssert(pVec[0] == j);
				if(d < pVec[2])
				{
					pVec[1] = i;
					pVec[2] = d;
				}
			}
			posj++;
		}
		posi++;
	}

	// Add the best merges to the queue
	data.sort(2);
	int nMax = MAX(1, m_nCurrentClusterCount / 4);
	for(i = 0; i < m_nCurrentClusterCount && (int)m_q.size() < nMax; i++)
	{
		pVec = data.row(i);
		if(pVec[2] >= 1e200)
			break;
		m_q.push_back((int)pVec[0]);
		m_q.push_back((int)pVec[1]);
	}
}

bool GAgglomerativeTransducer::merge()
{
	// Determine which two clusters should be merged
	if(m_nCurrentClusterCount < 2)
		ThrowError("There's only one cluster");
	if(m_q.size() <= 0)
		findMerges();
	if(m_q.size() == 0)
		return false;
	int a1 = m_q.front();
	m_q.pop_front();
	int b1 = m_q.front();
	m_q.pop_front();
	return mergeRows(a1, b1);
}

bool GAgglomerativeTransducer::mergeRows(int a1, int b1)
{
	GAssert(a1 != b1);
	GAssert(a1 >= 0 && b1 >= 0);

	// Merge the classification
	double* pPatA = m_pData->row(a1);
	double* pPatB = m_pData->row(b1);
	if(pPatA[m_featureDims] != pPatB[m_featureDims])
	{
		int nHeadOfUnknownCluster;
		double classification;
		if(pPatA[m_featureDims] == UNKNOWN_DISCRETE_VALUE)
		{
			nHeadOfUnknownCluster = a1;
			classification = pPatB[m_featureDims];
		}
		else if(pPatB[m_featureDims] == UNKNOWN_DISCRETE_VALUE)
		{
			nHeadOfUnknownCluster = b1;
			classification = pPatA[m_featureDims];
		}
		else
			return true; // don't merge opposing classes
		GAssert(classification >= 0);
		int j = nHeadOfUnknownCluster;
		while(true)
		{
			pPatA = m_pData->row(j);
			pPatA[m_featureDims] = classification;
			j = m_pNextNeighbor[j];
			if(j == nHeadOfUnknownCluster)
				break;
		}
	}

	// Merge them
	if(m_pNextCluster[a1] == -2 || m_pNextCluster[b1] == -2)
		return true; // one of them isn't a cluster head anymore
	int a2 = m_pNextNeighbor[a1];
	int b2 = m_pNextNeighbor[b1];
	m_pNextNeighbor[a1] = b2;
	m_pNextNeighbor[b1] = a2;

	// Drop one of the clusters from the list
	m_nCurrentClusterCount--;
	if(m_pPrevCluster[b1] >= 0)
	{
		GAssert(m_pNextCluster[b1] >= -1);
		m_pNextCluster[m_pPrevCluster[b1]] = m_pNextCluster[b1];
	}
	else
	{
		GAssert(m_pNextCluster[b1] >= -1);
		m_nFirstCluster = m_pNextCluster[b1];
	}
	if(m_pNextCluster[b1] >= 0)
	{
		GAssert(m_pPrevCluster[b1] >= -1);
		m_pPrevCluster[m_pNextCluster[b1]] = m_pPrevCluster[b1];
	}
	m_pPrevCluster[b1] = -2;
	m_pNextCluster[b1] = -2;
	GAssert(m_nFirstCluster >= 0);

	return true;
}

void GAgglomerativeTransducer::mergeLikeClasses(GData* pDataLabeled)
{
	int n;
	double* pPat;
	int classes = pDataLabeled->relation()->valueCount(m_featureDims);
	for(int i = 0; i < classes; i++)
	{
		n = -1;
		for(size_t j = 0; j < m_pData->rows(); j++)
		{
			pPat = m_pData->row(j);
			if((int)pPat[m_featureDims] != i)
				continue;
			if(n < 0)
				n = (int)j;
			else
				mergeRows(n, (int)j);
		}
	}
}

// virtual
void GAgglomerativeTransducer::transduce(GData* pDataLabeled, GData* pDataUnlabeled, int labelDims)
{
	if(labelDims != 1)
		ThrowError("Sorry, only 1 nominal label is supported");
	m_featureDims = pDataLabeled->cols() - labelDims;
	sp_relation pRelation = pDataLabeled->relation();
	if(!pRelation->areContinuous(0, m_featureDims))
		ThrowError("GAgglomerativeTransducer only supports real features.");
	if(!pRelation->areNominal(m_featureDims, 1))
		ThrowError("GAgglomerativeTransducer only supports one nominal label.");
	if(pDataLabeled->cols() != pDataUnlabeled->cols())
		ThrowError("Expected both datasets to have the same number of columns");
	m_featureDims = pRelation->size() - labelDims;

	// Merge into one dataset with unlabeled rows marked with unknown values
	GData both(pDataLabeled->relation());
	both.reserve(pDataLabeled->rows() + pDataUnlabeled->rows());
	for(size_t i = 0; i < pDataLabeled->rows(); i++)
		both.takeRow(pDataLabeled->row(i));
	for(size_t i = 0; i < pDataUnlabeled->rows(); i++)
	{
		double* pPat = pDataUnlabeled->row(i);
		both.takeRow(pPat);
		for(int j = 0; j < labelDims; j++)
		{
			if(pRelation->valueCount(m_featureDims + j) == 0)
				pPat[m_featureDims + j] = UNKNOWN_REAL_VALUE;
			else
				pPat[m_featureDims + j] = UNKNOWN_DISCRETE_VALUE;
		}
	}
	GReleaseDataHolder hBoth(&both);

	m_q.clear();
	m_pData = &both;
	m_nVectorCount = (int)m_pData->rows();
	m_nCurrentClusterCount = (int)m_pData->rows();
	m_nFirstCluster = 0;
	delete[] m_pNextCluster;
	m_pNextCluster = new int[3 * m_nVectorCount];
	m_pPrevCluster = m_pNextCluster + m_nVectorCount;
	m_pNextNeighbor = m_pPrevCluster + m_nVectorCount;
	int i;
	for(i = 0; i < m_nVectorCount; i++)
	{
		m_pNextCluster[i] = i + 1;
		m_pPrevCluster[i] = i - 1;
		m_pNextNeighbor[i] = i;
	}
	m_pNextCluster[m_nVectorCount - 1] = -1;
	int nTargetClusterCount = pRelation->valueCount(m_featureDims);
	mergeLikeClasses(pDataLabeled);
	while(m_nCurrentClusterCount > nTargetClusterCount)
	{
		if(!merge())
			break;
	}
}

#ifndef NO_TEST_CODE
void GAgglomerativeTransducer::makeStateImage(const char* szFilename)
{
	GImage image;
	image.setSize(800, 800);
	image.clear(0xff000000);

	// Measure the ranges
	double xmin, xrange, ymin, yrange;
	m_pData->minAndRange(0, &xmin, &xrange);
	m_pData->minAndRange(1, &ymin, &yrange);

	// Draw the connecting lines
	int x1, y1, x2, y2;
	double* pPat;
	int head;
	int prev = 0;
	int nCluster = 0;
	for(head = m_nFirstCluster; head >= 0; head = m_pNextCluster[head])
	{
		int j = head;
		while(true)
		{
			if(j != head)
			{
				pPat = m_pData->row(prev);
				x1 = (int)floor(GData::normalize(pPat[0], xmin, xrange, 0, 800));
				y1 = (int)floor(GData::normalize(pPat[1], ymin, yrange, 0, 800));
				pPat = m_pData->row(j);
				x2 = (int)floor(GData::normalize(pPat[0], xmin, xrange, 0, 800));
				y2 = (int)floor(GData::normalize(pPat[1], ymin, yrange, 0, 800));
				image.line(x1, y1, x2, y2, 0xff808080);
			}
			prev = j;
			j = m_pNextNeighbor[j];
			if(j == head)
				break;
		}
		nCluster++;
	}

	// Draw the points
	int output;
	int outputCount = m_pData->relation()->valueCount(m_featureDims);
	for(size_t j = 0; j < m_pData->rows(); j++)
	{
		pPat = m_pData->row(j);
		output = (int)pPat[m_featureDims] + 1;
		x1 = (int)floor(GData::normalize(pPat[0], xmin, xrange, 0, 800));
		y1 = (int)floor(GData::normalize(pPat[1], ymin, yrange, 0, 800));
		image.circleFill(x1, y1, 2, gAHSV(0xff, (float)output / (outputCount + 1), 1.0f, 1.0f));
	}

	image.savePng(szFilename);
}
#endif // !NO_TEST_CODE

// -----------------------------------------------------------------------------------------

GKMeans::GKMeans(size_t nClusters, GRand* pRand)
: GClusterer(nClusters)
{
	m_pRand = pRand;
	m_nClusters = nClusters;
	m_pClusters = NULL;
}

GKMeans::~GKMeans()
{
	delete[] m_pClusters;
}

bool GKMeans::selectSeeds(GData* pSeeds)
{
	// Randomly select the seed points
	size_t i, j, k, index;
	double* pVector;
	for(i = 0; i < m_nClusters; i++)
	{
		for(j = 0; j < m_nClusters; j++)
		{
			// Pick a point
			index = (size_t)m_pRand->next(m_pData->rows());
			pVector = m_pData->row(index);

			// Make sure we didn't pick the same point already
			bool bOK = true;
			for(k = 0; k < i; k++)
			{
				if(GVec::squaredDistance(pVector, pSeeds->row(k), m_nDims) <= 0)
				{
					bOK = false;
					break;
				}
			}

			// Accept this seed point
			if(bOK)
			{
				pSeeds->copyRow(pVector);
				break;
			}
		}
		if(j >= m_nClusters)
			return false; // failed to find "m_nClusters" unique seed points
	}
	return true;
}

bool GKMeans::clusterAttempt(int nMaxIterations)
{
	// Pick the seeds
	GHeap heap(1000);
	GData means(m_nDims, &heap);
	means.reserve(m_nClusters);
	if(!selectSeeds(&means))
		return false;
	GAssert(means.rows() == m_nClusters);

	// Do the clustering
	GTEMPBUF(int, pCounts, means.rows());
	int i, nCluster;
	double d, dBestDist;
	double* pVector;
	bool bChanged;
	for(i = 0; i < nMaxIterations; i++)
	{
		// Assign new cluster to each point
		bChanged = false;
		for(size_t j = 0; j < m_pData->rows(); j++)
		{
			pVector = m_pData->row(j);
			dBestDist = 1e200;
			nCluster = -1;
			for(size_t k = 0; k < m_nClusters; k++)
			{
				d = GVec::squaredDistance(pVector, means.row(k), m_nDims);
				if(d < dBestDist)
				{
					dBestDist = d;
					nCluster = k;
				}
			}
			if(m_pClusters[j] != nCluster)
				bChanged = true;
			m_pClusters[j] = nCluster;
		}
		if(!bChanged)
			break;

		// Update the seeds
		for(size_t j = 0; j < means.rows(); j++)
			GVec::setAll(means.row(j), 0.0, m_nDims);
		memset(pCounts, '\0', sizeof(int) * means.rows());
		for(size_t j = 0; j < m_pData->rows(); j++)
		{
			GVec::add(means.row(m_pClusters[j]), m_pData->row(j), m_nDims);
			pCounts[m_pClusters[j]]++;
		}
		for(size_t j = 0; j < means.rows(); j++)
			GVec::multiply(means.row(j), 1.0 / pCounts[j], m_nDims);
	}
	return (i < nMaxIterations);
}

// virtual
void GKMeans::cluster(GData* pData)
{
	if(!pData->relation()->areContinuous(0, pData->cols()))
		ThrowError("GKMeans doesn't support nominal attributes. You should filter with the NominalToCat transform to convert nominal values to reals.");
	m_nDims = pData->relation()->size();
	if(pData->rows() < m_nClusters)
		throw "The number of clusters must be less than the number of data points";
	delete[] m_pClusters;
	m_pData = pData;
	m_pClusters = new int[pData->rows()];
	memset(m_pClusters, 0xff, sizeof(int) * pData->rows());
	int i;
	for(i = 0; i < 5; i++)
	{
		if(clusterAttempt(m_nDims * (int)pData->rows()))
			break;
	}
	GAssert(i < 5);
}

// virtual
int GKMeans::whichCluster(size_t nVector)
{
	GAssert(nVector < m_pData->rows());
	return m_pClusters[nVector];
}


// -----------------------------------------------------------------------------------------

GKMedoids::GKMedoids(int clusters)
: GClusterer(clusters)
{
	m_pMedoids = new size_t[clusters];
	m_pMetric = NULL;
	m_pData = NULL;
}

// virtual
GKMedoids::~GKMedoids()
{
	delete[] m_pMedoids;
	delete(m_pMetric);
}

void GKMedoids::setMetric(GDissimilarityMetric* pMetric)
{
	delete(m_pMetric);
	m_pMetric = pMetric;
}

double GKMedoids::curErr()
{
	double err = 0;
	for(size_t i = 0; i < m_pData->rows(); i++)
	{
		whichCluster(i);
		err += m_d;
	}
	return err;
}

// virtual
void GKMedoids::cluster(GData* pData)
{
	m_pData = pData;
	if(!m_pMetric)
		setMetric(new GRowDistance());
	m_pMetric->init(pData->relation());
	if(pData->rows() < (size_t)m_clusterCount)
		ThrowError("Fewer data point than clusters");
	for(int i = 0; i < m_clusterCount; i++)
		m_pMedoids[i] = i;
	double err = curErr();
	while(true)
	{
		bool improved = false;
		for(size_t i = 0; i < pData->rows(); i++)
		{
			// See if it's already a medoid
			int j;
			for(j = 0; j < m_clusterCount; j++)
			{
				if(m_pMedoids[j] == i)
					break;
			}
			if(j < m_clusterCount)
				continue;

			// Try this point in place of each medoid
			for(j = 0; j < m_clusterCount; j++)
			{
				size_t old = m_pMedoids[j];
				m_pMedoids[j] = i;
				double cand = curErr();
				if(cand < err)
				{
					err = cand;
					break;
				}
				else
					m_pMedoids[j] = old;
			}
		}
		if(!improved)
			break;
	}
}

// virtual
int GKMedoids::whichCluster(size_t nVector)
{
	double* pVec = m_pData->row(nVector);
	int cluster = 0;
	m_d = m_pMetric->dissimilarity(pVec, m_pData->row(m_pMedoids[0]));
	for(int i = 1; i < m_clusterCount; i++)
	{
		double d = m_pMetric->dissimilarity(pVec, m_pData->row(m_pMedoids[i]));
		if(d < m_d)
		{
			m_d = d;
			cluster = i;
		}
	}
	return cluster;
}


// -----------------------------------------------------------------------------------------
/*
void BlurVector(int nDims, double* pInput, double* pOutput, double dAmount)
{
	double dWeight, dSumWeight;
	int i, j;
	for(i = 0; i < nDims; i++)
	{
		pOutput[i] = 0;
		dSumWeight = 0;
		for(j = 0; j < nDims; j++)
		{
			dWeight = GMath::gaussian((double)(j - i) / dAmount);
			dSumWeight += dWeight;
			pOutput[i] += dWeight * pInput[j];
		}
		pOutput[i] /= dSumWeight;
	}
}
*/
/*
void MakeHistogramWithGaussianParzenWindow(int nElements, double* pInput, double* pOutput, double dBlurAmount)
{
	int i, j;
	for(i = 0; i < nElements; i++)
	{
		pOutput[i] = 0;
		for(j = 0; j < nElements; j++)
			pOutput[i] += GMath::gaussian((pOutput[j] - pOutput[i]) / dBlurAmount);
	}
}

int CountLocalMaximums(int nElements, double* pData)
{
	if(nElements < 2)
		return nElements;
	int nCount = 0;
	if(pData[0] > pData[1])
		nCount++;
	int i;
	nElements--;
	for(i = 1; i < nElements; i++)
	{
		if(pData[i] > pData[i - 1] && pData[i] > pData[i + 1])
			nCount++;
	}
	if(pData[nElements] > pData[nElements - 1])
		nCount++;
	return nCount;
}
*/

// -----------------------------------------------------------------------------------------

GGraphCutTransducer::GGraphCutTransducer(int neighborCount, GRand* pRand)
: GTransducer(), m_neighborCount(neighborCount), m_pRand(pRand)
{
	m_pNeighbors = new size_t[neighborCount];
	m_pDistances = new double[neighborCount];
}

// virtual
GGraphCutTransducer::~GGraphCutTransducer()
{
	delete[] m_pNeighbors;
	delete[] m_pDistances;
}

// virtual
void GGraphCutTransducer::transduce(GData* pDataLabeled, GData* pDataUnlabeled, int labelDims)
{
	if(labelDims != 1)
		ThrowError("Only 1 nominal label dim is supported");

	// Use k-NN to compute a distance metric with good scale factors for prediction
	int featureDims = pDataLabeled->cols() - labelDims;
	GKNN knn(m_neighborCount, m_pRand);
	knn.setOptimizeScaleFactors(true);
	knn.train(pDataLabeled, labelDims);
	GRowDistanceScaled* pMetric = knn.metric();

	// Merge into one dataset and build a kd-tree
	GData both(pDataLabeled->relation());
	both.reserve(pDataLabeled->rows() + pDataUnlabeled->rows());
	for(size_t i = 0; i < pDataLabeled->rows(); i++)
		both.takeRow(pDataLabeled->row(i));
	for(size_t i = 0; i < pDataUnlabeled->rows(); i++)
	{
		double* pPat = pDataUnlabeled->row(i);
		pPat[featureDims] = 0;
		both.takeRow(pPat);
	}
	GReleaseDataHolder hBoth(&both);
	GMixedRelation* pRelInput = new GMixedRelation();
	sp_relation relInput;
	relInput = pRelInput;
	pRelInput->addAttrs(pDataLabeled->relation().get(), 0, featureDims);
	GRowDistanceScaled metric2;
	GKdTree neighborFinder(&both, labelDims, m_neighborCount, &metric2, false);
	GVec::copy(metric2.scaleFactors(), pMetric->scaleFactors(), featureDims);

	// Use Graph-cut to separate out each label value
	int valueCount = pDataLabeled->relation()->valueCount(featureDims);
	for(int val = 1; val < valueCount; val++)
	{
		// Add neighborhood edges
		GGraphCut gc(2 + (int)pDataLabeled->rows() + (int)pDataUnlabeled->rows());
		for(size_t i = 0; i < both.rows(); i++)
		{
			neighborFinder.neighbors(m_pNeighbors, m_pDistances, i);
			for(int j = 0; j < m_neighborCount; j++)
			{
				if(m_pNeighbors[j] >= both.rows())
					continue;
				gc.addEdge(2 + (int)i, 2 + (int)m_pNeighbors[j], (float)(1.0 / MAX(sqrt(m_pDistances[j]), 1e-9))); // connect neighbors
			}
		}

		// Add source and sink edges
		for(size_t i = 0; i < pDataLabeled->rows(); i++)
		{
			double* pPat = pDataLabeled->row(i);
			if((int)pPat[featureDims] == val)
				gc.addEdge(0, 2 + (int)i, 1e12f); // connect to source
			else
				gc.addEdge(1, 2 + (int)i, 1e12f); // connect to sink
		}

		// Cut
		gc.cut(0, 1);

		// Label the unlabeled rows
		for(size_t i = 0; i < pDataUnlabeled->rows(); i++)
		{
			double* pPat = pDataUnlabeled->row(i);
			if(gc.isSource(2 + (int)pDataLabeled->rows() + (int)i))
				pPat[featureDims] = (double)val;
		}
	}
}






/*
GNeuralTransducer::GNeuralTransducer(GRand* pRand)
: m_pRand(pRand)
{
	m_pNN = new GNeuralNet(m_pRand);
	m_pNN->setLearningRate(0.005);
}

// virtual
GNeuralTransducer::~GNeuralTransducer()
{
	delete(m_pNN);
}

void GNeuralTransducer::setParams(std::vector<size_t>& ranges)
{
	m_paramRanges.resize(ranges.size());
	vector<size_t>::iterator itSrc = ranges.begin();
	vector<size_t>::iterator itDst = m_paramRanges.begin();
	while(itSrc != ranges.end())
		*(itDst++) = *(itSrc++);
}

// virtual
void GNeuralTransducer::transduce(GData* pDataLabeled, GData* pDataUnlabeled, int labelDims)
{
	if(labelDims != 1)
		ThrowError("Sorry, only one label dim is currently supported");
	if(pDataLabeled->cols() != pDataUnlabeled->cols())
		ThrowError("mismatching number of columns");
	size_t featureDims = pDataLabeled->cols() - 1;
	int labelValues = pDataLabeled->relation()->valueCount(featureDims);
	if(labelValues < 1)
		ThrowError("expected the labels to be nominal");

	// Compute the number of pixels and channels
	size_t pixels = 1;
	for(size_t i = 0; i < m_paramRanges.size(); i++)
		pixels *= m_paramRanges[i];
	size_t channels = featureDims / pixels;
	if((channels * pixels) != featureDims)
		ThrowError("params don't align with the number of feature dims");
	size_t paramDims = m_paramRanges.size();

	// Make the initial cluster data
	GData outData(paramDims + labelValues);
	size_t totalRows = pDataLabeled->rows() + pDataUnlabeled->rows();
	outData.newRows(pDataLabeled->rows() + pDataUnlabeled->rows());
	for(size_t i = 0; i < pDataLabeled->rows(); i++)
	{
		double* pVec = outData.row(i);
		int label = (int)pDataLabeled->row(i)[featureDims];
		if(label >= 0 && label < labelValues)
		{
			GVec::setAll(pVec, 0.0, labelValues);
			pVec[label] = 1.0;
		}
		else
			GVec::setAll(pVec, 1.0 / labelValues, labelValues);
	}
	{
		double dist;
		size_t neigh;
		GKdTree neighborFinder(pDataLabeled, 1, 1, NULL, false);
		for(size_t i = 0; i < pDataUnlabeled->rows(); i++)
		{
			double* pVec = outData.row(pDataLabeled->rows() + i);
			neighborFinder.neighbors(&neigh, &dist, pDataUnlabeled->row(i));
			int label = (int)pDataLabeled->row(neigh)[featureDims];
			if(label >= 0 && label < labelValues)
			{
				GVec::setAll(pVec, 0.0, labelValues);
				pVec[label] = 1.0;
			}
			else
				GVec::setAll(pVec, 1.0 / labelValues, labelValues);
		}
	}

	// Prepare for incremental learning
	sp_relation pRel = new GUniformRelation(paramDims + labelValues + channels);
	m_pNN->enableIncrementalLearning(pRel, channels, NULL, NULL);

	// Iterate
	GBackProp* pBP = m_pNN->backProp();
	GBackPropLayer& bpLayer = pBP->layer(0);
	GNeuralNetLayer& nnLayer = m_pNN->getLayer(0);
	GCoordVectorIterator cvi(m_paramRanges);
	double startTime = GTime::seconds();
	while(true)
	{
		size_t index = (size_t)m_pRand->next(totalRows);
		if(index == 0 && GTime::seconds() - startTime > 18 * 60 * 60)
			break;
		double* pRow = index < pDataLabeled->rows() ? pDataLabeled->row(index) : pDataUnlabeled->row(index - pDataLabeled->rows());
		double* pVec = outData.row(index);
		cvi.setRandom(m_pRand);
		cvi.currentNormalized(pVec);

		// Train the weights
		m_pNN->trainIncremental(pVec, pRow + channels * cvi.currentIndex());

		// Train the contexts
		for(int i = 0; i < labelValues; i++)
		{
			for(size_t j = 0; j < nnLayer.m_neurons.size(); j++)
				pVec[paramDims + i] += m_pNN->learningRate() * bpLayer.m_neurons[j].m_error * nnLayer.m_neurons[j].m_weights[1 + i];
		}

		// Use semi-supervision with the labeled contexts
		if(index < pDataLabeled->rows())
		{
			int mi = GVec::indexOfMax(pVec + paramDims, labelValues, m_pRand);
			int ti = (int)pRow[featureDims];
			if(ti != mi)
			{
				for(int i = 0; i < labelValues; i++)
				{
					if(i == ti)
						pVec[paramDims + i] += m_pNN->learningRate() * (1.0 - pVec[paramDims + i]);
					else
						pVec[paramDims + i] += m_pNN->learningRate() * (0.0 - pVec[paramDims + i]);
				}
			}
		}

		// Regularize the weights
		m_pNN->decayWeights(0.001);
	}

GTwtDoc doc;
doc.setRoot(m_pNN->toTwt(&doc));
doc.save("cluster_model.twt");

	// Deterimine the most likely label for each row
	for(size_t i = 0; i < pDataUnlabeled->rows(); i++)
	{
		double* pVec = outData.row(pDataLabeled->rows() + i);
		int mi = GVec::indexOfMax(pVec + paramDims, labelValues, m_pRand);
		pDataUnlabeled->row(i)[featureDims] = (double)mi;
	}
}
*/
