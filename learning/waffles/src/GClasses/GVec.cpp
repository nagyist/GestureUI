/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#include "GVec.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "GRand.h"
#include "GMacros.h"
#include "GData.h"
#include "GBits.h"
#include "GTwt.h"
#include "GMath.h"
#include "GImage.h"
#include "GBitTable.h"

namespace GClasses {

using std::vector;

// static
bool GVec::doesContainUnknowns(const double* pVector, int nSize)
{
	int n;
	for(n = 0; n < nSize; n++)
	{
		if(*pVector == UNKNOWN_REAL_VALUE)
			return true;
		pVector++;
	}
	return false;
}

// static
void GVec::copy(double* pDest, const double* pSource, int nDims)
{
	memcpy(pDest, pSource, sizeof(double) * nDims);
}

// static
double GVec::dotProduct(const double* pA, const double* pB, int nSize)
{
	double d = 0;
	while(nSize > 0)
	{
		d += *(pA++) * *(pB++);
		nSize--;
	}
	return d;
}

// static
double GVec::dotProduct(const double* pOrigin, const double* pTarget, const double* pVector, int nSize)
{
	double d = 0;
	while(nSize > 0)
	{
		d += (*(pTarget++) - *(pOrigin++)) * (*(pVector++));
		nSize--;
	}
	return d;
}

// static
double GVec::dotProduct(const double* pOriginA, const double* pTargetA, const double* pOriginB, const double* pTargetB, int nSize)
{
	double dVal = 0;
	for(int n = 0; n < nSize; n++)
	{
		dVal += (*pTargetA - *pOriginA) * (*pTargetB - *pOriginB);
		pTargetA++;
		pOriginA++;
		pTargetB++;
		pOriginB++;
	}
	return dVal;
}

// static
double GVec::dotProductIgnoringUnknowns(const double* pOrigin, const double* pTarget, const double* pVector, int nSize)
{
	double dVal = 0;
	for(int n = 0; n < nSize; n++)
	{
		GAssert(pOrigin[n] != UNKNOWN_REAL_VALUE && pVector[n] != UNKNOWN_REAL_VALUE); // unknowns in pOrigin or pVector not supported
		if(pTarget[n] != UNKNOWN_REAL_VALUE)
			dVal += (pTarget[n] - pOrigin[n]) * pVector[n];
	}
	return dVal;
}

// static
double GVec::squaredDistance(const double* pA, const double* pB, int nDims)
{
	double dist = 0;
	double d;
	for(int n = 0; n < nDims; n++)
	{
		d = (*pA) - (*pB);
		dist += (d * d);
		pA++;
		pB++;
	}
	return dist;
}

// static
double GVec::estimateSquaredDistanceWithUnknowns(const double* pA, const double* pB, int nDims)
{
	double dist = 0;
	double d;
	int nMissing = 0;
	int n;
	for(n = 0; n < nDims; n++)
	{
		if(pA[n] == UNKNOWN_REAL_VALUE || pB[n] == UNKNOWN_REAL_VALUE)
			nMissing++;
		else
		{
			d = pA[n] - pB[n];
			dist += (d * d);
		}
	}
	if(nMissing >= nDims)
		return 1e50; // we have no info, so let's make a wild guess
	else
		return dist * nDims / (nDims - nMissing);
}

// static
double GVec::squaredMagnitude(const double* pVector, int nSize)
{
	double dMag = 0;
	while(nSize > 0)
	{
		dMag += ((*pVector) * (*pVector));
		pVector++;
		nSize--;
	}
	return dMag;
}

// static
double GVec::minkowskiMagnitude(double norm, const double* pVector, int nSize)
{
	double dMag = 0;
	int i;
	for(i = 0; i < nSize; i++)
		dMag += pow(ABS(pVector[i]), norm);
	return pow(dMag, 1.0 / norm);
}

// static
double GVec::minkowskiDistance(double norm, const double* pA, const double* pB, int dims)
{
	double dist = 0;
	int i;
	for(i = 0; i < dims; i++)
	{
		dist += pow(ABS(*pA - *pB), norm);
		pA++;
		pB++;
	}
	return pow(dist, 1.0 / norm);
}

// static
void GVec::minkowskiNormalize(double norm, double* pVector, int nSize)
{
	double dMag = minkowskiMagnitude(norm, pVector, nSize);
	int i;
	for(i = 0; i < nSize; i++)
		pVector[i] /= dMag;
}


// static
double GVec::correlation(const double* pA, const double* pB, int nDims)
{
	double dDotProd = dotProduct(pA, pB, nDims);
	if(dDotProd == 0)
		return 0;
	return dDotProd / (sqrt(squaredMagnitude(pA, nDims) * squaredMagnitude(pB, nDims)));
}

// static
double GVec::correlation(const double* pOriginA, const double* pTargetA, const double* pB, int nDims)
{
	double dDotProd = dotProduct(pOriginA, pTargetA, pB, nDims);
	if(dDotProd == 0)
		return 0;
	return dDotProd / (sqrt(squaredDistance(pOriginA, pTargetA, nDims) * squaredMagnitude(pB, nDims)));
}

// static
double GVec::correlation(const double* pOriginA, const double* pTargetA, const double* pOriginB, const double* pTargetB, int nDims)
{
	double dDotProd = dotProduct(pOriginA, pTargetA, pOriginB, pTargetB, nDims);
	if(dDotProd == 0)
		return 0;
	return dDotProd / (sqrt(squaredDistance(pOriginA, pTargetA, nDims) * squaredDistance(pOriginB, pTargetB, nDims)));
}

// static
void GVec::normalize(double* pVector, int nSize)
{
	double dMag = squaredMagnitude(pVector, nSize);
	if(dMag <= 0)
		ThrowError("Can't normalize a vector with zero magnitude");
	GVec::multiply(pVector, 1.0  / sqrt(dMag), nSize);
}

// static
void GVec::safeNormalize(double* pVector, int nSize, GRand* pRand)
{
	double dMag = squaredMagnitude(pVector, nSize);
	if(dMag <= 0)
		pRand->spherical(pVector, nSize);
	else
		GVec::multiply(pVector, 1.0  / sqrt(dMag), nSize);
}

// static
void GVec::sumToOne(double* pVector, int size)
{
	double sum = GVec::sumElements(pVector, size);
	if(sum == 0)
		GVec::setAll(pVector, 1.0 / size, size);
	else
		GVec::multiply(pVector, 1.0 / sum, size);
}

// static
int GVec::indexOfMin(const double* pVector, int dims, GRand* pRand)
{
	int index = 0;
	int count = 1;
	int n;
	for(n = 1; n < dims; n++)
	{
		if(pVector[n] <= pVector[index])
		{
			if(pVector[n] == pVector[index])
			{
				count++;
				if(pRand->next(count) == 0)
					index = n;
			}
			else
			{
				index = n;
				count = 1;
			}
		}
	}
	return index;
}

// static
int GVec::indexOfMax(const double* pVector, int dims, GRand* pRand)
{
	int index = 0;
	int count = 1;
	int n;
	for(n = 1; n < dims; n++)
	{
		if(pVector[n] >= pVector[index])
		{
			if(pVector[n] == pVector[index])
			{
				count++;
				if(pRand->next(count) == 0)
					index = n;
			}
			else
			{
				index = n;
				count = 1;
			}
		}
	}
	return index;
}

// static
int GVec::indexOfMaxMagnitude(const double* pVector, int dims, GRand* pRand)
{
	int index = 0;
	int count = 1;
	int n;
	for(n = 1; n < dims; n++)
	{
		if(ABS(pVector[n]) >= ABS(pVector[index]))
		{
			if(ABS(pVector[n]) == ABS(pVector[index]))
			{
				count++;
				if(pRand->next(count) == 0)
					index = n;
			}
			else
			{
				index = n;
				count = 1;
			}
		}
	}
	return index;
}

// static
void GVec::add(double* pDest, const double* pSource, int nDims)
{
	for(int i = 0; i < nDims; i++)
	{
		*pDest += *pSource;
		pDest++;
		pSource++;
	}
}

// static
void GVec::addScaled(double* pDest, double dMag, const double* pSource, int nDims)
{
	for(int i = 0; i < nDims; i++)
	{
		*pDest += (dMag * (*pSource));
		pDest++;
		pSource++;
	}
}

// static
void GVec::addLog(double* pDest, const double* pSource, int nDims)
{
	int i;
	for(i = 0; i < nDims; i++)
		pDest[i] += log(pSource[i]);
}

// static
void GVec::subtract(double* pDest, const double* pSource, int nDims)
{
	for(int i = 0; i < nDims; i++)
	{
		*pDest -= *pSource;
		pDest++;
		pSource++;
	}
}

// static
void GVec::multiply(double* pVector, double dScalar, int nDims)
{
	for(int i = 0; i < nDims; i++)
	{
		*pVector *= dScalar;
		pVector++;
	}
}

// static
void GVec::pairwiseMultiply(double* pDest, double* pOther, int dims)
{
	while(dims > 0)
	{
		*(pDest++) *= *(pOther++);
		dims--;
	}
}

// static
void GVec::pairwiseDivide(double* pDest, double* pOther, int dims)
{
	while(dims > 0)
	{
		*(pDest++) /= *(pOther++);
		dims--;
	}
}

// static
void GVec::setAll(double* pVector, double value, int dims)
{
	for(int i = 0; i < dims; i++)
	{
		*pVector = value;
		pVector++;
	}
}

void GVec::mostEccentricMaxs(int* pOutPoints, int nPoints, bool bMax, double* pVector, int nDims)
{
	if(nPoints < 1)
		return;

	// Handle the case where we don't have enough points
	int i;
	if(nDims <= nPoints)
	{
		for(i = 0; i < nPoints; i++)
			pOutPoints[i] = i * nDims / nPoints;
		return;
	}

	// Find the most eccentric of the eccentric points until we have fewer than we need
	int nSign = bMax ? 1 : -1;
	GTEMPBUF(int, pWorkBuffers, nDims * 2);
	int* pPrev = &pWorkBuffers[nDims];
	int* pCurrent = pWorkBuffers;
	for(i = 0; i < nDims; i++)
		pCurrent[i] = i;
	int nPointCount = nDims;
	int nPrevPointCount = 0;
	int a, b, c, nPos, maxGap;
	while(nPointCount > nPoints)
	{
		// Swap the buffers
		int* pTemp = pPrev;
		pPrev = pCurrent;
		pCurrent = pTemp;

		// Compute the max gap
		maxGap = nDims * 5 / nPointCount; // 5 is the max allowed gap ratio. Smaller values restrict to more even-ness (and increase computation time)

		// Find all the local maxes
		nPos = 0;
		a = pPrev[0];
		b = pPrev[1];
		if(b - a > maxGap || (pVector[a] - pVector[b]) * nSign > 0)
			pCurrent[nPos++] = a;
		for(i = 1; i < nPointCount - 1; i++)
		{
			a = pPrev[i - 1];
			b = pPrev[i];
			c = pPrev[i + 1];
			if(b - a > maxGap || c - b > maxGap || GBits::sign((pVector[b] - pVector[a]) * nSign) + GBits::sign((pVector[b] - pVector[c]) * nSign) > 0)
				pCurrent[nPos++] = b;
		}
		a = pPrev[i - 1];
		b = pPrev[i];
		if(b - a > maxGap || (pVector[b] - pVector[a]) * nSign > 0)
			pCurrent[nPos++] = b;

		nPrevPointCount = nPointCount;
		nPointCount = nPos;
	}

	// Use all the points in pCurrent and some of the points in pPrev
	GAssert(nPointCount <= nPoints && nPrevPointCount > nPoints); // unexpected array sizes
	nPos = 0;
	int nPrevPos = 0;
	int nNeeded;
	for(i = 0; i < nPoints; i++)
	{
		// Skip excess points in pPrev
		nNeeded = nPoints - i;
		while(		(nPos >= nPointCount || pPrev[nPrevPos] != pCurrent[nPos]) &&
				nPrevPointCount - nPrevPos > nNeeded &&
				i * nPrevPointCount / nPoints > nPrevPos
			)
			nPrevPos++;

		// Take the next point
		pOutPoints[i] = pPrev[nPrevPos];
		if(nPos < nPointCount && pCurrent[nPos] == pPrev[nPrevPos])
			nPos++;
		nPrevPos++;
	}
}

void GVec::mostEccentricPoints(int* pOutPoints, int nPoints, double* pVector, int nDims)
{
	// Get a set of maxs and a set of mins
	int nMins = nPoints / 2;
	int nMaxs = nMins;
	if(nMaxs + nMins < nPoints)
		nMaxs++;
	GTEMPBUF(int, pMaxs, nPoints);
	int* pMins = &pMaxs[nMaxs];
	GVec::mostEccentricMaxs(pMaxs, nMaxs, true, pVector, nDims);
	GVec::mostEccentricMaxs(pMins, nMins, false, pVector, nDims);

	// Merge the points
	int i;
	int nMaxPos = 0;
	int nMinPos = 0;
	for(i = 0; i < nPoints; i++)
	{
		if(nMaxPos >= nMaxs)
			pOutPoints[i] = pMins[nMinPos++];
		else if(nMinPos >= nMins)
			pOutPoints[i] = pMaxs[nMaxPos++];
		else if(pMins[nMinPos] < pMaxs[nMaxPos])
			pOutPoints[i] = pMins[nMinPos++];
		else
			pOutPoints[i] = pMaxs[nMaxPos++];
	}
}

void GVec::interpolateIndexes(int nIndexes, double* pInIndexes, double* pOutIndexes, float fRatio, int nCorrIndexes, double* pCorrIndexes1, double* pCorrIndexes2)
{
	GAssert(nCorrIndexes >= 2); // need at least two correlated indexes (at least the two extremes)
	int nCorr = 0;
	double fInvRatio = (float)1 - fRatio;
	double fIndex, fWeight, f0, f1;
	int i;
	for(i = 0; i < nIndexes; i++)
	{
		fIndex = pInIndexes[i];
		while(nCorr < nCorrIndexes - 2 && fIndex >= pCorrIndexes1[nCorr + 1])
			nCorr++;
		fWeight = (fIndex - pCorrIndexes1[nCorr]) / (pCorrIndexes1[nCorr + 1] - pCorrIndexes1[nCorr]);
		f0 = fInvRatio * pCorrIndexes1[nCorr] + fRatio * pCorrIndexes2[nCorr];
		f1 = fInvRatio * pCorrIndexes1[nCorr + 1] + fRatio * pCorrIndexes2[nCorr + 1];
		pOutIndexes[i] = ((float)1 - fWeight) * f0 + fWeight * f1;
	}
}

void GVec::rotate(double* pVector, int nDims, double dAngle, const double* pA, const double* pB)
{
	// Check that the vectors are orthogonal
	GAssert(pVector != pA && pVector != pB); // expected different vectors
	GAssert(ABS(GVec::dotProduct(pA, pB, nDims)) < 1e-4); // expected orthogonal plane axes

	// Remove old planar component
	double x = GVec::dotProduct(pVector, pA, nDims);
	double y = GVec::dotProduct(pVector, pB, nDims);
	GVec::addScaled(pVector, -x, pA, nDims);
	GVec::addScaled(pVector, -y, pB, nDims);

	// Rotate
	double dRadius = sqrt(x * x + y * y);
	double dTheta = atan2(y, x);
	dTheta += dAngle;
	x = dRadius * cos(dTheta);
	y = dRadius * sin(dTheta);

	// Add new planar component
	GVec::addScaled(pVector, x, pA, nDims);
	GVec::addScaled(pVector, y, pB, nDims);
}

int GVec::valueIndex(double* pVector, int nDims, double dVal)
{
	for(int i = 0; i < nDims; i++)
	{
		if(*pVector == dVal)
			return i;
		pVector++;
	}
	return -1;
}

void GVec::addInterpolatedFunction(double* pOut, int nOutVals, double* pIn, int nInVals)
{
	if(nInVals > nOutVals)
	{
		int inPos = 0;
		int outPos, n, count;
		double d;
		for(outPos = 0; outPos < nOutVals; outPos++)
		{
			n = outPos * nInVals / nOutVals;
			d = 0;
			count = 0;
			while(inPos <= n)
			{
				d += pIn[inPos++];
				count++;
			}
			pOut[outPos] += d / count;
		}
	}
	else if(nInVals < nOutVals)
	{
		double d;
		int n, i, j;
		for(n = 0; n < nOutVals; n++)
		{
			d = (double)n * nInVals / nOutVals;
			i = (int)d;
			j = MIN(i + 1, nInVals - 1);
			d -= i;
			pOut[n] += ((1.0 - d) * pIn[i] + d * pIn[j]);
		}
	}
	else
	{
		for(int n = 0; n < nOutVals; n++)
			pOut[n] += pIn[n];
	}
}

// static
GTwtNode* GVec::toTwt(GTwtDoc* pDoc, const double* pVec, int dims)
{
	GTwtNode* pNode = pDoc->newList(dims);
	int i;
	for(i = 0; i < dims; i++)
		pNode->setItem(i, pDoc->newDouble(pVec[i]));
	return pNode;
}

// static
void GVec::fromTwt(double* pVec, int dims, GTwtNode* pNode)
{
	if(dims != (int)pNode->itemCount())
		ThrowError("Expected ", gformat(dims), " dims, but the twt node specified ", gformat(pNode->itemCount()), " dims");
	int i;
	for(i = 0; i < dims; i++)
		pVec[i] = pNode->item(i)->asDouble();
}

// static
void GVec::print(std::ostream& stream, int precision, double* pVec, int dims)
{
	if(dims == 0)
		return;
	stream.precision(precision);
	stream << *pVec;
	pVec++;
	for(int i = 1; i < dims; i++)
	{
		stream << ", ";
		stream << *pVec;
		pVec++;
	}
}

void GVec::project(double* pDest, const double* pPoint, const double* pOrigin, const double* pBasis, int basisCount, int dims)
{
	GVec::copy(pDest, pOrigin, dims);
	for(int j = 0; j < basisCount; j++)
	{
		GVec::addScaled(pDest, GVec::dotProduct(pOrigin, pPoint, pBasis, dims), pBasis, dims);
		pBasis += dims;
	}
}

void GVec::subtractComponent(double* pInOut, const double* pBasis, int dims)
{
	double component = dotProduct(pInOut, pBasis, dims);
	for(int i = 0; i < dims; i++)
	{
		*pInOut -= *pBasis * component;
		pBasis++;
		pInOut++;
	}
}

void GVec::subtractComponent(double* pInOut, const double* pOrigin, const double* pTarget, int dims)
{
	double component = dotProduct(pInOut, pOrigin, pTarget, dims) / squaredDistance(pOrigin, pTarget, dims);
	for(int i = 0; i < dims; i++)
	{
		*pInOut -= (*pTarget - *pOrigin) * component;
		pTarget++;
		pOrigin++;
		pInOut++;
	}
}

double GVec::sumElements(const double* pVec, int dims)
{
	double sum = 0;
	while(dims > 0)
	{
		sum += *pVec;
		pVec++;
		dims--;
	}
	return sum;
}


void GVec_InsertionSort(double* pVec, int size, double* pParallel1, size_t* pParallel2)
{
	for(int i = 1; i < size; i++)
	{
		for(int j = i; j > 0; j--)
		{
			if(pVec[j] >= pVec[j - 1])
				break;

			// Swap
			std::swap(pVec[j - 1], pVec[j]);
			if(pParallel1)
				std::swap(pParallel1[j - 1], pParallel1[j]);
			if(pParallel2)
				std::swap(pParallel2[j - 1], pParallel2[j]);
		}
	}
}

// static
void GVec::smallestToFront(double* pVec, int k, int size, double* pParallel1, size_t* pParallel2)
{
	// Use insertion sort if the list is small
	if(size < 7)
	{
		if(k < size)
			GVec_InsertionSort(pVec, size, pParallel1, pParallel2);
		return;
	}
	int beg = 0;
	int end = size - 1;

	// Pick a pivot (using the median of 3 technique)
	double pivA = pVec[0];
	double pivB = pVec[size / 2];
	double pivC = pVec[size - 1];
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
		while(beg < end && pVec[beg] < pivot)
			beg++;
		while(end > beg && pVec[end] > pivot)
			end--;
		if(beg >= end)
			break;
		std::swap(pVec[beg], pVec[end]);
		if(pParallel1)
			std::swap(pParallel1[beg], pParallel1[end]);
		if(pParallel2)
			std::swap(pParallel2[beg], pParallel2[end]);
		beg++;
		end--;
	}

	// Recurse
	if(pVec[beg] < pivot)
		beg++;
	if(k < beg)
		GVec::smallestToFront(pVec, k, beg, pParallel1, pParallel2);
	if(k > beg)
		GVec::smallestToFront(pVec + beg, k - beg, size - beg, pParallel1 ? pParallel1 + beg : NULL, pParallel2 ? pParallel2 + beg : NULL);
}

// static
double GVec::refinePoint(double* pPoint, double* pNeighbor, int dims, double distance, double learningRate, GRand* pRand)
{
	GTEMPBUF(double, buf, dims);
	GVec::copy(buf, pPoint, dims);
	GVec::subtract(buf, pNeighbor, dims);
	double mag = squaredMagnitude(buf, dims);
	GVec::safeNormalize(buf, dims, pRand);
	GVec::multiply(buf, distance, dims);
	GVec::add(buf, pNeighbor, dims);
	GVec::subtract(buf, pPoint, dims);
	GVec::multiply(buf, learningRate, dims);
	GVec::add(pPoint, buf, dims);
	return mag;
}

// static
void GVec::toImage(const double* pVec, GImage* pImage, int width, int height, int channels, double range)
{
	pImage->setSize(width, height);
	unsigned int* pix = pImage->pixels();
	if(channels == 3)
	{
		for(int y = 0; y < height; y++)
		{
			for(int x = 0; x < width; x++)
			{
				int r = ClipChan((int)(*(pVec++) * 255 / range));
				int g = ClipChan((int)(*(pVec++) * 255 / range));
				int b = ClipChan((int)(*(pVec++) * 255 / range));
				*(pix++) = gARGB(0xff, r, g, b);
			}
		}
	}
	else if(channels == 1)
	{
		for(int y = 0; y < height; y++)
		{
			for(int x = 0; x < width; x++)
			{
				int v = MAX(0, MIN(MAX_GRAY_VALUE, (int)(*(pVec++) * MAX_GRAY_VALUE / range)));
				*(pix++) = gGray(v);
			}
		}
	}
	else
		ThrowError("unsupported value for channels");
}

// static
void GVec::fromImage(GImage* pImage, double* pVec, int width, int height, int channels, double range)
{
	unsigned int* pix = pImage->pixels();
	if(channels == 3)
	{
		for(int y = 0; y < height; y++)
		{
			for(int x = 0; x < width; x++)
			{
				*(pVec++) = gRed(*pix) * range / 255;
				*(pVec++) = gGreen(*pix) * range / 255;
				*(pVec++) = gBlue(*pix) * range / 255;
				pix++;
			}
		}
	}
	else if(channels == 1)
	{
		for(int y = 0; y < height; y++)
		{
			for(int x = 0; x < width; x++)
			{
				*(pVec++) = gGray(*pix) * range / MAX_GRAY_VALUE;
				pix++;
			}
		}
	}
	else
		ThrowError("unsupported value for channels");
}

// static
void GVec::capValues(double* pVec, double cap, int dims)
{
	while(dims >= 0)
	{
		*pVec = MIN(*pVec, cap);
		pVec++;
		dims--;
	}
}

// static
void GVec::floorValues(double* pVec, double floor, int dims)
{
	while(dims >= 0)
	{
		*pVec = MAX(*pVec, floor);
		pVec++;
		dims--;
	}
}

#ifndef NO_TEST_CODE
// static
void GVec::test()
{
	GRand prng(0);
	GTEMPBUF(double, v1, 200);
	double* v2 = v1 + 100;
	for(int i = 0; i < 10; i++)
	{
		prng.spherical(v1, 100);
		prng.spherical(v2, 100);
		GVec::subtractComponent(v2, v1, 100);
		GVec::normalize(v2, 100);
		if(ABS(GVec::correlation(v1, v2, 100)) > 1e-4)
			ThrowError("Failed");
		if(ABS(GVec::squaredMagnitude(v1, 100) - 1) > 1e-4)
			ThrowError("Failed");
		if(ABS(GVec::squaredMagnitude(v2, 100) - 1) > 1e-4)
			ThrowError("Failed");
	}
}
#endif // NO_TEST_CODE







// static
void GIndexVec::makeIndexVec(size_t* pVec, size_t size)
{
	for(size_t i = 0; i < size; i++)
	{
		*pVec = i;
		pVec++;
	}
}

// static
void GIndexVec::shuffle(size_t* pVec, size_t size, GRand* pRand)
{
	for(size_t i = size; i > 1; i--)
	{
		size_t r = (size_t)pRand->next(i);
		size_t t = pVec[i - 1];
		pVec[i - 1] = pVec[r];
		pVec[r] = t;
	}
}

// static
void GIndexVec::setAll(size_t* pVec, size_t value, size_t size)
{
	while(size > 0)
	{
		*pVec = value;
		pVec++;
		size--;
	}
}

// static
void GIndexVec::copy(size_t* pDest, const size_t* pSource, size_t nDims)
{
	memcpy(pDest, pSource, sizeof(size_t) * nDims);
}

// static
size_t GIndexVec::maxValue(size_t* pVec, size_t size)
{
	size_t m = *(pVec++);
	size--;
	while(size > 0)
	{
		m = MAX(m, *(pVec++));
		size--;
	}
	return m;
}











GCoordVectorIterator::GCoordVectorIterator(size_t dims, size_t* pRanges)
{
	m_pCoords = NULL;
	reset(dims, pRanges);
}

GCoordVectorIterator::GCoordVectorIterator(vector<size_t>& ranges)
{
	m_pCoords = NULL;
	reset(ranges);
}

GCoordVectorIterator::~GCoordVectorIterator()
{
	delete[] m_pCoords;
}

void GCoordVectorIterator::reset()
{
	memset(m_pCoords, '\0', sizeof(size_t) * m_dims);
	m_sampleShift = 0xffffffff;
}

void GCoordVectorIterator::reset(size_t dims, size_t* pRanges)
{
	m_dims = dims;
	delete[] m_pCoords;
	if(dims > 0)
	{
		m_pCoords = new size_t[2 * dims];
		m_pRanges = m_pCoords + dims;
		if(pRanges)
			memcpy(m_pRanges, pRanges, sizeof(size_t) * dims);
		else
		{
			for(size_t i = 0; i < dims; i++)
				m_pRanges[i] = 1;
		}
	}
	else
	{
		m_pCoords = NULL;
		m_pRanges = NULL;
	}
	reset();
}

void GCoordVectorIterator::reset(vector<size_t>& ranges)
{
	m_dims = ranges.size();
	delete[] m_pCoords;
	if(m_dims > 0)
	{
		m_pCoords = new size_t[2 * m_dims];
		m_pRanges = m_pCoords + m_dims;
		for(size_t i = 0; i < m_dims; i++)
			m_pRanges[i] = ranges[i];
	}
	else
	{
		m_pCoords = NULL;
		m_pRanges = NULL;
	}
	reset();
}

size_t GCoordVectorIterator::coordCount()
{
	size_t n = 1;
	size_t* pR = m_pRanges;
	for(size_t i = 0; i < m_dims; i++)
		n *= (*(pR++));
	return n;
}

bool GCoordVectorIterator::advance()
{
	size_t j;
	for(j = 0; j < m_dims; j++)
	{
		if(++m_pCoords[j] >= m_pRanges[j])
			m_pCoords[j] = 0;
		else
			break;
	}

	// Test if we're done
	if(j >= m_dims)
		return false;
	return true;
}

bool GCoordVectorIterator::advance(size_t steps)
{
	size_t j;
	for(j = 0; j < m_dims; j++)
	{
		size_t t = m_pCoords[j] + steps;
		m_pCoords[j] = t % m_pRanges[j];
		steps = t / m_pRanges[j];
		if(t == 0)
			break;
	}

	// Test if we're done
	if(j >= m_dims)
		return false;
	return true;
}

bool GCoordVectorIterator::advanceSampling()
{
	if(m_sampleShift == 0xffffffff) // if we have not yet computed the step size
	{
		size_t r = m_pRanges[0];
		for(size_t i = 1; i < m_dims; i++)
			r = MAX(r, m_pRanges[i]);
		m_sampleShift = GBits::boundingShift(r);
		m_sampleMask = 0;
	}

	m_pCoords[0] += (1 << (m_sampleShift + (m_sampleMask ? 0 : 1)));
	if(m_pCoords[0] >= m_pRanges[0])
	{
		m_pCoords[0] = 0;
		size_t j = 1;
		for( ; j < m_dims; j++)
		{
			m_pCoords[j] += (1u << m_sampleShift);
			m_sampleMask ^= (1u << j);
			if(m_pCoords[j] < m_pRanges[j])
				break;
			m_pCoords[j] = 0;
			m_sampleMask &= ~(1u << j);
		}
		if(j >= m_dims)
		{
			if(--m_sampleShift == 0xffffffff) // if we're all done
				return false;
		}
		if(m_sampleMask == 0)
		{
			m_pCoords[0] -= (1u << m_sampleShift);
			return advanceSampling();
		}
	}
	return true;
}

size_t* GCoordVectorIterator::current()
{
	return m_pCoords;
}

void GCoordVectorIterator::currentNormalized(double* pCoords)
{
	for(size_t i = 0; i < m_dims; i++)
	{
		*pCoords = ((double)m_pCoords[i] + 0.5) / m_pRanges[i];
		pCoords++;
	}
}

size_t GCoordVectorIterator::currentIndex()
{
	size_t index = 0;
	size_t n = 1;
	for(size_t i = 0; i < m_dims; i++)
	{
		index += n * m_pCoords[i];
		n *= m_pRanges[i];
	}
	return index;
}

void GCoordVectorIterator::setRandom(GRand* pRand)
{
	for(size_t i = 0; i < m_dims; i++)
		m_pCoords[i] = (size_t)pRand->next(m_pRanges[i]);
}

#ifndef NO_TEST_CODE
#define TEST_DIMS 4
// static
void GCoordVectorIterator::test()
{
	size_t r = 11;
	size_t size = 1;
	for(size_t i = 0; i < TEST_DIMS; i++)
		size *= r;
	GBitTable bt(size);
	size_t ranges[TEST_DIMS];
	for(size_t i = 0; i < TEST_DIMS; i++)
		ranges[i] = r;
	GCoordVectorIterator cvi(TEST_DIMS, ranges);
	size_t count = 0;
	while(true)
	{
		size_t index = cvi.currentIndex();
		if(bt.bit(index))
			ThrowError("already got this one");
		bt.set(index);
		count++;
		if(!cvi.advanceSampling())
			break;
	}
	if(count != size)
		ThrowError("didn't get them all");
}
#endif


} // namespace GClasses

