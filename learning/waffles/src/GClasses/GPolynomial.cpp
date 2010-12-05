/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#include "GMacros.h"
#include "GPolynomial.h"
#include "GData.h"
#include "GVec.h"
#include <math.h>
#include <stdlib.h>
#include "GDistribution.h"
#include "GMath.h"
#include "GHillClimber.h"
#include "GTwt.h"

namespace GClasses {

class GPolynomialLatticeIterator
{
protected:
	int* m_pCoords;
	int m_nDimensions;
	int m_nControlPoints;
	int m_nSkipDimension;

public:
	GPolynomialLatticeIterator(int* pCoords, int nDimensions, int nControlPoints, int nSkipDimension)
	{
		m_pCoords = pCoords;
		m_nDimensions = nDimensions;
		m_nControlPoints = nControlPoints;
		m_nSkipDimension = nSkipDimension;
		int i;
		for(i = 0; i < nDimensions; i++)
		{
			if(i == nSkipDimension)
				continue;
			m_pCoords[i] = nControlPoints - 1;
		}
	}

	~GPolynomialLatticeIterator()
	{
	}

	bool Advance()
	{
		// Move on to the next point in the lattice
		int i = 0;
		if(i == m_nSkipDimension)
			i++;
		while(true)
		{
			if(m_pCoords[i]-- == 0)
			{
				m_pCoords[i] = m_nControlPoints - 1;
				if(++i == m_nSkipDimension)
					++i;
				if(i < m_nDimensions)
					continue;
				else
					return false;
			}
			return true;
		}
	}
};

// ---------------------------------------------------------------------------

GPolynomial::GPolynomial(int nControlPoints)
: GSupervisedLearner(), m_featureDims(0), m_nControlPoints(nControlPoints), m_nCoefficients(0), m_pCoefficients(NULL)
{
	GAssert(nControlPoints > 0);
}

GPolynomial::GPolynomial(GTwtNode* pNode)
 : GSupervisedLearner(pNode)
{
	m_nControlPoints = (int)pNode->field("controlPoints")->asInt();
	m_nCoefficients = 1;
	m_featureDims = (int)pNode->field("featureDims")->asInt();
	int i = m_featureDims;
	while(i > 0)
	{
		m_nCoefficients *= m_nControlPoints;
		i--;
	}
	m_pCoefficients = new double[m_nCoefficients];
	GVec::fromTwt(m_pCoefficients, m_nCoefficients, pNode->field("coefficients"));
}

GPolynomial::~GPolynomial()
{
	clear();
}

int GPolynomial::featureDims()
{
	return m_featureDims;
}

// virtual
GTwtNode* GPolynomial::toTwt(GTwtDoc* pDoc)
{
	if(m_featureDims == 0)
		ThrowError("train has not been called");
	GTwtNode* pNode = baseTwtNode(pDoc, "GPolynomial");
	pNode->addField(pDoc, "featureDims", pDoc->newInt(m_featureDims));
	pNode->addField(pDoc, "controlPoints", pDoc->newInt(m_nControlPoints));
	pNode->addField(pDoc, "coefficients", GVec::toTwt(pDoc, m_pCoefficients, m_nCoefficients));
	return pNode;
}

int GPolynomial::calcIndex(int* pCoords)
{
	int nIndex = 0;
	int n;
	for(n = m_featureDims - 1; n >= 0; n--)
	{
		nIndex *= m_nControlPoints;
		GAssert(pCoords[n] >= 0 && pCoords[n] < m_nControlPoints); // out of range
		nIndex += pCoords[n];
	}
	return nIndex;
}

double GPolynomial::coefficient(int* pCoords)
{
	if(m_featureDims == 0)
		ThrowError("init has not been called");
	return m_pCoefficients[calcIndex(pCoords)];
}

void GPolynomial::setCoefficient(int* pCoords, double dVal)
{
	if(m_featureDims == 0)
		ThrowError("init has not been called");
	m_pCoefficients[calcIndex(pCoords)] = dVal;
}

class GPolynomialRegressCritic : public GTargetFunction
{
protected:
	GPolynomial* m_pPolynomial;
	GData* m_pData;

public:
	GPolynomialRegressCritic(GPolynomial* pPolynomial, GData* pData)
	: GTargetFunction(pPolynomial->coefficientCount())
	{
		m_pPolynomial = pPolynomial;
		m_pData = pData;
	}

	virtual ~GPolynomialRegressCritic()
	{
	}

	virtual bool isStable() { return true; }
	virtual bool isConstrained() { return false; }

protected:
	virtual void initVector(double* pVector)
	{
		GVec::setAll(pVector, 0.0, relation()->size());
	}

	virtual double computeError(const double* pVector)
	{
		m_pPolynomial->setCoefficients(pVector);
		m_pPolynomial->fromBezierCoefficients();
		GPrediction out;
		double d;
		double dSumSquaredError = 0;
		double* pVec;
		int nInputCount = m_pPolynomial->featureDims();
		for(size_t i = 0; i < m_pData->rows(); i++)
		{
			pVec = m_pData->row(i);
			m_pPolynomial->predictDistribution(pVec, &out);
			if(out.isContinuous())
			{
				if(pVec[nInputCount] != UNKNOWN_REAL_VALUE)
				{
					d = out.asNormal()->mean() - pVec[nInputCount];
					dSumSquaredError += (d * d);
				}
			}
			else
			{
				if((int)pVec[nInputCount] >= 0 && (int)out.asCategorical()->mode() != (int)pVec[nInputCount])
					dSumSquaredError += 1;
			}
		}
		return dSumSquaredError / m_pData->rows();
	}
};

void GPolynomial::setCoefficients(const double* pVector)
{
	if(m_featureDims == 0)
		ThrowError("init has not been called");
	GVec::copy(m_pCoefficients, pVector, m_nCoefficients);
}

// virtual
void GPolynomial::clear()
{
	delete[] m_pCoefficients;
	m_pCoefficients = NULL;
}

void GPolynomial::init(GRelation* pRelation)
{
	m_featureDims = pRelation->size() - 1;
	if(pRelation->valueCount(m_featureDims) > 0)
		ThrowError("Only continuous labels are supported");
	m_nCoefficients = 1;
	int i = m_featureDims;
	while(i > 0)
	{
		m_nCoefficients *= m_nControlPoints;
		i--;
	}
	m_pCoefficients = new double[m_nCoefficients];
	GVec::setAll(m_pCoefficients, 0.0, m_nCoefficients);
}

// virtual
void GPolynomial::train(GData* pData, int labelDims)
{
	if(labelDims != 1)
		ThrowError("Sorry, only 1 label dim is supported");
	init(pData->relation().get());
	GPolynomialRegressCritic critic(this, pData);
	//GStochasticGreedySearch search(&critic);
	GMomentumGreedySearch search(&critic);
	search.searchUntil(100, 30, .01);
	setCoefficients(search.currentVector());
	fromBezierCoefficients();
}

// (Warning: this method relies on the order in which GPolynomialLatticeIterator
// visits coefficients in the lattice)
// virtual
void GPolynomial::predictDistribution(const double* pIn, GPrediction* pOut)
{
	if(m_featureDims == 0)
		ThrowError("init has not been called");
	GTEMPBUF(int, pCoords, m_featureDims);
	GPolynomialLatticeIterator iter(pCoords, m_featureDims, m_nControlPoints, -1);
	int n, i;
	double dSum = 0;
	double dVar;
	int nCoeff;
	for(nCoeff = m_nCoefficients - 1; nCoeff >= 0; nCoeff--)
	{
		dVar = 1;
		for(n = 0; n < m_featureDims; n++)
		{
			for(i = pCoords[n]; i > 0; i--)
				dVar *= pIn[n];
		}
		dVar *= m_pCoefficients[nCoeff];
		dSum += dVar;
		iter.Advance();
	}
	pOut->makeNormal()->setMeanAndVariance(dSum, 1);
}

// virtual
void GPolynomial::predict(const double* pIn, double* pOut)
{
	GPrediction out;
	predictDistribution(pIn, &out);
	GPrediction::predictionArrayToVector(1, &out, pOut);
}

void GPolynomial::toBezierCoefficients()
{
	if(m_featureDims == 0)
		ThrowError("init has not been called");
	// Make Pascal's triangle
	GTEMPBUF(int, pCoords, m_featureDims);
	GTEMPBUF(int, pPascalsTriangle, m_nControlPoints);

	// In each dimensional direction...
	GMath::pascalsTriangle(pPascalsTriangle, m_nControlPoints - 1);
	int n, i, j;
	for(n = 0; n < m_featureDims; n++)
	{
		// Across that dimension...
		for(j = 0; j < m_nControlPoints; j++)
		{
			// Iterate over the entire lattice of coefficients (except in dimension n)
			GPolynomialLatticeIterator iter(pCoords, m_featureDims, m_nControlPoints, n);
			while(true)
			{
				// Divide by the corresponding row of Pascal's triangle
				pCoords[n] = j;
				m_pCoefficients[calcIndex(pCoords)] /= pPascalsTriangle[j];
				if(!iter.Advance())
					break;
			}
		}
	}

	// Forward sum the coefficients
	double d;
	for(i = m_nControlPoints - 1; i >= 1; i--)
	{
		// In each dimensional direction...
		for(n = 0; n < m_featureDims; n++)
		{
			// Subtract the neighbor of lesser-significance from each coefficient
			for(j = i; j < m_nControlPoints; j++)
			{
				// Iterate over the entire lattice of coefficients (except in dimension n)
				GPolynomialLatticeIterator iter(pCoords, m_featureDims, m_nControlPoints, n);
				while(true)
				{
					// Subtract the neighbor of lesser-significance from this coefficient
					pCoords[n] = j - 1;
					d = m_pCoefficients[calcIndex(pCoords)];
					pCoords[n] = j;
					m_pCoefficients[calcIndex(pCoords)] += d;
					if(!iter.Advance())
						break;
				}
			}
		}
	}
}

void GPolynomial::fromBezierCoefficients()
{
	if(m_featureDims == 0)
		ThrowError("init has not been called");
	// Forward difference the coefficients
	GTEMPBUF(int, pCoords, m_featureDims);
	GTEMPBUF(int, pPascalsTriangle, m_nControlPoints);
	int n, i, j;
	double d;
	for(i = 1; i < m_nControlPoints; i++)
	{
		// In each dimensional direction...
		for(n = 0; n < m_featureDims; n++)
		{
			// Subtract the neighbor of lesser-significance from each coefficient
			for(j = m_nControlPoints - 1; j >= i; j--)
			{
				// Iterate over the entire lattice of coefficients (except in dimension n)
				GPolynomialLatticeIterator iter(pCoords, m_featureDims, m_nControlPoints, n);
				while(true)
				{
					// Subtract the neighbor of lesser-significance from this coefficient
					pCoords[n] = j - 1;
					d = m_pCoefficients[calcIndex(pCoords)];
					pCoords[n] = j;
					m_pCoefficients[calcIndex(pCoords)] -= d;
					if(!iter.Advance())
						break;
				}
			}
		}
	}

	// In each dimensional direction...
	GMath::pascalsTriangle(pPascalsTriangle, m_nControlPoints - 1);
	for(n = 0; n < m_featureDims; n++)
	{
		// Across that dimension...
		for(j = 0; j < m_nControlPoints; j++)
		{
			// Iterate over the entire lattice of coefficients (except in dimension n)
			GPolynomialLatticeIterator iter(pCoords, m_featureDims, m_nControlPoints, n);
			while(true)
			{
				// Multiply by the corresponding row of Pascal's triangle
				pCoords[n] = j;
				m_pCoefficients[calcIndex(pCoords)] *= pPascalsTriangle[j];
				if(!iter.Advance())
					break;
			}
		}
	}
}

void GPolynomial::differentiate()
{
	if(m_featureDims == 0)
		ThrowError("init has not been called");
	GTEMPBUF(int, pCoords, m_featureDims);
	int n, j;
	double d;
	for(n = 0; n < m_featureDims; n++)
	{
		// Iterate over the entire lattice of coefficients (except in dimension n)
		GPolynomialLatticeIterator iter(pCoords, m_featureDims, m_nControlPoints, n);
		while(true)
		{
			// Differentiate with respect to the n'th dimension
			for(j = 1; j < m_nControlPoints; j++)
			{
				pCoords[n] = j;
				d = m_pCoefficients[calcIndex(pCoords)];
				pCoords[n] = j - 1;
				m_pCoefficients[calcIndex(pCoords)] = d * j;
			}
			pCoords[n] = m_nControlPoints - 1;
			m_pCoefficients[calcIndex(pCoords)] = 0;
			if(!iter.Advance())
				break;
		}
	}
}

void GPolynomial::integrate()
{
	if(m_featureDims == 0)
		ThrowError("init has not been called");
	GTEMPBUF(int, pCoords, m_featureDims);
	int n, j, index;
	double d;
	for(n = 0; n < m_featureDims; n++)
	{
		// Iterate over the entire lattice of coefficients (except in dimension n)
		GPolynomialLatticeIterator iter(pCoords, m_featureDims, m_nControlPoints, n);
		while(true)
		{
			// Integrate in the n'th dimension
			pCoords[n] = 0;
			m_pCoefficients[calcIndex(pCoords)] = 0;
			for(j = m_nControlPoints - 1; j > 0; j--)
			{
				pCoords[n] = j - 1;
				d = m_pCoefficients[calcIndex(pCoords)];
				pCoords[n] = j;
				index = calcIndex(pCoords);
				GAssert(j < m_nControlPoints - 1 || m_pCoefficients[index] == 0); // There's a non-zero value in a highest-order coefficient. This polynomial, therefore, isn't big enough to hold the integral
				m_pCoefficients[index] = d / j;
			}
			if(!iter.Advance())
				break;
		}
	}
}

void GPolynomial::copy(GPolynomial* pOther)
{
	m_featureDims = pOther->m_featureDims;
	if(controlPointCount() >= pOther->controlPointCount())
		ThrowError("this polynomial must have at least as many control points as pOther");
	if(controlPointCount() > pOther->controlPointCount())
		GVec::setAll(m_pCoefficients, 0.0, m_nCoefficients);
	GTEMPBUF(int, pCoords, m_featureDims);
	GPolynomialLatticeIterator iter(pCoords, m_featureDims, pOther->m_nControlPoints, -1);
	while(true)
	{
		m_pCoefficients[calcIndex(pCoords)] = pOther->m_pCoefficients[pOther->calcIndex(pCoords)];
		if(!iter.Advance())
			break;
	}
}

#ifndef NO_TEST_CODE
// static
void GPolynomial::test()
{
	// This test involves a two-dimensional polynomial with three controll points in each dimension.
	// In other words, there is a 3x3 lattice of control points, so there are 9 total control points.
	// In this case, we arbitrarily use {1,2,3,4,5,6,7,8,9} for the coefficients.
	sp_relation rel;
	rel = new GUniformRelation(3, 0);
	GPolynomial gp(3);
	gp.init(rel.get());
	int degrees[2];
	degrees[0] = 0;
	degrees[1] = 0;
	gp.setCoefficient(degrees, 1);
	degrees[0] = 1;
	degrees[1] = 0;
	gp.setCoefficient(degrees, 2);
	degrees[0] = 2;
	degrees[1] = 0;
	gp.setCoefficient(degrees, 3);
	degrees[0] = 0;
	degrees[1] = 1;
	gp.setCoefficient(degrees, 4);
	degrees[0] = 1;
	degrees[1] = 1;
	gp.setCoefficient(degrees, 5);
	degrees[0] = 2;
	degrees[1] = 1;
	gp.setCoefficient(degrees, 6);
	degrees[0] = 0;
	degrees[1] = 2;
	gp.setCoefficient(degrees, 7);
	degrees[0] = 1;
	degrees[1] = 2;
	gp.setCoefficient(degrees, 8);
	degrees[0] = 2;
	degrees[1] = 2;
	gp.setCoefficient(degrees, 9);
	double vars[2];
	vars[0] = 7;
	vars[1] = 11;
	GPrediction val;
	gp.predictDistribution(vars, &val);
	// 1 + 2 * (7) + 3 * (7 * 7) +
	// 4 * (11) + 5 * (11 * 7) + 6 * (11 * 7 * 7) +
	// 7 * (11 * 11) + 8 * (11 * 11 * 7) + 9 * (11 * 11 * 7 * 7)
	// = 64809
	if(val.asNormal()->mean() != 64809)
		throw "wrong answer";
}
#endif // NO_TEST_CODE

} // namespace GClasses

