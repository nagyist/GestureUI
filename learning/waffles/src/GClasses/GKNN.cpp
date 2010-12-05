/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#include "GKNN.h"
#include <math.h>
#include "GMacros.h"
#include <stdlib.h>
#include "GTwt.h"
#include "GDistribution.h"
#include "GRand.h"
#include "GHeap.h"
#include "GNeighborFinder.h"
#include "GVec.h"
#include "GHillClimber.h"
#include "GCluster.h"
#include "GBitTable.h"

using namespace GClasses;

namespace GClasses {
class GKnnScaleFactorCritic : public GTargetFunction
{
protected:
	int m_labelDims;
	GKNN* m_pLearner;
	double* m_pAccuracy;

public:
	GKnnScaleFactorCritic(GKNN* pLearner, int featureDims, int labelDims)
	: GTargetFunction(featureDims), m_labelDims(labelDims)
	{
		m_pLearner = pLearner;
		m_pAccuracy = new double[labelDims];
	}

	virtual ~GKnnScaleFactorCritic()
	{
		delete[] m_pAccuracy;
	}

	virtual void initVector(double* pVector)
	{
		GRowDistanceScaled* pMetric = m_pLearner->metric();
		GVec::copy(pVector, pMetric->scaleFactors(), relation()->size());
	}

	virtual bool isStable() { return false; }
	virtual bool isConstrained() { return false; }

protected:
	virtual double computeError(const double* pVector)
	{
		GData* pData = m_pLearner->instances();
		if(!pData)
			return 1e308;
		GKNN temp(m_pLearner->neighborCount(), m_pLearner->getRand());
		temp.enableIncrementalLearning(pData->relation(), m_labelDims, NULL, NULL);
		GVec::copy(temp.metric()->scaleFactors(), pVector, relation()->size());
		return temp.heuristicValidate(pData, m_labelDims, m_pLearner->getRand());
	}
};
}


GKNN::GKNN(int nNeighbors, GRand* pRand)
: GIncrementalLearner(), m_pRand(pRand), m_featureDims(0), m_labelDims(0)
{
	m_eInterpolationMethod = Linear;
	m_pLearner = NULL;
	m_bOwnLearner = false;
	m_nNeighbors = nNeighbors;
	m_pInstances = NULL;
	m_pNeighborFinder = NULL;
	m_pNeighborFinder2 = NULL;
	m_pEvalNeighbors = new size_t[m_nNeighbors + 1];
	m_pEvalDistances = new double[m_nNeighbors + 1];
	m_optimizeScaleFactors = false;
	m_pDistanceMetric = new GRowDistanceScaled();
	m_pValueCounts = NULL;
	m_pCritic = NULL;
	m_pScaleFactorOptimizer = NULL;
	m_dElbowRoom = 0.01;
}

GKNN::GKNN(GTwtNode* pNode, GRand* pRand)
: GIncrementalLearner(pNode), m_pRand(pRand)
{
	m_pNeighborFinder = NULL;
	m_pNeighborFinder2 = NULL;
	m_pCritic = NULL;
	m_pScaleFactorOptimizer = NULL;
	m_pLearner = NULL;
	m_pValueCounts = NULL;
	m_bOwnLearner = false;
	m_eInterpolationMethod = (InterpolationMethod)pNode->field("interpMethod")->asInt();
	m_dElbowRoom = pNode->field("elbowRoom")->asDouble();
	m_optimizeScaleFactors = pNode->field("optimize")->asBool();
	sp_relation pRel = GRelation::fromTwt(pNode->field("relation"));
	m_nNeighbors = (int)pNode->field("neighbors")->asInt();
	m_labelDims = 0;
	int labelDims = (int)pNode->field("labelDims")->asInt();
	m_pDistanceMetric = new GRowDistanceScaled(pNode->field("metric"));
	m_pInstances = NULL;
	m_pEvalNeighbors = new size_t[m_nNeighbors + 1];
	m_pEvalDistances = new double[m_nNeighbors + 1];
	enableIncrementalLearning(pRel, labelDims, NULL, NULL);
	delete(m_pInstances);
	m_pInstances = new GData(pNode->field("instances"));
}

GKNN::~GKNN()
{
	delete(m_pNeighborFinder);
	delete(m_pInstances);
	delete[] m_pEvalNeighbors;
	delete[] m_pEvalDistances;
	delete[] m_pValueCounts;
	delete(m_pScaleFactorOptimizer);
	delete(m_pCritic);
	delete(m_pDistanceMetric);
}

// virtual
GTwtNode* GKNN::toTwt(GTwtDoc* pDoc)
{
	GTwtNode* pNode = baseTwtNode(pDoc, "GKNN");
	pNode->addField(pDoc, "relation", m_pRelation->toTwt(pDoc));
	pNode->addField(pDoc, "neighbors", pDoc->newInt(m_nNeighbors));
	pNode->addField(pDoc, "labelDims", pDoc->newInt(m_labelDims));
	if(m_eInterpolationMethod == Learner)
		ThrowError("Sorry, toTwt is not supported for the \"Learner\" interpolation method");
	pNode->addField(pDoc, "interpMethod", pDoc->newInt(m_eInterpolationMethod));
	pNode->addField(pDoc, "optimize", pDoc->newBool(m_optimizeScaleFactors));
	pNode->addField(pDoc, "elbowRoom", pDoc->newDouble(m_dElbowRoom));
	pNode->addField(pDoc, "metric", m_pDistanceMetric->toTwt(pDoc));
	pNode->addField(pDoc, "instances", m_pInstances->toTwt(pDoc));
	return pNode;
}

// virtual
int GKNN::featureDims()
{
	if(m_labelDims < 1)
		ThrowError("not yet trained");
	return m_pInstances->cols() - m_labelDims;
}

// virtual
int GKNN::labelDims()
{
	if(m_labelDims < 1)
		ThrowError("not yet trained");
	return m_labelDims;
}

void GKNN::setInterpolationMethod(InterpolationMethod eMethod)
{
	if(eMethod == Learner)
		ThrowError("You should call SetInterpolationLearner instead");
	m_eInterpolationMethod = eMethod;
}

void GKNN::setInterpolationLearner(GSupervisedLearner* pLearner, bool bOwnLearner)
{
	if(m_bOwnLearner)
		delete(m_pLearner);
//	GAssert(pLearner->relation() == GetRelation()); // Expected the learner to be constructed with the same relation as this object
	m_pLearner = pLearner;
	m_eInterpolationMethod = Learner;
	m_bOwnLearner = bOwnLearner;
}

size_t GKNN::addVector(const double* pIn, const double* pOut)
{
	if(!m_pInstances)
		m_pInstances = new GData(m_pRelation);

	// Make a copy of the vector
	double* pVec = new double[m_featureDims + m_labelDims];
	GVec::copy(pVec, pIn, m_featureDims);
	GVec::copy(pVec + m_featureDims, pOut, m_labelDims);

	// Add it to the known instances
	size_t index;
	if(m_pNeighborFinder)
		index = m_pNeighborFinder->addVector(pVec); // GNeighborFinder has a pointer to m_pInstances, and it will add pVec to m_pInstances when we make this call
	else
	{
		index = m_pInstances->rows();
		m_pInstances->takeRow(pVec);
	}
	return index;
}

void GKNN::setOptimizeScaleFactors(bool b)
{
	m_optimizeScaleFactors = b;
}

// virtual
void GKNN::enableIncrementalLearning(sp_relation& pRelation, int labelDims, double* pMins, double* pRanges)
{
	clear();
	m_pRelation = pRelation;
	m_labelDims = labelDims;
	m_featureDims = pRelation->size() - m_labelDims;
	m_pInstances = new GData(m_pRelation);
	GMixedRelation* pMixedRel = new GMixedRelation();
	pMixedRel->addAttrs(pRelation.get(), 0, pRelation->size() - labelDims);
	sp_relation pRelInputs = pMixedRel;
	m_pDistanceMetric->init(pRelInputs);

	// Allocate some other buffers
	int maxOutputValueCount = 0;
	for(int n = 0; n < m_labelDims; n++)
		maxOutputValueCount = MAX(maxOutputValueCount, pRelation->valueCount(m_featureDims + n));
	m_pValueCounts = new double[maxOutputValueCount];

	// Scale factor optimization
	if(m_optimizeScaleFactors)
	{
		m_pCritic = new GKnnScaleFactorCritic(this, m_labelDims, m_featureDims);
		m_pScaleFactorOptimizer = new GMomentumGreedySearch(m_pCritic);
	}
	else
	{
		m_pCritic = NULL;
		m_pScaleFactorOptimizer = NULL;
	}
}

// virtual
void GKNN::trainIncremental(const double* pIn, const double* pOut)
{
	// Make a copy of the vector
	size_t index = addVector(pIn, pOut);

	// Delete the closest neighbor if the (k+1)th neighbor is closer than the specified threshold
	if(!m_pNeighborFinder2)
	{
		m_pNeighborFinder2 = new GKdTree(m_pInstances, m_labelDims, m_nNeighbors + 1, m_pDistanceMetric, false);
		return;
	}
	m_pNeighborFinder2->neighbors(m_pEvalNeighbors, m_pEvalDistances, index);
	m_pNeighborFinder2->sortNeighbors(m_pEvalNeighbors, m_pEvalDistances);
	if(m_pEvalNeighbors[m_nNeighbors] >= 0 && m_pEvalDistances[m_nNeighbors] < m_dElbowRoom)
	{
		double* pClosest = m_pNeighborFinder->releaseVector(m_pEvalNeighbors[0]);
		delete[] pClosest;
	}

	// Learn how to scale the attributes
	if(m_pScaleFactorOptimizer && m_pInstances->rows() > 50)
	{
		m_pScaleFactorOptimizer->iterate();
		GVec::copy(m_pDistanceMetric->scaleFactors(), m_pScaleFactorOptimizer->currentVector(), m_featureDims);
	}
}

void GKNN::train(GData* pData, int labelDims)
{
	enableIncrementalLearning(pData->relation(), labelDims, NULL, NULL);
	m_pInstances->reserve(pData->rows());
	size_t nVectors = pData->rows();
	double* pVector;
	for(size_t n = 0; n < nVectors; n++)
	{
		pVector = pData->row(n);
		addVector(pVector, pVector + m_featureDims);
	}

	// Give each attribute an equal chance by scaling out the deviation
	double m, v;
	double* pScaleFactors = m_pDistanceMetric->scaleFactors();
	for(int i = 0; i < m_featureDims; i++)
	{
		if(pData->relation()->valueCount(i) == 0)
		{
			m = m_pInstances->mean(i);
			v = m_pInstances->variance(i, m);
			pScaleFactors[i] = 1.0 / (2.0 * sqrt(MAX(v, 1e-12)));
		}
		else
			pScaleFactors[i] = 1.0;
	}

	// Learn to scale the attributes
	if(m_pScaleFactorOptimizer)
	{
		if(!m_pNeighborFinder)
		{
			if(!m_pInstances)
				m_pInstances = new GData(m_pRelation);
			m_pNeighborFinder = new GKdTree(m_pInstances, m_labelDims, m_nNeighbors, m_pDistanceMetric, false);
		}
		int j;
		for(j = 0; j < 5; j++)
		{
			for(int i = 0; i < 20; i++)
				m_pScaleFactorOptimizer->iterate();
			m_pNeighborFinder->reoptimize();
		}
		GVec::copy(pScaleFactors, m_pScaleFactorOptimizer->currentVector(), m_featureDims);
	}
}

// virtual
void GKNN::trainSparse(GSparseMatrix* pData, int labelDims)
{
	ThrowError("Sorry, trainSparse is not implemented yet in GKNN");
}

void GKNN::findNeighbors(const double* pVector)
{
	if(!m_pNeighborFinder)
	{
		if(!m_pInstances)
			m_pInstances = new GData(m_pRelation);
		m_pNeighborFinder = new GKdTree(m_pInstances, m_labelDims, m_nNeighbors, m_pDistanceMetric, false);
	}
	m_pNeighborFinder->neighbors(m_pEvalNeighbors, m_pEvalDistances, pVector);
}

void GKNN::interpolateMean(const double* pIn, GPrediction* pOut, double* pOut2)
{
	int i, j, index;
	double* pVectorNeighbor;
	double mean;
	for(i = 0; i < m_labelDims; i++)
	{
		index = m_featureDims + i;
		if(m_pRelation->valueCount(index) == 0)
		{
			double dSum = 0;
			double dSumOfSquares = 0;
			int count = 0;
			for(j = 0; j < m_nNeighbors; j++)
			{
				size_t k = m_pEvalNeighbors[j];
				if(k < m_pInstances->rows())
				{
					pVectorNeighbor = m_pInstances->row(k);
					dSum += pVectorNeighbor[index];
					dSumOfSquares += (pVectorNeighbor[index] * pVectorNeighbor[index]);
					count++;
				}
			}
			if(pOut)
			{
				if(count > 0)
				{
					mean = dSum / count;
					pOut[i].makeNormal()->setMeanAndVariance(mean, dSumOfSquares / count - (mean * mean));
				}
				else
					pOut[i].makeNormal()->setMeanAndVariance(0, 1);
			}
			if(pOut2)
			{
				if(count > 0)
					pOut2[i] = dSum / count;
				else
					pOut2[i] = 0;
			}
		}
		else
		{
			int nValueCount = m_pRelation->valueCount(index);
			GVec::setAll(m_pValueCounts, 0.0, nValueCount);
			for(j = 0; j < m_nNeighbors; j++)
			{
				size_t k = m_pEvalNeighbors[j];
				if(k < m_pInstances->rows())
				{
					pVectorNeighbor = m_pInstances->row(k);
					m_pValueCounts[(int)pVectorNeighbor[index]]++;
				}
			}
			if(pOut)
				pOut[i].makeCategorical()->setValues(nValueCount, m_pValueCounts);
			if(pOut2)
				pOut2[i] = GVec::indexOfMax(m_pValueCounts, nValueCount, m_pRand);
		}
	}
}

void GKNN::interpolateLinear(const double* pIn, GPrediction* pOut, double* pOut2)
{
	int i, j, index, val;
	double d;
	double* pVectorNeighbor;
	for(i = 0; i < m_labelDims; i++)
	{
		index = m_featureDims + i;
		if(m_pRelation->valueCount(index) == 0)
		{
			double dSum = 0;
			double dSumOfSquares = 0;
			double dTot = 0;
			for(j = 0; j < m_nNeighbors; j++)
			{
				size_t k = m_pEvalNeighbors[j];
				if(k < m_pInstances->rows())
				{
					pVectorNeighbor = m_pInstances->row(k);
					if(pVectorNeighbor[index] == UNKNOWN_REAL_VALUE)
						ThrowError("GKNN doesn't support unknown output values");
					d = m_pEvalDistances[j];
					d = 1.0 / MAX(sqrt(d), 1e-9);
					dTot += d;
					d *= pVectorNeighbor[index];
					dSum += d;
					d *= pVectorNeighbor[index];
					dSumOfSquares += d;
				}
			}
			if(pOut)
			{
				if(dTot > 0)
				{
					d = dSum / dTot;
					pOut[i].makeNormal()->setMeanAndVariance(d, dSumOfSquares / dTot - (d * d));
				}
				else
					pOut[i].makeNormal()->setMeanAndVariance(0, 1);
			}
			if(pOut2)
			{
				if(dTot > 0)
					pOut2[i] = dSum / dTot;
				else
					pOut2[i] = 0;
			}
		}
		else
		{
			int nValueCount = m_pRelation->valueCount(index);
			GVec::setAll(m_pValueCounts, 0.0, nValueCount);
			double dSumWeight = 0;
			for(j = 0; j < m_nNeighbors; j++)
			{
				size_t k = m_pEvalNeighbors[j];
				if(k < m_pInstances->rows())
				{
					pVectorNeighbor = m_pInstances->row(k);
					d = m_pEvalDistances[j];
					d = 1 / MAX(d, 1e-9); // to be truly "linear", we should use sqrt(d) instead of d, but this is faster to compute and arguably better for nominal values anyway
					val = (int)pVectorNeighbor[index];
					if(val < 0 || val >= nValueCount)
					{
						GAssert(val == UNKNOWN_DISCRETE_VALUE);
						ThrowError("GKNN doesn't support unknown output values");
					}
					m_pValueCounts[(int)pVectorNeighbor[index]] += d;
					dSumWeight += d;
				}
			}
			if(pOut)
				pOut[i].makeCategorical()->setValues(nValueCount, m_pValueCounts);
			if(pOut2)
				pOut2[i] = GVec::indexOfMax(m_pValueCounts, nValueCount, m_pRand);
		}
	}
}

void GKNN::interpolateLearner(const double* pIn, GPrediction* pOut, double* pOut2)
{
	GAssert(m_pLearner); // no learner is set
	int i;
	GHeap heap(1000);
	GData data(m_pRelation, &heap);
	GReleaseDataHolder hData(&data);
	data.reserve(m_nNeighbors);
	for(i = 0; i < m_nNeighbors; i++)
	{
		size_t nNeighbor = m_pEvalNeighbors[i];
		if(nNeighbor < m_pInstances->rows())
			data.takeRow(m_pInstances->row(nNeighbor));
	}
	m_pLearner->train(&data, m_labelDims);
	if(pOut)
		m_pLearner->predictDistribution(pIn, pOut);
	if(pOut2)
		m_pLearner->predict(pIn, pOut2);
}

// virtual
void GKNN::predictDistribution(const double* pIn, GPrediction* pOut)
{
	findNeighbors(pIn);
	switch(m_eInterpolationMethod)
	{
		case Linear: interpolateLinear(pIn, pOut, NULL); break;
		case Mean: interpolateMean(pIn, pOut, NULL); break;
		case Learner: interpolateLearner(pIn, pOut, NULL); break;
		default:
			GAssert(false); // unexpected enumeration
			break;
	}
}

// virtual
void GKNN::predict(const double* pIn, double* pOut)
{
	findNeighbors(pIn);
	switch(m_eInterpolationMethod)
	{
		case Linear: interpolateLinear(pIn, NULL, pOut); break;
		case Mean: interpolateMean(pIn, NULL, pOut); break;
		case Learner: interpolateLearner(pIn, NULL, pOut); break;
		default:
			GAssert(false); // unexpected enumeration
			break;
	}
}

// virtual
void GKNN::clear()
{
	delete(m_pNeighborFinder);
	m_pNeighborFinder = NULL;
	if(m_pInstances)
		delete(m_pInstances);
	m_pInstances = NULL;
	delete(m_pScaleFactorOptimizer);
	m_pScaleFactorOptimizer = NULL;
	delete(m_pCritic);
	m_pCritic = NULL;
}

#ifndef NO_TEST_CODE
//static
void GKNN::test()
{
	GRand prng(0);
	GKNN knn(3, &prng);
	knn.basicTest(0.68, &prng);
}
#endif

// ---------------------------------------------------------------------------------------

namespace GClasses {

/*
class GAutoScaleInstance
{
protected:
	int m_nDims;
	double* m_pVec;
	double* m_pScales;

public:
	GAutoScaleInstance()
	{
		m_pVec = NULL;
	}

	~GAutoScaleInstance()
	{
		delete[] m_pVec;
	}

	void fromTwt(GTwtNode* pNode)
	{
		GTwtNode* pVec = pNode->field("vec");
		m_nDims = (int)pVec->itemCount();
		m_pVec = new double[2 * m_nDims];
		m_pScales = m_pVec + m_nDims;
		GTwtNode* pScales = pNode->field("scales");
		if((int)pScales->itemCount() != m_nDims)
			ThrowError("Wrong size");
		GVec::fromTwt(m_pVec, m_nDims, pVec);
		GVec::fromTwt(m_pScales, m_nDims, pScales);
	}

	GTwtNode* toTwt(GTwtDoc* pDoc)
	{
		GTwtNode* pNode = pDoc->newObj();
		pNode->addField(pDoc, "vec", GVec::toTwt(pDoc, m_pVec, m_nDims));
		pNode->addField(pDoc, "scales", GVec::toTwt(pDoc, m_pScales, m_nDims));
		return pNode;
	}

	double ComputeLikelihood(const double* pParams)
	{
		double dSum = 0;
		double d;
		int i;
		for(i = 0; i < m_nDims; i++)
		{
			d = (pParams[i] - m_pVec[i]) * m_pScales[i];
			dSum += (d * d);
		}
		return exp(-0.5 * dSum);
	}

	void EvalContinuous(int nOutput, const double* pIn, double* pWeight, double* pSum, double* pSumOfSquares)
	{
		double dWeight = ComputeLikelihood(pIn);
		(*pWeight) += dWeight;
		(*pSum) += dWeight * m_pVec[m_nDims + nOutput];
		(*pSumOfSquares) += dWeight * m_pVec[m_nDims + nOutput] * m_pVec[m_nDims + nOutput];
	}

	void EvalNominal(int nOutput, const double* pIn, double* pValues)
	{
		double dWeight = ComputeLikelihood(pIn);
		int nValue = (int)m_pVec[m_nDims + nOutput];
		pValues[nValue] += dWeight;
	}

	void Set(const double* pVars)
	{
		GVec::copy(m_pScales, pVars, m_nDims);
	}

	void GetVars(double* pOut)
	{
		GVec::copy(pOut, m_pScales, m_nDims);
	}

	void Init(GKernelInstanceLearner* pLearner, double* pVector, double* pInitialScales, int nInputCount, int nOutputCount)
	{
		m_nDims = nInputCount;
		m_pVec = new double[m_nDims + nOutputCount + m_nDims];
		m_pScales = m_pVec + m_nDims + nOutputCount;
		GVec::copy(m_pVec, pVector, m_nDims + nOutputCount);
		GVec::copy(m_pScales, pInitialScales, m_nDims);
	}
};

class GAutoScaleInstanceCritic : public GTargetFunction
{
protected:
	double* m_pStartState;
	GKernelInstanceLearner* m_pLearner;
	GData* m_pData;
	int m_featureDims;
	int m_labelDims;
	GAgglomerativeClusterer* m_pClust;
	int m_nCluster;

public:
	GAutoScaleInstanceCritic(GKernelInstanceLearner* pLearner, double* pStartState, GData* pData, int nInputCount, int nOutputCount, GAgglomerativeClusterer* pClust, int nCluster)
	: GTargetFunction(nInputCount)
	{
		m_pLearner = pLearner;
		m_pStartState = pStartState;
		m_pData = pData;
		m_labelDims = nOutputCount;
		m_featureDims = nInputCount;
		m_pClust = pClust;
		m_nCluster = nCluster;
	}

	virtual ~GAutoScaleInstanceCritic()
	{
	}

	virtual void initVector(double* pVector)
	{
		GVec::copy(pVector, m_pStartState, relation()->size());
	}

	virtual bool isStable() { return true; }
	virtual bool isConstrained() { return false; }

protected:
	virtual double computeError(const double* pVector)
	{
		if(m_pClust)
			m_pLearner->SetClusterVars(pVector, m_pData, m_pClust, m_nCluster);
		else
			m_pLearner->SetGlobalVars(pVector, m_pData);
		double* pVec;
#ifdef WIN32
		GPrediction* out = new GPrediction[m_labelDims];
		ArrayHolder<GPrediction> hOut(out);
#else
		GPrediction out[m_labelDims];
#endif
		double dSumError = 0;
		double d;
		for(size_t i = 0; i < m_pData->rows(); i++)
		{
			pVec = m_pData->row(i);
			m_pLearner->HalfEval(pVec, out, m_pData);
			for(int j = 0; j < m_labelDims; j++)
			{
				if(out[j].isContinuous())
					d = out[j].mode() - pVec[m_featureDims + j];
				else
					d = 1.0 - out[j].asCategorical()->likelihood((int)pVec[m_featureDims + j]);
				dSumError += (d * d);
			}
		}
		return dSumError;
	}
};
}


GKernelInstanceLearner::GKernelInstanceLearner(GRand* pRand, int labelDims)
: GSupervisedLearner(), m_labelDims(0)
{
	m_bOptimizeGlobally = true;
	m_bOptimizeLocally = false;
	m_pSet1 = NULL;
	m_pSet2 = NULL;
	m_pRand = pRand;
	m_pInstances = NULL;
	m_nInstanceCount = 0;
}

GKernelInstanceLearner::GKernelInstanceLearner(GTwtNode* pNode, GRand* pRand)
 : GSupervisedLearner(pNode), m_pRand(pRand)
{
	m_pRelation = GRelation::fromTwt(pNode->field("relation"));
	m_bOptimizeGlobally = pNode->field("global")->asBool();
	m_bOptimizeLocally = pNode->field("local")->asBool();
	GTwtNode* pInstances = pNode->field("instances");
	m_nInstanceCount = pInstances->itemCount();
	m_pInstances = new GAutoScaleInstance[m_nInstanceCount];
	for(size_t i = 0; i < m_nInstanceCount; i++)
		m_pInstances[i].fromTwt(pInstances->item(i));
	m_pSet1 = NULL;
	m_pSet2 = NULL;
}

// virtual
GKernelInstanceLearner::~GKernelInstanceLearner()
{
	clear();
}

// virtual
GTwtNode* GKernelInstanceLearner::toTwt(GTwtDoc* pDoc)
{
	GTwtNode* pNode = baseTwtNode(pDoc, "GKernelInstanceLearner");
	pNode->addField(pDoc, "relation", m_pRelation->toTwt(pDoc));
	pNode->addField(pDoc, "global", pDoc->newBool(m_bOptimizeGlobally));
	pNode->addField(pDoc, "local", pDoc->newBool(m_bOptimizeLocally));
	GTwtNode* pInstances = pNode->addField(pDoc, "instances", pDoc->newList(m_nInstanceCount));
	for(size_t i = 0; i < m_nInstanceCount; i++)
		pInstances->setItem(i, m_pInstances[i].toTwt(pDoc));
	return pNode;
}

// virtual
void GKernelInstanceLearner::train(GData* pData, int labelDims)
{
	m_labelDims = labelDims;
	clear();
	int featureDims = pData->cols() - labelDims;
	if(!pData->relation()->areContinuous(0, featureDims))
		ThrowError("GKernelInstanceLearner doesn't support nominal features. You should filter with the NominalToCat transform to convert nominals to reals.");
	m_pRelation = pData->relation();

	// Set the mean and initial covariance for each instance
	double* pInitialSingle = new double[featureDims];
	ArrayHolder<double> hInitialSingle(pInitialSingle);
	double mean, var;
	for(int i = 0; i < featureDims; i++)
	{
		mean = pData->mean(i);
		var = pData->variance(i, mean);
		pInitialSingle[i] = 1.0 / MAX(1e-6, sqrt(var));
	}
	m_nInstanceCount = pData->rows();
	m_pInstances = new GAutoScaleInstance[m_nInstanceCount];
	for(size_t i = 0; i < m_nInstanceCount; i++)
		m_pInstances[i].Init(this, pData->row(i), pInitialSingle, featureDims, labelDims);
	if(m_bOptimizeLocally)
		m_bOptimizeGlobally = true;
	if(!m_bOptimizeGlobally)
		return;

	// Cluster the data
	int nClusterCount = 100;
	GAgglomerativeClusterer clust(nClusterCount);
	if(m_bOptimizeLocally)
		clust.cluster(pData);

	// Split the data into two parts
	pData->shuffle(m_pRand);
	size_t nSize2 = pData->rows() / 2;
	GData dataOther(pData->relation(), pData->heap());
	pData->splitBySize(&dataOther, nSize2);
	m_pSet1 = pData;
	m_pSet2 = &dataOther;
	size_t nSize1 = m_pSet1->rows();
	for(size_t i = 0; i < nSize1; i++)
		m_pInstances[i].Init(this, m_pSet1->row(i), pInitialSingle, featureDims, labelDims);
	for(size_t i = 0; i < nSize2; i++)
		m_pInstances[nSize1 + i].Init(this, m_pSet2->row(i), pInitialSingle, featureDims, labelDims);

	// Find global scale parameters
	{
		GAutoScaleInstanceCritic critic1(this, pInitialSingle, m_pSet2, featureDims, labelDims, NULL, 0);
		GAutoScaleInstanceCritic critic2(this, pInitialSingle, m_pSet1, featureDims, labelDims, NULL, 0);
		GMomentumGreedySearch search1(&critic1);
		GMomentumGreedySearch search2(&critic2);
		search1.setAllStepSizes(.1);
		search2.setAllStepSizes(.1);
		for(int j = 0; j < 10; j++)
		{
			int nIters = 50;
			for(int i = 0; i < nIters; i++)
				search1.iterate();
			for(int i = 0; i < nIters; i++)
				search2.iterate();
		}
	}

	// Find local scale parameters
	if(m_bOptimizeLocally)
	{
		int nCluster, nPass;
		for(nPass = 0; nPass < 4; nPass++)
		{
			for(nCluster = 0; nCluster < nClusterCount; nCluster++)
			{
				// Make sure there's enough data points in this cluster to do something useful
				int count = 0;
				for(size_t i = 0; i < pData->rows(); i++)
				{
					if(clust.whichCluster((int)i) == nCluster)
						count++;
				}
				if(count < 12)
					continue;

				// Find scale parameters for the current cluster
				GAutoScaleInstanceCritic critic1(this, pInitialSingle, m_pSet2, featureDims, labelDims, &clust, nCluster);
				GAutoScaleInstanceCritic critic2(this, pInitialSingle, m_pSet1, featureDims, labelDims, &clust, nCluster);
				GMomentumGreedySearch search1(&critic1);
				GMomentumGreedySearch search2(&critic2);
				search1.setAllStepSizes(.01);
				search2.setAllStepSizes(.01);
				int j;
				for(j = 0; j < 6; j++)
				{
					int nIters = 25;
					for(int i = 0; i < nIters; i++)
						search1.iterate();
					for(int i = 0; i < nIters; i++)
						search2.iterate();
				}
			}
		}
	}

	pData->mergeVert(&dataOther);
	m_pSet1 = NULL;
	m_pSet2 = NULL;
}

void GKernelInstanceLearner::SetGlobalVars(const double* pVector, GData* pValidationSet)
{
	if(pValidationSet == m_pSet2)
	{
		size_t nEnd = m_pSet1->rows();
		for(size_t i = 0; i < nEnd; i++)
			m_pInstances[i].Set(pVector);
	}
	else
	{
		GAssert(pValidationSet == m_pSet1); // unrecognized set
		size_t nStart = m_pSet1->rows();
		size_t nEnd = m_nInstanceCount - nStart;
		for(size_t i = 0; i < nEnd; i++)
			m_pInstances[nStart + i].Set(pVector);
	}
}

void GKernelInstanceLearner::SetClusterVars(const double* pVector, GData* pValidationSet, GAgglomerativeClusterer* pClust, int nCluster)
{
	if(pValidationSet == m_pSet2)
	{
		size_t nEnd = m_pSet1->rows();
		for(size_t i = 0; i < nEnd; i++)
		{
			if(pClust->whichCluster((int)i) == nCluster)
				m_pInstances[i].Set(pVector);
		}
	}
	else
	{
		GAssert(pValidationSet == m_pSet1); // unrecognized set
		size_t nStart = m_pSet1->rows();
		size_t nEnd = m_nInstanceCount - nStart;
		for(size_t i = 0; i < nEnd; i++)
		{
			if(pClust->whichCluster((int)nStart + (int)i) == nCluster)
				m_pInstances[nStart + i].Set(pVector);
		}
	}
}

// virtual
void GKernelInstanceLearner::predictDistribution(const double* pIn, GPrediction* pOut)
{
	int i;
	int nValues;
	int featureDims = m_pRelation->size() - m_labelDims;
	for(i = 0; i < m_labelDims; i++)
	{
		nValues = m_pRelation->valueCount(featureDims + i);
		if(nValues == 0)
		{
			double dWeight = 0;
			double dSum = 0;
			double dSumOfSquares = 0;
			for(size_t j = 0; j < m_nInstanceCount; j++)
				m_pInstances[j].EvalContinuous(i, pIn, &dWeight, &dSum, &dSumOfSquares);
			if(dWeight > 0)
			{
				double mean = dSum / dWeight;
				double var = dSumOfSquares / dWeight - (mean * mean);
				pOut[i].makeNormal()->setMeanAndVariance(mean, var);
			}
			else
				pOut[i].makeNormal()->setMeanAndVariance(0, 1e50);
		}
		else
		{
			GCategoricalDistribution* pCat = pOut[i].makeCategorical();
			double* pValues = pCat->values(nValues);
			GVec::setAll(pValues, 0.0, nValues);
			for(size_t j = 0; j < m_nInstanceCount; j++)
				m_pInstances[j].EvalNominal(i, pIn, pValues);
			pCat->normalize();
		}
	}
}

void GKernelInstanceLearner::HalfEval(const double* pIn, GPrediction* pOut, GData* pSet)
{
	int featureDims = m_pRelation->size() - m_labelDims;
	size_t nStart, nEnd;
	if(pSet == m_pSet1)
	{
		nStart = m_pSet1->rows();
		nEnd = m_nInstanceCount;
	}
	else
	{
		GAssert(pSet == m_pSet2); // unexpected set
		nStart = 0;
		nEnd = m_pSet1->rows();
	}
	int nValues;
	for(int i = 0; i < m_labelDims; i++)
	{
		nValues = m_pRelation->valueCount(featureDims + i);
		if(nValues == 0)
		{
			double dWeight = 0;
			double dSum = 0;
			double dSumOfSquares = 0;
			for(size_t j = nStart; j < nEnd; j++)
				m_pInstances[j].EvalContinuous(i, pIn, &dWeight, &dSum, &dSumOfSquares);
			double mean = dSum / dWeight;
			double var = dSumOfSquares / dWeight - (mean * mean);
			pOut[i].makeNormal()->setMeanAndVariance(mean, var);
		}
		else
		{
			GCategoricalDistribution* pCat = pOut[i].makeCategorical();
			double* pValues = pCat->values(nValues);
			GVec::setAll(pValues, 0.0, nValues);
			for(size_t j = nStart; j < nEnd; j++)
				m_pInstances[j].EvalNominal(i, pIn, pValues);
			pCat->normalize();
		}
	}
}

// virtual
void GKernelInstanceLearner::clear()
{
	delete[] m_pInstances;
	m_pInstances = NULL;
	m_nInstanceCount = 0;
}

*/







GNeighborTransducer::GNeighborTransducer(int neighborCount, GRand* pRand)
: GTransducer(), m_friendCount(neighborCount), m_pRand(pRand)
{
	m_intrinsicDims = 0;
	m_alpha = 0;
	m_beta = 0;
	m_prune = false;
}

void GNeighborTransducer::useFriends(int intrinsicDims, double alpha, double beta, bool prune)
{
	m_intrinsicDims = intrinsicDims;
	m_alpha = alpha;
	m_beta = beta;
	m_prune = prune;
}

// virtual
void GNeighborTransducer::transduce(GData* pDataLabeled, GData* pDataUnlabeled, int labelDims)
{
	if(labelDims != 1)
		ThrowError("Only 1 nominal label is supported");
	if(!pDataLabeled->relation()->areNominal(pDataLabeled->relation()->size() - 1, 1))
		ThrowError("Only nominal labels are supported");
	if(!pDataLabeled->relation()->areContinuous(0, pDataLabeled->relation()->size() - 1))
		ThrowError("Only continuous features are supported");
	if(pDataLabeled->cols() != pDataUnlabeled->cols())
		ThrowError("relations don't match");

	// Make a dataset containing all rows
	GData dataAll(pDataLabeled->relation());
	dataAll.reserve(pDataLabeled->rows() + pDataUnlabeled->rows());
	GReleaseDataHolder hDataAll(&dataAll);
	for(size_t i = 0; i < pDataUnlabeled->rows(); i++)
		dataAll.takeRow(pDataUnlabeled->row(i));
	for(size_t i = 0; i < pDataLabeled->rows(); i++)
		dataAll.takeRow(pDataLabeled->row(i));
	int featureDims = pDataLabeled->cols() - labelDims;
	sp_relation pRelInputs = new GUniformRelation(featureDims, 0);
	dataAll.setRelation(pRelInputs);

	// Find friends
	GNeighborFinder* pNF;
	if(m_intrinsicDims == 0)
		pNF = new GNeighborFinderCacheWrapper(new GKdTree(&dataAll, 0, m_friendCount, NULL, true), true);
	else
		pNF = new GManifoldNeighborFinder(
			&dataAll,
			m_friendCount, // littleK
			m_friendCount * 4, // bigK
			m_intrinsicDims, // intrinsicDims
			m_alpha, // alpha
			m_beta, // beta
			false, // prune?
			m_pRand);
	Holder<GNeighborFinder> hNF(pNF);
	GTEMPBUF(size_t, neighbors, m_friendCount);
	int labelValues = pDataLabeled->relation()->valueCount(featureDims);
	GTEMPBUF(double, tallys, labelValues);

	// Label the unlabeled patterns
	GBitTable labeled(pDataUnlabeled->rows());
	GData labelList(3); // pattern index, most likely label, confidence
	labelList.newRows(pDataUnlabeled->rows());
	for(size_t i = 0; i < pDataUnlabeled->rows(); i++)
		labelList.row(i)[0] = i;
	while(labelList.rows() > 0)
	{
		// Compute the most likely label and the confidence for each pattern
		for(size_t i = 0; i < labelList.rows(); i++)
		{
			// Find the most common label
			double* pRow = labelList.row(i);
			size_t index = (size_t)pRow[0];
			pNF->neighbors(neighbors, index);
			GVec::setAll(tallys, 0.0, labelValues);
			for(int j = 0; j < m_friendCount; j++)
			{
				if(neighbors[j] >= dataAll.rows())
					continue;
				double* pFriend = dataAll.row(neighbors[j]);
				if(neighbors[j] >= pDataUnlabeled->rows())
				{
					if((int)pFriend[featureDims] >= 0 && (int)pFriend[featureDims] < labelValues)
						tallys[(int)pFriend[featureDims]] += 1.0;
				}
				else if(labeled.bit(neighbors[j]))
				{
					if((int)pFriend[featureDims] >= 0 && (int)pFriend[featureDims] < labelValues)
						tallys[(int)pFriend[featureDims]] += 0.6;
				}
			}
			int label = GVec::indexOfMax(tallys, labelValues, m_pRand);
			double conf = tallys[label];

			// Penalize for dissenting votes
			for(int j = 0; j < m_friendCount; j++)
			{
				if(neighbors[j] >= dataAll.rows())
					continue;
				double* pFriend = dataAll.row(neighbors[j]);
				if(neighbors[j] >= pDataUnlabeled->rows())
				{
					if((int)pFriend[featureDims] != label)
						conf *= 0.5;
				}
				else if(labeled.bit(neighbors[j]))
				{
					if((int)pFriend[featureDims] != label)
						conf *= 0.8;
				}
			}
			pRow[1] = label;
			pRow[2] = conf;
		}
		labelList.sort(2);

		// Assign the labels to the patterns we are most confident about
		size_t maxCount = MAX((size_t)5, pDataLabeled->rows() / 5);
		size_t count = 0;
		for(size_t i = labelList.rows() - 1; i < labelList.rows(); i--)
		{
			double* pRow = labelList.row(i);
			size_t index = (size_t)pRow[0];
			int label = (int)pRow[1];
			pDataUnlabeled->row(index)[featureDims] = label;
			labeled.set(index);
			labelList.deleteRow(i);
			if(count >= maxCount)
				break;
			count++;
		}
	}
}









GInstanceTable::GInstanceTable(int dims, size_t* pDims, GRand* pRand)
: GIncrementalLearner(), m_dims(dims), m_pRand(pRand)
{
	m_pDims = new size_t[dims];
	memcpy(m_pDims, pDims, sizeof(size_t) * dims);
	m_product = 1;
	m_pScales = new size_t[dims];
	for(int i = 0; i < dims; i++)
	{
		m_pScales[i] = m_product;
		m_product *= pDims[i];
		m_pDims[i] = pDims[i];
	}
	m_pTable = new double[m_product];
	clear();
}

// virtual
GInstanceTable::~GInstanceTable()
{
	delete[] m_pDims;
	delete[] m_pScales;
	delete[] m_pTable;
}

// virtual
GTwtNode* GInstanceTable::toTwt(GTwtDoc* pDoc)
{
	ThrowError("not implemented yet");
	return NULL;
}

// virtual
void GInstanceTable::trainSparse(GSparseMatrix* pData, int labelDims)
{
	ThrowError("Sorry, trainSparse is not implemented yet in GInstanceTable");
}

// virtual
void GInstanceTable::train(GData* pData, int labelDims)
{
	enableIncrementalLearning(pData->relation(), labelDims, NULL, NULL);
	int dims = pData->relation()->size() - 1;
	for(size_t i = 0; i < pData->rows(); i++)
	{
		double* pRow = pData->row(i);
		trainIncremental(pRow, pRow + dims);
	}
}

// virtual
void GInstanceTable::predictDistribution(const double* pIn, GPrediction* pOut)
{
	size_t pos = 0;
	for(int i = 0; i < m_dims; i++)
	{
		size_t n = (size_t)floor(pIn[i] + 0.5);
		if(n >= m_pDims[i])
			ThrowError("dim=", gformat(i), ", index=", gformat(pIn[i]), ", out of range. Expected >= 0 and < ", gformat(m_pDims[i]));
		pos += n * m_pScales[i];
	}
	int values = m_pRelation->valueCount(m_dims);
	if(values == 0)
	{
		GNormalDistribution* pNorm = pOut->makeNormal();
		pNorm->setMeanAndVariance(m_pTable[pos], 1.0);
	}
	else
	{
		GCategoricalDistribution* pCat = pOut->makeCategorical();
		pCat->setSpike(values, (int)m_pTable[pos], 1);
	}
}

// virtual
void GInstanceTable::predict(const double* pIn, double* pOut)
{
	size_t pos = 0;
	for(int i = 0; i < m_dims; i++)
	{
		size_t n = (size_t)floor(pIn[i] + 0.5);
		if(n >= m_pDims[i])
			ThrowError("dim=", gformat(i), ", index=", gformat(pIn[i]), ", out of range. Expected >= 0 and < ", gformat(m_pDims[i]));
		pos += n * m_pScales[i];
	}
	*pOut = m_pTable[pos];
}

// virtual
void GInstanceTable::clear()
{
	double* p = m_pTable;
	for(size_t i = 0; i < m_product; i++)
		*(p++) = m_pRand->uniform() * 0.1;
}

// virtual
void GInstanceTable::enableIncrementalLearning(sp_relation& pRelation, int labelDims, double* pMins, double* pRanges)
{
	if(labelDims != 1)
		ThrowError("only 1 label dim is supported");
	if(pRelation->size() - 1 != m_dims)
		ThrowError("mismatching number of attributes");
	m_pRelation = pRelation;
}

// virtual
void GInstanceTable::trainIncremental(const double* pIn, const double* pOut)
{
	size_t pos = 0;
	for(int i = 0; i < m_dims; i++)
	{
		size_t n = (size_t)floor(pIn[i] + 0.5);
		if(n >= m_pDims[i])
			ThrowError("dim=", gformat(i), ", index=", gformat(pIn[i]), ", out of range. Expected >= 0 and < ", gformat(m_pDims[i]));
		pos += n * m_pScales[i];
	}
	m_pTable[pos] = *pOut;
}

}
