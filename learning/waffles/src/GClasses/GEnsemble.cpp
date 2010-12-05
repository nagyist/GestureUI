/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#include "GEnsemble.h"
#include "GVec.h"
#include <stdlib.h>
#include "GDistribution.h"
#include "GTwt.h"
#include "GRand.h"

using namespace GClasses;
using std::vector;

GBag::GBag(GRand* pRand)
: GSupervisedLearner(), m_labelDims(0), m_pRand(pRand)
{
	m_nAccumulatorDims = 0;
	m_pAccumulator = NULL;
	m_pCB = NULL;
	m_pThis = NULL;
}

GBag::GBag(GTwtNode* pNode, GRand* pRand, GLearnerLoader* pLoader)
: GSupervisedLearner(pNode), m_pRand(pRand)
{
	m_pRelation = GRelation::fromTwt(pNode->field("relation"));
	m_labelDims = (int)pNode->field("labelDims")->asInt();
	m_nAccumulatorDims = (int)pNode->field("accum")->asInt();
	m_pAccumulator = new double[m_nAccumulatorDims];
	m_pCB = NULL;
	m_pThis = NULL;
	GTwtNode* pModels = pNode->field("models");
	size_t modelCount = pModels->itemCount();
	for(size_t i = 0; i < modelCount; i++)
		m_models.push_back(pLoader->loadModeler(pModels->item(i), pRand));
}

GBag::~GBag()
{
	for(vector<GSupervisedLearner*>::iterator it = m_models.begin(); it != m_models.end(); it++)
		delete(*it);
	delete[] m_pAccumulator;
}

// virtual
GTwtNode* GBag::toTwt(GTwtDoc* pDoc)
{
	GTwtNode* pNode = baseTwtNode(pDoc, "GBag");
	pNode->addField(pDoc, "relation", m_pRelation->toTwt(pDoc));
	pNode->addField(pDoc, "labelDims", pDoc->newInt(m_labelDims));
	pNode->addField(pDoc, "accum", pDoc->newInt(m_nAccumulatorDims));
	GTwtNode* pModels = pNode->addField(pDoc, "models", pDoc->newList(m_models.size()));
	for(size_t i = 0; i < m_models.size(); i++)
		pModels->setItem((int)i, m_models[i]->toTwt(pDoc));
	return pNode;
}

// virtual
int GBag::featureDims()
{
	if(m_labelDims < 1)
		ThrowError("not yet trained");
	return m_pRelation->size() - m_labelDims;
}

// virtual
int GBag::labelDims()
{
	if(m_labelDims < 1)
		ThrowError("not yet trained");
	return m_labelDims;
}

void GBag::clear()
{
	for(vector<GSupervisedLearner*>::iterator it = m_models.begin(); it != m_models.end(); it++)
		(*it)->clear();
}

void GBag::flush()
{
	for(vector<GSupervisedLearner*>::iterator it = m_models.begin(); it != m_models.end(); it++)
		delete(*it);
	m_models.clear();
}

void GBag::addLearner(GSupervisedLearner* pLearner)
{
	m_models.push_back(pLearner);
}

// virtual
void GBag::train(GData* pData, int labelDims)
{
	m_pRelation = pData->relation();
	m_labelDims = labelDims;

	// Make the accumulator buffer
	int featureDims = pData->cols() - labelDims;
	delete[] m_pAccumulator;
	m_nAccumulatorDims = 0;
	int i, nValues;
	for(i = 0; i < labelDims; i++)
	{
		nValues = pData->relation()->valueCount(featureDims + i);
		if(nValues > 0)
			m_nAccumulatorDims += nValues;
		else
			m_nAccumulatorDims += 2; // mean and variance
	}
	m_pAccumulator = new double[m_nAccumulatorDims];

	// Train all the models
	int nLearnerCount = (int)m_models.size();
	size_t nVectorCount = pData->rows();
	GData drawnData(pData->relation(), pData->heap());
	drawnData.reserve(nVectorCount);
	{
		GReleaseDataHolder hDrawnData(&drawnData);
		int i;
		for(i = 0; i < nLearnerCount; i++)
		{
			if(m_pCB)
				m_pCB(m_pThis, i, nLearnerCount);

			// Randomly draw some data (with replacement)
			for(size_t j = 0; j < nVectorCount; j++)
				drawnData.takeRow(pData->row((size_t)m_pRand->next(nVectorCount)));

			// Train the learner with the drawn data
			m_models[i]->train(&drawnData, labelDims);
			drawnData.releaseAllRows();
		}
		if(m_pCB)
			m_pCB(m_pThis, nLearnerCount, nLearnerCount);
	}
}

void GBag::accumulate(const double* pOut)
{
	int featureDims = m_pRelation->size() - m_labelDims;
	int nDims = 0;
	int i, nValues, nVal;
	double dVal;
	for(i = 0; i < m_labelDims; i++)
	{
		nValues = m_pRelation->valueCount(featureDims + i);
		if(nValues > 0)
		{
			nVal = (int)pOut[i];
			if(nVal >= 0 && nVal < nValues)
				m_pAccumulator[nDims + nVal]++;
			nDims += nValues;
		}
		else
		{
			dVal = pOut[i];
			m_pAccumulator[nDims++] += dVal;
			m_pAccumulator[nDims++] += (dVal * dVal);
		}
	}
	GAssert(nDims == m_nAccumulatorDims); // invalid dim count
}

void GBag::tally(int nCount, GPrediction* pOut)
{
	int featureDims = m_pRelation->size() - m_labelDims;
	int nDims = 0;
	int i, nValues;
	double mean;
	for(i = 0; i < m_labelDims; i++)
	{
		nValues = m_pRelation->valueCount(featureDims + i);
		if(nValues > 0)
		{
			pOut[i].makeCategorical()->setValues(nValues, &m_pAccumulator[nDims]);
			nDims += nValues;
		}
		else
		{
			mean = m_pAccumulator[nDims] / nCount;
			pOut[i].makeNormal()->setMeanAndVariance(mean, m_pAccumulator[nDims + 1] / nCount - (mean * mean));
			nDims += 2;
		}
	}
	GAssert(nDims == m_nAccumulatorDims); // invalid dim count
}

void GBag::tally(int nCount, double* pOut)
{
	int featureDims = m_pRelation->size() - m_labelDims;
	int nDims = 0;
	int i, nValues;
	for(i = 0; i < m_labelDims; i++)
	{
		nValues = m_pRelation->valueCount(featureDims + i);
		if(nValues > 0)
		{
			pOut[i] = (double)GVec::indexOfMax(m_pAccumulator + nDims, nValues, m_pRand);
			nDims += nValues;
		}
		else
		{
			pOut[i] = m_pAccumulator[nDims] / nCount;
			nDims += 2;
		}
	}
	GAssert(nDims == m_nAccumulatorDims); // invalid dim count
}

// virtual
void GBag::predictDistribution(const double* pIn, GPrediction* pOut)
{
	GTEMPBUF(double, pTmp, m_labelDims);
	GVec::setAll(m_pAccumulator, 0.0, m_nAccumulatorDims);
	for(vector<GSupervisedLearner*>::iterator it = m_models.begin(); it != m_models.end(); it++)
	{
		(*it)->predict(pIn, pTmp);
		accumulate(pTmp);
	}
	tally((int)m_models.size(), pOut);
}

// virtual
void GBag::predict(const double* pIn, double* pOut)
{
	GVec::setAll(m_pAccumulator, 0.0, m_nAccumulatorDims);
	for(vector<GSupervisedLearner*>::iterator it = m_models.begin(); it != m_models.end(); it++)
	{
		(*it)->predict(pIn, pOut);
		accumulate(pOut);
	}
	tally((int)m_models.size(), pOut);
}
/*
// virtual
double GBag::CrossValidate(GData* pData, int nFolds, bool bRegression)
{
	// Split the data into parts
	GTEMPBUF(GData*, pSets, nFolds);
	int nSize = pData->size() / nFolds + nFolds;
	int n, i, j, nLearner;
	for(n = 0; n < nFolds; n++)
		pSets[n] = new GData(nSize);
	int nRowCount = pData->size();
	double* pRow;
	for(n = 0; n < nRowCount; n++)
	{
		pRow = pData->row(n);
		pSets[n % nFolds]->AddVector(pRow);
	}

	// Do the training and testing
	double d;
	double dScore = 0;
	GData trainingSet(pData->size());
	for(n = 0; n < nFolds; n++)
	{
		// Train with all of the sub-sets except one
		{
			GReleaseDataHolder hReleaseData(&trainingSet);
			for(i = 0; i < nFolds; i++)
			{
				if(i == n)
					continue;
				int nCount = pSets[i]->size();
				for(j = 0; j < nCount; j++)
				{
					pRow = pSets[i]->row(j);
					trainingSet.AddVector(pRow);
				}
			}


			initialize the accumulator
			for(nLearner = 0; nLearner < nLearner->size(); nLearner++)
			{
				pLearner = m_models[i];
				pLearner->train(&pTrainer);

				eval each row in pSets[n], and accumulate the results

				// Free the model
				pLearner->Clear();
			}
			pAverageSet = generate a set of average results

			// Measure accuracy
			if(bRegression)
				d = MeasureMeanSquaredError(pAverageSet);
			else
				d = MeasurePredictiveAccuracy(pAverageSet);
			dScore += d;
		}
	}
	dScore /= nFolds;

	// Clean up
	for(n = 0; n < nFolds; n++)
	{
		pSets[n]->releaseAllRows();
		delete(pSets[n]);
	}

	return dScore;
}
*/

#ifndef NO_TEST_CODE
#include "GDecisionTree.h"
// static
void GBag::test()
{
	GRand prng(0);
	GBag bag(&prng);
	for(int i = 0; i < 64; i++)
	{
		GDecisionTree* pTree = new GDecisionTree(&prng);
		pTree->useRandomDivisions();
		bag.addLearner(pTree);
	}
	bag.basicTest(0.74, &prng, 0.01);
}
#endif

// -------------------------------------------------------------------------

GBucket::GBucket(GRand* pRand)
: GSupervisedLearner(), m_labelDims(0)
{
	m_nBestLearner = -1;
	m_pRand = pRand;
}

GBucket::GBucket(GTwtNode* pNode, GRand* pRand, GLearnerLoader* pLoader)
: GSupervisedLearner(pNode), m_pRand(pRand)
{
	m_featureDims = (int)pNode->field("featureDims")->asInt();
	m_labelDims = (int)pNode->field("labelDims")->asInt();
	GTwtNode* pModels = pNode->field("models");
	size_t modelCount = pModels->itemCount();
	for(size_t i = 0; i < modelCount; i++)
		m_models.push_back(pLoader->loadModeler(pModels->item(i), pRand));
	m_nBestLearner = (int)pNode->field("best")->asInt();
}

GBucket::~GBucket()
{
	for(vector<GSupervisedLearner*>::iterator it = m_models.begin(); it != m_models.end(); it++)
		delete(*it);
}

// virtual
GTwtNode* GBucket::toTwt(GTwtDoc* pDoc)
{
	GTwtNode* pNode = baseTwtNode(pDoc, "GBucket");
	pNode->addField(pDoc, "featureDims", pDoc->newInt(m_featureDims));
	pNode->addField(pDoc, "labelDims", pDoc->newInt(m_labelDims));
	GTwtNode* pModels = pNode->addField(pDoc, "models", pDoc->newList(m_models.size()));
	for(size_t i = 0; i < m_models.size(); i++)
		pModels->setItem((int)i, m_models[i]->toTwt(pDoc));
	pNode->addField(pDoc, "best", pDoc->newInt(m_nBestLearner));
	return pNode;
}

void GBucket::clear()
{
	for(vector<GSupervisedLearner*>::iterator it = m_models.begin(); it != m_models.end(); it++)
		(*it)->clear();
}

void GBucket::flush()
{
	for(vector<GSupervisedLearner*>::iterator it = m_models.begin(); it != m_models.end(); it++)
		delete(*it);
	m_models.clear();
}

void GBucket::addLearner(GSupervisedLearner* pLearner)
{
	m_models.push_back(pLearner);
}

// virtual
int GBucket::featureDims()
{
	if(m_labelDims < 1)
		ThrowError("not yet trained");
	return m_featureDims;
}

// virtual
int GBucket::labelDims()
{
	if(m_labelDims < 1)
		ThrowError("not yet trained");
	return m_labelDims;
}

// virtual
void GBucket::train(GData* pData, int labelDims)
{
	m_labelDims = labelDims;
	m_featureDims = pData->cols() - labelDims;
	int nLearnerCount = (int)m_models.size();
	double dBestError = 1e200;
	GSupervisedLearner* pLearner;
	m_nBestLearner = (int)m_pRand->next(nLearnerCount);
	int i;
	double err;
	for(i = 0; i < nLearnerCount; i++)
	{
		pLearner = m_models[i];
		try
		{
			err = pLearner->heuristicValidate(pData, labelDims, m_pRand);
		}
		catch(std::exception& e)
		{
			onError(e);
			continue;
		}
		if(err < dBestError)
		{
			dBestError = err;
			m_nBestLearner = i;
		}
		pLearner->clear();
	}
	pLearner = m_models[m_nBestLearner];
	pLearner->train(pData, labelDims);
}

GSupervisedLearner* GBucket::releaseBestModeler()
{
	if(m_nBestLearner < 0)
		ThrowError("Not trained yet");
	GSupervisedLearner* pModeler = m_models[m_nBestLearner];
	m_models[m_nBestLearner] = m_models[m_models.size() - 1];
	m_models.pop_back();
	m_nBestLearner = -1;
	return pModeler;
}

// virtual
void GBucket::predictDistribution(const double* pIn, GPrediction* pOut)
{
	if(m_nBestLearner < 0)
		ThrowError("not trained yet");
	m_models[m_nBestLearner]->predictDistribution(pIn, pOut);
}

// virtual
void GBucket::predict(const double* pIn, double* pOut)
{
	if(m_nBestLearner < 0)
		ThrowError("not trained yet");
	m_models[m_nBestLearner]->predict(pIn, pOut);
}

// virtual
void GBucket::onError(std::exception& e)
{
	//cout << e.what() << "\n";
}

#ifndef NO_TEST_CODE
#include "GDecisionTree.h"
// static
void GBucket::test()
{
	GRand prng(0);
	GBucket bucket(&prng);
	bucket.addLearner(new GBaselineLearner(&prng));
	bucket.addLearner(new GDecisionTree(&prng));
	bucket.basicTest(0.66, &prng);
}
#endif
