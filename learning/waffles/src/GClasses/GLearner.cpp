/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#include "GLearner.h"
#include <stdlib.h>
#include <string.h>
#include "GMacros.h"
#include "GVec.h"
#include "GHeap.h"
#include "GTwt.h"
#include "GImage.h"
#include "GNeuralNet.h"
#include "GKNN.h"
#include "GDecisionTree.h"
#include "GNaiveInstance.h"
#include "GLinear.h"
#include "GNaiveBayes.h"
#include "GEnsemble.h"
#include "GPolynomial.h"
#include "GTransform.h"
#include "GRand.h"
#include "GPlot.h"
#include "GDistribution.h"

using namespace GClasses;

GPrediction::~GPrediction()
{
	delete(m_pDistribution);
}

bool GPrediction::isContinuous()
{
	return m_pDistribution->type() == GUnivariateDistribution::normal;
}

// static
void GPrediction::predictionArrayToVector(int nOutputCount, GPrediction* pOutputs, double* pVector)
{
	int i;
	for(i = 0; i < nOutputCount; i++)
		pVector[i] = pOutputs[i].mode();
}

// static
void GPrediction::vectorToPredictionArray(GRelation* pRelation, int nOutputCount, double* pVector, GPrediction* pOutputs)
{
	int nInputs = pRelation->size() - nOutputCount;
	int i, nValueCount;
	for(i = 0; i < nOutputCount; i++)
	{
		nValueCount = pRelation->valueCount(nInputs + i);
		if(nValueCount == 0)
			pOutputs[i].makeNormal()->setMeanAndVariance(pVector[i], 1);
		else
			pOutputs[i].makeCategorical()->setSpike(nValueCount, (int)pVector[i], 1);
	}
}

double GPrediction::mode()
{
	return m_pDistribution->mode();
}

GCategoricalDistribution* GPrediction::makeCategorical()
{
	if(!m_pDistribution || m_pDistribution->type() != GUnivariateDistribution::categorical)
	{
		delete(m_pDistribution);
		m_pDistribution = new GCategoricalDistribution();
	}
	return (GCategoricalDistribution*)m_pDistribution;
}

GNormalDistribution* GPrediction::makeNormal()
{
	if(!m_pDistribution || m_pDistribution->type() != GUnivariateDistribution::normal)
	{
		delete(m_pDistribution);
		m_pDistribution = new GNormalDistribution();
	}
	return (GNormalDistribution*)m_pDistribution;
}

GCategoricalDistribution* GPrediction::asCategorical()
{
	if(!m_pDistribution || m_pDistribution->type() != GUnivariateDistribution::categorical)
		ThrowError("The current distribution is not a categorical distribution");
	return (GCategoricalDistribution*)m_pDistribution;
}

GNormalDistribution* GPrediction::asNormal()
{
	if(!m_pDistribution || m_pDistribution->type() != GUnivariateDistribution::normal)
		ThrowError("The current distribution is not a normal distribution");
	return (GNormalDistribution*)m_pDistribution;
}

// ---------------------------------------------------------------

GTransducer::GTransducer()
{
}

GTransducer::GTransducer(GTwtNode* pLearner)
{
}

GTransducer::~GTransducer()
{
}

GTwtNode* GTransducer::baseTwtNode(GTwtDoc* pDoc, const char* szClassName)
{
	GTwtNode* pNode = pDoc->newObj();
	pNode->addField(pDoc, "class", pDoc->newString(szClassName));
	return pNode;
}

class GTransducerTrainAndTestCleanUpper
{
protected:
	GData* m_pData;
	size_t m_nTestSize;

public:
	GTransducerTrainAndTestCleanUpper(GData* pData, size_t nTestSize)
	: m_pData(pData), m_nTestSize(nTestSize)
	{
	}

	~GTransducerTrainAndTestCleanUpper()
	{
		while(m_pData->rows() > m_nTestSize)
			m_pData->releaseRow(m_pData->rows() - 1);
	}
};

// virtual
void GTransducer::trainAndTest(GData* pTrainingSet, GData* pTestSet, int labelDims, double* pOutResults)
{
	// Transduce
	int featureDims = pTrainingSet->cols() - labelDims;
	GData dataPredicted(pTrainingSet->relation());
	dataPredicted.newRows(pTestSet->rows());
	dataPredicted.copyColumns(0, pTestSet, 0, featureDims);
	transduce(pTrainingSet, &dataPredicted, labelDims);

	// Evaluate the results
	GVec::setAll(pOutResults, 0.0, labelDims);
	double d, w;
	for(size_t i = 0; i < pTestSet->rows(); i++)
	{
		w = 1.0 / (i + 1);
		double* pMaster = pTestSet->row(i);
		double* pPat = dataPredicted.row(i);
		for(int j = 0; j < labelDims; j++)
		{
			int nIndex = featureDims + j;
			if(pTrainingSet->relation()->valueCount(nIndex) == 0)
			{
				d = pMaster[featureDims + j] - pPat[featureDims + j];
				d *= d;
			}
			else
			{
				if((int)pPat[featureDims + j] == (int)pMaster[featureDims + j])
					d = 1;
				else
					d = 0;
			}
			pOutResults[j] *= (1.0 - w);
			pOutResults[j] += w * d;
		}
	}
}

class GMultipleDataSetsHolder
{
protected:
	GData** m_pSets;
	int m_count;

public:
	GMultipleDataSetsHolder(GData** pSets, int count)
	{
		m_pSets = pSets;
		m_count = count;
	}

	~GMultipleDataSetsHolder()
	{
		int n;
		for(n = 0; n < m_count; n++)
		{
			m_pSets[n]->releaseAllRows();
			delete(m_pSets[n]);
		}
	}
};

double GTransducer::heuristicValidate(GData* pData, int labelDims, GRand* pRand)
{
	GData setA(pData->relation());
	GReleaseDataHolder hA(&setA);
	GData setB(pData->relation());
	GReleaseDataHolder hB(&setB);
	for(size_t i = 0; i < pData->rows(); i++)
	{
		if(pRand->next(2) == 0)
			setA.takeRow(pData->row(i));
		else
			setB.takeRow(pData->row(i));
	}
	GTEMPBUF(double, pResults1, 2 * labelDims);
	double* pResults2 = pResults1 + labelDims;
	trainAndTest(&setA, &setB, labelDims, pResults1);
	trainAndTest(&setB, &setA, labelDims, pResults2);
	int featureDims = pData->cols() - labelDims;
	double err = 0;
	for(int i = 0; i < labelDims; i++)
	{
		if(pData->relation()->valueCount(featureDims + i) == 0)
		{
			err += pResults1[i];
			err += pResults2[i];
		}
		else
		{
			double d = 1.000001 - pResults1[i];
			err += (d * d);
			d = 1.000001 - pResults2[i];
			err += (d * d);
		}
	}
	return err;
}

GData* GTransducer::crossValidate(GData* pData, int nFolds, int labelDims, RepValidateCallback pCB, int nRep, void* pThis)
{
	// Make a place to store the results
	GData* pResults = new GData(labelDims);
	Holder<GData> hResults(pResults);

	// Split the data into parts
	GTEMPBUF(GData*, pSets, nFolds);
	size_t nSize = pData->rows() / nFolds + nFolds;
	for(int n = 0; n < nFolds; n++)
	{
		pSets[n] = new GData(pData->relation(), pData->heap());
		pSets[n]->reserve(nSize);
	}
	GMultipleDataSetsHolder hSets(pSets, nFolds);
	size_t nRowCount = pData->rows();
	double* pRow;
	for(size_t n = 0; n < nRowCount; n++)
	{
		pRow = pData->row(n);
		pSets[n % nFolds]->takeRow(pRow);
	}

	// Do the training and testing
	GData trainingSet(pData->relation(), pData->heap());
	trainingSet.reserve(pData->rows());
	for(int n = 0; n < nFolds; n++)
	{
		double* pFoldResults = pResults->newRow();
		{
			GReleaseDataHolder hReleaseData(&trainingSet);
			for(int i = 0; i < nFolds; i++)
			{
				if(i == n)
					continue;
				size_t nCount = pSets[i]->rows();
				for(size_t j = 0; j < nCount; j++)
				{
					pRow = pSets[i]->row(j);
					trainingSet.takeRow(pRow);
				}
			}

			trainAndTest(&trainingSet, pSets[n], labelDims, pFoldResults);
		}

		// Call the callback
		if(pCB)
			pCB(pThis, nRep, n, labelDims, pFoldResults);
	}
	return hResults.release();
}

GData* GTransducer::repValidate(GData* pData, int nReps, int nFolds, int labelDims, GRand* pRand, RepValidateCallback pCB, void* pThis)
{
	GData* pResults = new GData(labelDims);
	Holder<GData> hResults(pResults);
	for(int i = 0; i < nReps; i++)
	{
		pData->shuffle(pRand);
		GData* pRepResults = crossValidate(pData, nFolds, labelDims, pCB, i, pThis);
		pResults->mergeVert(pRepResults);
		delete(pRepResults);
	}
	return hResults.release();
}

// ---------------------------------------------------------------

GSupervisedLearner::GSupervisedLearner()
: GTransducer()
{
}

GSupervisedLearner::GSupervisedLearner(GTwtNode* pLearner)
: GTransducer(pLearner)
{
}

GSupervisedLearner::~GSupervisedLearner()
{
}

void GSupervisedLearner::accuracy(GData* pData, double* pOutResults)
{
	int labels = labelDims();
	GTEMPBUF(double, out, labels);
	int nIndex;
	double* pIn;
	size_t nRowCount = pData->rows();
	GVec::setAll(pOutResults, 0.0, labels);
	int featureDims = pData->cols() - labels;
	double d, w;
	for(size_t n = 0; n < nRowCount; n++)
	{
		pIn = pData->row(n);
		predict(pIn, out);

		// Check the answer
		w = 1.0 / (n + 1);
		for(int i = 0; i < labels; i++)
		{
			nIndex = featureDims + i;
			if(pData->relation()->valueCount(nIndex) == 0)
			{
				// Squared error
				d = pIn[nIndex] - out[i];
				d *= d;
			}
			else
			{
				// Predictive accuracy
				if((int)pIn[nIndex] == (int)out[i])
					d = 1.0;
				else
					d = 0.0;
			}
			pOutResults[i] *= (1.0 - w);
			pOutResults[i] += w * d;
		}
	}
}

// virtual
void GSupervisedLearner::transduce(GData* pDataLabeled, GData* pDataUnlabeled, int labelDims)
{
	// Train
	train(pDataLabeled, labelDims);

	// Evaluate
	int featureDims = pDataLabeled->cols() - labelDims;
	for(size_t n = 0; n < pDataUnlabeled->rows(); n++)
	{
		double* pPat = pDataUnlabeled->row(n);
		predict(pPat, pPat + featureDims);
	}
}

// virtual
void GSupervisedLearner::trainAndTest(GData* pTrainingSet, GData* pTestSet, int labelDims, double* pOutResults)
{
	train(pTrainingSet, labelDims);
	accuracy(pTestSet, pOutResults);
}

int GSupervisedLearner::precisionRecallContinuous(GPrediction* pOutput, double* pFunc, GData* pTrain, GData* pTest, int labelDims, int label)
{
	// Predict the variance for each pattern
	train(pTrain, labelDims);
	int featureDims = pTrain->cols() - labelDims;
	GData stats(3);
	stats.newRows(pTest->rows());
	for(size_t i = 0; i < pTest->rows(); i++)
	{
		double* pTestVec = pTest->row(i);
		predictDistribution(pTestVec, pOutput);
		double* pResultsVec = stats.row(i);
		pResultsVec[0] = pTestVec[featureDims + label]; // actual label
		GNormalDistribution* pDist = pOutput[label].asNormal();
		pResultsVec[1] = pDist->mean();
		pResultsVec[2] = pDist->variance();
	}

	// Make the precision/recall data
	stats.sort(2); // biggest variance last
	double sse = 0.0;
	for(size_t i = 0; i < stats.rows(); i++)
	{
		double* pVecIn = stats.row(i);
		double d = pVecIn[0] - pVecIn[1];
		sse += (d * d);
		pFunc[i] = sqrt(sse / (i + 1));
	}
	return stats.rows();
}

int GSupervisedLearner::precisionRecallNominal(GPrediction* pOutput, double* pFunc, GData* pTrain, GData* pTest, int labelDims, int label, int value)
{
	// Predict the likelihood that each pattern is relevant
	train(pTrain, labelDims);
	int featureDims = pTrain->cols() - labelDims;
	GData stats(2);
	stats.newRows(pTest->rows());
	int nActualRelevant = 0;
	for(size_t i = 0; i < pTest->rows(); i++)
	{
		double* pTestVec = pTest->row(i);
		predictDistribution(pTestVec, pOutput);
		double* pStatsVec = stats.row(i);
		pStatsVec[0] = pTestVec[featureDims + label]; // actual label
		if((int)pStatsVec[0] == label)
			nActualRelevant++;
		GCategoricalDistribution* pDist = pOutput[label].asCategorical();
		pStatsVec[1] = pDist->likelihood(label); // predicted confidence that it is relevant
	}

	// Make the precision/recall data
	stats.sort(1); // most confident last
	int nFoundRelevant = 0;
	int nFoundTotal = 0;
	for(size_t i = stats.rows() - 1; i < stats.rows(); i--)
	{
		double* pVecIn = stats.row(i);
		nFoundTotal++;
		if((int)pVecIn[0] == label) // if actually relevant
		{
			nFoundRelevant++;
			if(nFoundTotal <= 1)
				pFunc[nFoundRelevant - 1] = 1.0;
			else
				pFunc[nFoundRelevant - 1] = (double)(nFoundRelevant - 1) / (nFoundTotal - 1);
		}
	}
	GAssert(nFoundRelevant == nActualRelevant);
	return nActualRelevant;
}

void GSupervisedLearner::precisionRecall(double* pOutPrecision, int nPrecisionSize, GData* pData, int labelDims, int nOutput, int nReps, GRand* pRand)
{
	int featureDims = pData->cols() - labelDims;
	int nFuncs = MAX(1, pData->relation()->valueCount(featureDims + nOutput));
	GVec::setAll(pOutPrecision, 0.0, nFuncs * nPrecisionSize);
	double* pFunc = new double[pData->rows()];
	ArrayHolder<double> hFunc(pFunc);
#ifdef WIN32
	GPrediction* out = new GPrediction[labelDims];
	ArrayHolder<GPrediction> hOut(out);
#else
	GPrediction out[labelDims];
#endif
	GData dataOtherHalf(pData->relation(), pData->heap());
	int valueCount = pData->relation()->valueCount(featureDims + nOutput);
	for(int nRep = 0; nRep < nReps; nRep++)
	{
		pData->shuffle(pRand);
		pData->splitBySize(&dataOtherHalf, pData->rows() / 2);

		if(valueCount == 0)
		{
			int relevant = precisionRecallContinuous(out, pFunc, pData, &dataOtherHalf, labelDims, nOutput);
			GVec::addInterpolatedFunction(pOutPrecision, nPrecisionSize, pFunc, relevant);
			relevant = precisionRecallContinuous(out, pFunc, &dataOtherHalf, pData, labelDims, nOutput);
			GVec::addInterpolatedFunction(pOutPrecision, nPrecisionSize, pFunc, relevant);
		}
		else
		{
			for(int i = 0; i < valueCount; i++)
			{
				int relevant = precisionRecallNominal(out, pFunc, pData, &dataOtherHalf, labelDims, nOutput, i);
				GVec::addInterpolatedFunction(pOutPrecision + nPrecisionSize * i, nPrecisionSize, pFunc, relevant);
				relevant = precisionRecallNominal(out, pFunc, &dataOtherHalf, pData, labelDims, nOutput, i);
				GVec::addInterpolatedFunction(pOutPrecision + nPrecisionSize * i, nPrecisionSize, pFunc, relevant);
			}
		}
		pData->mergeVert(&dataOtherHalf);
	}
	GVec::multiply(pOutPrecision, 1.0 / (2 * nReps), nFuncs * nPrecisionSize);
}

#ifndef NO_TEST_CODE
void GSupervisedLearner::basicTest(double minAccuracy, GRand* pRand, double deviation)
{
	// Make the training and test data
	GMixedRelation* pMixedRel = new GMixedRelation();
	pMixedRel->addAttr(0);
	pMixedRel->addAttr(0);
	pMixedRel->addAttr(3);
	sp_relation pRel = pMixedRel;
	GData dataTrain(pRel);
	for(size_t i = 0; i < 2000; i++)
	{
		int c = (int)pRand->next(3);
		double* pRow = dataTrain.newRow();
		pRow[0] = pRand->normal() + (c == 1 ? 2.0 : 0.0);
		pRow[1] = pRand->normal() + (c == 2 ? 2.0 : 0.0);
		pRow[2] = (double)c;
	}
	GData dataTest(pRel);
	dataTrain.splitBySize(&dataTest, dataTrain.rows() / 2);

	// Train the model
	train(&dataTrain, 1);
	dataTrain.flush(); // free some memory, just because we can

	// Test the accuracy. A good algorithm will score about 0.7.
	double resultsBefore;
	accuracy(&dataTest, &resultsBefore);
	if(resultsBefore < minAccuracy)
		ThrowError("accuracy has regressed");
	if(resultsBefore > 0.9)
		ThrowError("impossible accuracy");

	// Roundtrip the model through serialization
	GTwtDoc doc;
	doc.setRoot(toTwt(&doc));
	clear(); // free some memory, just because we can
	GLearnerLoader ll;
	GSupervisedLearner* pModel = ll.loadModeler(doc.root(), pRand);
	Holder<GSupervisedLearner> hModel(pModel);

	// Test the accuracy again
	double resultsAfter;
	pModel->accuracy(&dataTest, &resultsAfter);
	if(ABS(resultsAfter - resultsBefore) > deviation)
		ThrowError("serialization shouldn't influence accuracy this much");
}
#endif

// ---------------------------------------------------------------

// virtual
GIncrementalTransform* GLearnerLoader::loadIncrementalTransform(GTwtNode* pNode, GRand* pRand)
{
	const char* szClass = pNode->field("class")->asString();
	if(szClass[0] == 'G')
	{
		if(szClass[1] < 'P')
		{
			if(strcmp(szClass, "GAttributeSelector") == 0)
				return new GAttributeSelector(pNode, pRand);
			else if(strcmp(szClass, "GNoiseGenerator") == 0)
				return new GNoiseGenerator(pNode, pRand);
		}
		else
		{
			if(strcmp(szClass, "GPairProduct") == 0)
				return new GPairProduct(pNode);
			else if(strcmp(szClass, "GPCA") == 0)
				return new GPCA(pNode, pRand);
		}
	}
	return loadTwoWayIncrementalTransform(pNode, pRand);
}

// virtual
GTwoWayIncrementalTransform* GLearnerLoader::loadTwoWayIncrementalTransform(GTwtNode* pNode, GRand* pRand)
{
	const char* szClass = pNode->field("class")->asString();
	if(szClass[0] == 'G')
	{
		if(strcmp(szClass, "GNominalToCat") == 0)
			return new GNominalToCat(pNode);
		else if(strcmp(szClass, "GDiscretize") == 0)
			return new GDiscretize(pNode);
		else if(strcmp(szClass, "GNormalize") == 0)
			return new GNormalize(pNode);
	}
	if(m_throwIfClassNotFound)
		ThrowError("Unrecognized class: ", szClass);
	return NULL;
}

// virtual
GSupervisedLearner* GLearnerLoader::loadModeler(GTwtNode* pNode, GRand* pRand)
{
	const char* szClass = pNode->field("class")->asString();
	if(szClass[0] == 'G')
	{
		if(szClass[1] < 'J')
		{
			if(szClass[1] < 'C')
			{
				if(strcmp(szClass, "GBag") == 0)
					return new GBag(pNode, pRand, this);
				else if(strcmp(szClass, "GBaselineLearner") == 0)
					return new GBaselineLearner(pNode, pRand);
				else if(strcmp(szClass, "GBucket") == 0)
					return new GBucket(pNode, pRand, this);
			}
			else
			{
				if(strcmp(szClass, "GDecisionTree") == 0)
					return new GDecisionTree(pNode, pRand);
				else if(strcmp(szClass, "GFilter") == 0)
					return new GFilter(pNode, pRand, this);
				else if(strcmp(szClass, "GIdentityFunction") == 0)
					return new GIdentityFunction(pNode);
			}
		}
		else
		{
			if(szClass[1] < 'N')
			{
				if(strcmp(szClass, "GLinearRegressor") == 0)
					return new GLinearRegressor(pNode, pRand);
				else if(strcmp(szClass, "GMeanMarginsTree") == 0)
					return new GMeanMarginsTree(pNode, pRand);
			}
			else
			{
				if(strcmp(szClass, "GNaiveMLE") == 0)
					return new GNaiveMLE(pNode);
				else if(strcmp(szClass, "GPolynomial") == 0)
					return new GPolynomial(pNode);
			}
		}
	}
	return loadIncrementalLearner(pNode, pRand);
}

// virtual
GIncrementalLearner* GLearnerLoader::loadIncrementalLearner(GTwtNode* pNode, GRand* pRand)
{
	const char* szClass = pNode->field("class")->asString();
	if(szClass[0] == 'G')
	{
		if(strcmp(szClass, "GKNN") == 0)
			return new GKNN(pNode, pRand);
		else if(strcmp(szClass, "GNaiveBayes") == 0)
			return new GNaiveBayes(pNode, pRand);
		else if(strcmp(szClass, "GNaiveInstance") == 0)
			return new GNaiveInstance(pNode);
		else if(strcmp(szClass, "GNeuralNet") == 0)
			return new GNeuralNet(pNode, pRand);
	}
	if(m_throwIfClassNotFound)
		ThrowError("Unrecognized class: ", szClass);
	return NULL;
}

// ---------------------------------------------------------------

GBaselineLearner::GBaselineLearner(GRand* pRand)
: GSupervisedLearner(), m_featureDims(0), m_labelDims(0), m_pRand(pRand)
{
	m_pLabels = NULL;
}

GBaselineLearner::GBaselineLearner(GTwtNode* pNode, GRand* pRand)
: GSupervisedLearner(pNode), m_pRand(pRand)
{
	GTwtNode* pOutputs = pNode->field("labels");
	int nRealSpaceDims = (int)pOutputs->itemCount();
	m_pLabels = new double[nRealSpaceDims];
	GVec::fromTwt(m_pLabels, nRealSpaceDims, pOutputs);
	m_labelDims = (int)pNode->field("labelDims")->asInt();
	m_pRelation = GRelation::fromTwt(pNode->field("relation"));
}

// virtual
GBaselineLearner::~GBaselineLearner()
{
	clear();
}

// virtual
void GBaselineLearner::clear()
{
	delete[] m_pLabels;
	m_pLabels = NULL;
	m_labelDims = 0;
}

// virtual
GTwtNode* GBaselineLearner::toTwt(GTwtDoc* pDoc)
{
	GTwtNode* pNode = baseTwtNode(pDoc, "GBaselineLearner");
	int nRealSpaceDims = m_pRelation->countRealSpaceDims(m_featureDims, m_labelDims);
	pNode->addField(pDoc, "labels", GVec::toTwt(pDoc, m_pLabels, nRealSpaceDims));
	pNode->addField(pDoc, "relation", m_pRelation->toTwt(pDoc));
	pNode->addField(pDoc, "featureDims", pDoc->newInt(m_featureDims));
	pNode->addField(pDoc, "labelDims", pDoc->newInt(m_labelDims));
	return pNode;
}

// virtual
int GBaselineLearner::featureDims()
{
	if(m_labelDims < 1)
		ThrowError("not yet trained");
	return m_featureDims;
}

// virtual
int GBaselineLearner::labelDims()
{
	if(m_labelDims < 1)
		ThrowError("not yet trained");
	return m_labelDims;
}

// virtual
void GBaselineLearner::train(GData* pData, int labelDims)
{
	clear();
	m_labelDims = labelDims;
	m_pRelation = pData->relation();
	m_featureDims = m_pRelation->size() - m_labelDims;
	int nRealSpaceDims = m_pRelation->countRealSpaceDims(m_featureDims, m_labelDims);
	GTEMPBUF(double, pTmp, nRealSpaceDims);
	delete[] m_pLabels;
	m_pLabels = new double[nRealSpaceDims];
	GVec::setAll(m_pLabels, 0.0, nRealSpaceDims);
	for(size_t i = 0; i < pData->rows(); i++)
	{
		m_pRelation->toRealSpace(pData->row(i) + m_featureDims, pTmp, m_featureDims, m_labelDims);
		GVec::add(m_pLabels, pTmp, nRealSpaceDims);
	}
	GVec::multiply(m_pLabels, 1.0 / pData->rows(), nRealSpaceDims);
}

// virtual
void GBaselineLearner::predictDistribution(const double* pIn, GPrediction* pOut)
{
	m_pRelation->fromRealSpace(m_pLabels, pOut, m_pRelation->size() - m_labelDims, m_labelDims);
}

// virtual
void GBaselineLearner::predict(const double* pIn, double* pOut)
{
	m_pRelation->fromRealSpace(m_pLabels, pOut, m_pRelation->size() - m_labelDims, m_labelDims, m_pRand);
}

#ifndef NO_TEST_CODE
// static
void GBaselineLearner::test()
{
	GRand prng(0);
	GBaselineLearner bl(&prng);
	bl.basicTest(0.31, &prng);
}
#endif

// ---------------------------------------------------------------

GIdentityFunction::GIdentityFunction()
: GSupervisedLearner(), m_labelDims(0), m_featureDims(0)
{
}

GIdentityFunction::GIdentityFunction(GTwtNode* pNode)
: GSupervisedLearner(pNode)
{
	m_labelDims = (int)pNode->field("labels")->asInt();
	m_featureDims = (int)pNode->field("features")->asInt();
}

// virtual
GIdentityFunction::~GIdentityFunction()
{
}

// virtual
void GIdentityFunction::clear()
{
	m_labelDims = 0;
	m_featureDims = 0;
}

// virtual
GTwtNode* GIdentityFunction::toTwt(GTwtDoc* pDoc)
{
	GTwtNode* pNode = baseTwtNode(pDoc, "GIdentityFunction");
	pNode->addField(pDoc, "labels", pDoc->newInt(m_labelDims));
	pNode->addField(pDoc, "features", pDoc->newInt(m_featureDims));
	return pNode;
}

// virtual
int GIdentityFunction::featureDims()
{
	if(m_labelDims < 1)
		ThrowError("not yet trained");
	return m_featureDims;
}

// virtual
int GIdentityFunction::labelDims()
{
	if(m_labelDims < 1)
		ThrowError("not yet trained");
	return m_labelDims;
}

// virtual
void GIdentityFunction::train(GData* pData, int labelDims)
{
	m_labelDims = labelDims;
	m_featureDims = pData->cols() - labelDims;
	if(m_featureDims < 0)
		ThrowError("cannot have more labelDims than columns");
}

// virtual
void GIdentityFunction::predictDistribution(const double* pIn, GPrediction* pOut)
{
	ThrowError("Sorry, not implemented yet");
}

// virtual
void GIdentityFunction::predict(const double* pIn, double* pOut)
{
	if(m_labelDims <= m_featureDims)
		GVec::copy(pOut, pIn, m_labelDims);
	else
	{
		GVec::copy(pOut, pIn, m_featureDims);
		GVec::setAll(pOut + m_featureDims, 0.0, m_labelDims - m_featureDims);
	}
}

