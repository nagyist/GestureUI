/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#include "GNaiveBayes.h"
#include "GMacros.h"
#include <math.h>
#include <stdlib.h>
#include "GVec.h"
#include "GTwt.h"
#include "GDistribution.h"
#include "GRand.h"
#include "GTransform.h"
#include "GSparseMatrix.h"

namespace GClasses {

struct GNaiveBayesInputAttr
{
	size_t m_nValues;
	size_t* m_pValueCounts;

	GNaiveBayesInputAttr(size_t nValues)
	{
		m_nValues = nValues;
		m_pValueCounts = new size_t[m_nValues];
		memset(m_pValueCounts, '\0', sizeof(size_t) * m_nValues);
	}

	GNaiveBayesInputAttr(GTwtNode* pNode)
	{
		m_nValues = pNode->itemCount();
		m_pValueCounts = new size_t[m_nValues];
		for(size_t i = 0; i < m_nValues; i++)
			m_pValueCounts[i] = (size_t)pNode->item(i)->asInt();
	}

	~GNaiveBayesInputAttr()
	{
		delete[] m_pValueCounts;
	}

	GTwtNode* toTwt(GTwtDoc* pDoc)
	{
		GTwtNode* pNode = pDoc->newList(m_nValues);
		for(size_t i = 0; i < m_nValues; i++)
			pNode->setItem(i, pDoc->newInt(m_pValueCounts[i]));
		return pNode;
	}

	void AddTrainingSample(const int inputValue)
	{
		if(inputValue >= 0 && (size_t)inputValue < m_nValues)
			m_pValueCounts[inputValue]++;
		else
			GAssert(inputValue == UNKNOWN_DISCRETE_VALUE);
	}

	int eval(const int inputValue)
	{
		if(inputValue >= 0 && (size_t)inputValue < m_nValues)
			return m_pValueCounts[inputValue];
		else
			return 0;
	}
};

// --------------------------------------------------------------------

struct GNaiveBayesOutputValue
{
	size_t m_nCount;
	size_t m_featureDims;
	struct GNaiveBayesInputAttr** m_pInputs;

	GNaiveBayesOutputValue(GRelation* pRelation, size_t nInputCount)
	{
		m_nCount = 0;
		m_featureDims = nInputCount;
		m_pInputs = new struct GNaiveBayesInputAttr*[nInputCount];
		for(size_t n = 0; n < m_featureDims; n++)
			m_pInputs[n] = new struct GNaiveBayesInputAttr(pRelation->valueCount(n));
	}

	GNaiveBayesOutputValue(GTwtNode* pNode, size_t nInputCount)
	{
		if(pNode->itemCount() != nInputCount + 1)
			ThrowError("Unexpected number of inputs");
		m_nCount = (size_t)pNode->item(0)->asInt();
		m_featureDims = nInputCount;
		m_pInputs = new struct GNaiveBayesInputAttr*[m_featureDims];
		for(size_t n = 0; n < m_featureDims; n++)
			m_pInputs[n] = new struct GNaiveBayesInputAttr(pNode->item(n + 1));
	}

	~GNaiveBayesOutputValue()
	{
		for(size_t n = 0; n < m_featureDims; n++)
			delete(m_pInputs[n]);
		delete[] m_pInputs;
	}

	GTwtNode* toTwt(GTwtDoc* pDoc)
	{
		GTwtNode* pNode = pDoc->newList(m_featureDims + 1);
		pNode->setItem(0, pDoc->newInt(m_nCount));
		for(size_t i = 0; i < m_featureDims; i++)
			pNode->setItem(i + 1, m_pInputs[i]->toTwt(pDoc));
		return pNode;
	}

	void AddTrainingSample(const double* pIn)
	{
		for(size_t n = 0; n < m_featureDims; n++)
			m_pInputs[n]->AddTrainingSample((int)pIn[n]);
		m_nCount++;
	}

	double eval(const double* pInputVector, double equivalentSampleSize)
	{
		// The prior output probability
		double dLogProb = log((double)m_nCount);

		// The probability of inputs given this output
		for(size_t n = 0; n < m_featureDims; n++)
		{
			dLogProb += log(MAX(1e-300,
					(
						(double)m_pInputs[n]->eval((int)pInputVector[n]) + 
						(equivalentSampleSize / m_pInputs[n]->m_nValues)
					) / 
					(equivalentSampleSize + m_nCount)
				));
		}
		return dLogProb;
	}
};

// --------------------------------------------------------------------

struct GNaiveBayesOutputAttr
{
	size_t m_nValueCount;
	struct GNaiveBayesOutputValue** m_pValues;

	GNaiveBayesOutputAttr(GRelation* pRelation, size_t nInputCount, size_t nValueCount)
	{
		m_nValueCount = nValueCount;
		m_pValues = new struct GNaiveBayesOutputValue*[m_nValueCount];
		for(size_t n = 0; n < m_nValueCount; n++)
			m_pValues[n] = new struct GNaiveBayesOutputValue(pRelation, nInputCount);
	}

	GNaiveBayesOutputAttr(GTwtNode* pNode, size_t nInputCount, size_t nValueCount)
	{
		if(pNode->itemCount() != nValueCount)
			ThrowError("Unexpected number of values");
		m_nValueCount = nValueCount;
		m_pValues = new struct GNaiveBayesOutputValue*[m_nValueCount];
		for(size_t n = 0; n < m_nValueCount; n++)
			m_pValues[n] = new struct GNaiveBayesOutputValue(pNode->item(n), nInputCount);
	}

	~GNaiveBayesOutputAttr()
	{
		for(size_t n = 0; n < m_nValueCount; n++)
			delete(m_pValues[n]);
		delete[] m_pValues;
	}

	GTwtNode* toTwt(GTwtDoc* pDoc)
	{
		GTwtNode* pNode = pDoc->newList(m_nValueCount);
		for(size_t i = 0; i < m_nValueCount; i++)
			pNode->setItem(i, m_pValues[i]->toTwt(pDoc));
		return pNode;
	}

	void AddTrainingSample(const double* pIn, int out)
	{
		if(out >= 0 && (size_t)out < m_nValueCount)
			m_pValues[out]->AddTrainingSample(pIn);
	}

	void eval(const double* pIn, GPrediction* pOut, double equivalentSampleSize)
	{
		GCategoricalDistribution* pDist = pOut->makeCategorical();
		double* pValues = pDist->values(m_nValueCount);
		for(size_t n = 0; n < m_nValueCount; n++)
			pValues[n] = m_pValues[n]->eval(pIn, equivalentSampleSize);
		pDist->normalizeFromLogSpace();
	}

	double predict(const double* pIn, double equivalentSampleSize, GRand* pRand)
	{
		GTEMPBUF(double, pValues, m_nValueCount);
		for(size_t n = 0; n < m_nValueCount; n++)
			pValues[n] = m_pValues[n]->eval(pIn, equivalentSampleSize);
		return (double)GVec::indexOfMax(pValues, m_nValueCount, pRand);
	}
};

// --------------------------------------------------------------------

GNaiveBayes::GNaiveBayes(GRand* pRand)
: GIncrementalLearner(), m_labelDims(0), m_pRand(pRand)
{
	m_pOutputs = NULL;
	m_equivalentSampleSize = 0.5;
	m_nSampleCount = 0;
}

GNaiveBayes::GNaiveBayes(GTwtNode* pNode, GRand* pRand)
: GIncrementalLearner(pNode), m_pRand(pRand)
{
	m_pRelation = GRelation::fromTwt(pNode->field("relation"));
	m_labelDims = (size_t)pNode->field("labelDims")->asInt();
	m_nSampleCount = (size_t)pNode->field("sampleCount")->asInt();
	m_equivalentSampleSize = pNode->field("ess")->asDouble();
	GTwtNode* pOutputs = pNode->field("outputs");
	if(pOutputs->itemCount() != m_labelDims)
		ThrowError("Wrong number of outputs");
	m_pOutputs = new struct GNaiveBayesOutputAttr*[m_labelDims];
	int featureDims = m_pRelation->size() - m_labelDims;
	for(size_t i = 0; i < m_labelDims; i++)
		m_pOutputs[i] = new struct GNaiveBayesOutputAttr(pOutputs->item(i), featureDims, m_pRelation->valueCount(featureDims + i));
}

GNaiveBayes::~GNaiveBayes()
{
	clear();
}

// virtual
GTwtNode* GNaiveBayes::toTwt(GTwtDoc* pDoc)
{
	GTwtNode* pNode = baseTwtNode(pDoc, "GNaiveBayes");
	pNode->addField(pDoc, "relation", m_pRelation->toTwt(pDoc));
	pNode->addField(pDoc, "labelDims", pDoc->newInt(m_labelDims));
	pNode->addField(pDoc, "sampleCount", pDoc->newInt(m_nSampleCount));
	pNode->addField(pDoc, "ess", pDoc->newDouble(m_equivalentSampleSize));
	GTwtNode* pOutputs = pNode->addField(pDoc, "outputs", pDoc->newList(m_labelDims));
	for(size_t i = 0; i < m_labelDims; i++)
		pOutputs->setItem(i, m_pOutputs[i]->toTwt(pDoc));
	return pNode;
}

// virtual
int GNaiveBayes::featureDims()
{
	if(m_labelDims < 1)
		ThrowError("not yet trained");
	return m_pRelation->size() - m_labelDims;
}

// virtual
int GNaiveBayes::labelDims()
{
	if(m_labelDims < 1)
		ThrowError("not yet trained");
	return m_labelDims;
}

// virtual
void GNaiveBayes::clear()
{
	m_nSampleCount = 0;
	if(m_pOutputs)
	{
		for(size_t n = 0; n < m_labelDims; n++)
			delete(m_pOutputs[n]);
		delete[] m_pOutputs;
	}
	m_pOutputs = NULL;
}

// virtual
void GNaiveBayes::enableIncrementalLearning(sp_relation& pRelation, int labelDims, double* pMins, double* pRanges)
{
	clear();
	if(!pRelation->areNominal(0, pRelation->size()))
		ThrowError("GNaiveBayes does not support continuous attributes. You should discretize first to convert real values to nominals.");
	m_labelDims = labelDims;
	m_pRelation = pRelation;
	m_pOutputs = new struct GNaiveBayesOutputAttr*[m_labelDims];
	int featureDims = pRelation->size() - labelDims;
	for(size_t n = 0; n < m_labelDims; n++)
		m_pOutputs[n] = new struct GNaiveBayesOutputAttr(m_pRelation.get(), featureDims, m_pRelation->valueCount(featureDims + n));
}

// virtual
void GNaiveBayes::trainIncremental(const double* pIn, const double* pOut)
{
	for(size_t n = 0; n < m_labelDims; n++)
		m_pOutputs[n]->AddTrainingSample(pIn, (int)pOut[n]);
	m_nSampleCount++;
}

// virtual
void GNaiveBayes::train(GData* pData, int labelDims)
{
	enableIncrementalLearning(pData->relation(), labelDims, NULL, NULL);
	int featureDims = pData->cols() - m_labelDims;
	for(size_t n = 0; n < pData->rows(); n++)
	{
		double* pRow = pData->row(n);
		trainIncremental(pRow, pRow + featureDims);
	}
}

// virtual
void GNaiveBayes::trainSparse(GSparseMatrix* pData, int labelDims)
{
	sp_relation pRel = new GUniformRelation(pData->cols(), 2);
	enableIncrementalLearning(pRel, labelDims, NULL, NULL);
	int featureDims = pData->cols() - m_labelDims;
	double* pFullRow = new double[pData->cols()];
	ArrayHolder<double> hFullRow(pFullRow);
	for(unsigned int n = 0; n < pData->rows(); n++)
	{
		pData->fullRow(pFullRow, n);
		trainIncremental(pFullRow, pFullRow + featureDims);
	}
}

void GNaiveBayes::predictDistribution(const double* pIn, GPrediction* pOut)
{
	if(m_nSampleCount <= 0)
		ThrowError("You must call train before you call eval");
	for(size_t n = 0; n < m_labelDims; n++)
		m_pOutputs[n]->eval(pIn, &pOut[n], m_equivalentSampleSize);
}

void GNaiveBayes::predict(const double* pIn, double* pOut)
{
	if(m_nSampleCount <= 0)
		ThrowError("You must call train before you call eval");
	for(size_t n = 0; n < m_labelDims; n++)
		pOut[n] = m_pOutputs[n]->predict(pIn, m_equivalentSampleSize, m_pRand);
}

#ifndef NO_TEST_CODE
void GNaiveBayes_CheckResults(double yprior, double ycond, double nprior, double ncond, GPrediction* out)
{
	double py = yprior * ycond;
	double pn = nprior * ncond;
	double sum = py + pn;
	py /= sum;
	pn /= sum;
	GCategoricalDistribution* pCat = out->asCategorical();
	double* pVals = pCat->values(2);
	if(ABS(pVals[0] - py) > 1e-8)
		ThrowError("wrong");
	if(ABS(pVals[1] - pn) > 1e-8)
		ThrowError("wrong");
}

void GNaiveBayes_testMath(GRand* pRand)
{
	const char* trainFile =
	"@RELATION test\n"
	"@ATTRIBUTE a {t,f}\n"
	"@ATTRIBUTE b {r,g,b}\n"
	"@ATTRIBUTE c {y,n}\n"
	"@DATA\n"
	"t,r,y\n"
	"f,r,n\n"
	"t,g,y\n"
	"f,g,y\n"
	"f,g,n\n"
	"t,r,n\n"
	"t,r,y\n"
	"t,b,y\n"
	"f,r,y\n"
	"f,g,n\n"
	"f,b,y\n"
	"t,r,n\n";
	GData* pTrain = GData::parseArff(trainFile, strlen(trainFile));
	Holder<GData> hTrain(pTrain);
	GNaiveBayes nb(pRand);
	nb.setEquivalentSampleSize(0.0);
	nb.train(pTrain, 1);
	GPrediction out;
	double pat[2];
	pat[0] = 0; pat[1] = 0;
	nb.predictDistribution(pat, &out);
	GNaiveBayes_CheckResults(7.0/12.0, 4.0/7.0*3.0/7.0, 5.0/12.0, 2.0/5.0*3.0/5.0, &out);
	pat[0] = 0; pat[1] = 1;
	nb.predictDistribution(pat, &out);
	GNaiveBayes_CheckResults(7.0/12.0, 4.0/7.0*2.0/7.0, 5.0/12.0, 2.0/5.0*2.0/5.0, &out);
	pat[0] = 0; pat[1] = 2;
	nb.predictDistribution(pat, &out);
	GNaiveBayes_CheckResults(7.0/12.0, 4.0/7.0*2.0/7.0, 5.0/12.0, 2.0/5.0*0.0/5.0, &out);
	pat[0] = 1; pat[1] = 0;
	nb.predictDistribution(pat, &out);
	GNaiveBayes_CheckResults(7.0/12.0, 3.0/7.0*3.0/7.0, 5.0/12.0, 3.0/5.0*3.0/5.0, &out);
	pat[0] = 1; pat[1] = 1;
	nb.predictDistribution(pat, &out);
	GNaiveBayes_CheckResults(7.0/12.0, 3.0/7.0*2.0/7.0, 5.0/12.0, 3.0/5.0*2.0/5.0, &out);
	pat[0] = 1; pat[1] = 2;
	nb.predictDistribution(pat, &out);
	GNaiveBayes_CheckResults(7.0/12.0, 3.0/7.0*2.0/7.0, 5.0/12.0, 3.0/5.0*0.0/5.0, &out);
}

// static
void GNaiveBayes::test()
{
	GRand prng(0);
	GNaiveBayes_testMath(&prng);
	GNaiveBayes nb(&prng);
	GFilter tl(&nb, false);
	tl.setFeatureTransform(new GDiscretize(), true);
	tl.basicTest(0.76, &prng);
}
#endif // !NO_TEST_CODE

// -----------------------------------------------------------------------------

GNaiveMLE::GNaiveMLE(sp_relation& pRelation)
: GSupervisedLearner(), m_labelDims(0)
{
	m_pPredictions = NULL;
	m_dLimboValue = 1.0;
	m_dEquivalentSampleSize = .001;
}

GNaiveMLE::GNaiveMLE(GTwtNode* pNode)
: GSupervisedLearner(pNode)
{
	m_pRelation = GRelation::fromTwt(pNode->field("relation"));
	m_labelDims = (int)pNode->field("labelDims")->asInt();
	m_dLimboValue = pNode->field("limbo")->asDouble();
	m_dEquivalentSampleSize = pNode->field("ess")->asDouble();
	GTwtNode* pPredictions = pNode->field("predictions");
	m_nValues = (int)pPredictions->itemCount();
	m_pPredictions = new GCategoricalDistribution[m_nValues];
	int i;
	for(i = 0; i < m_nValues; i++)
		m_pPredictions[i].fromTwt(pPredictions->item(i));
}

// virtual
GNaiveMLE::~GNaiveMLE()
{
	delete[] m_pPredictions;
}

GTwtNode* GNaiveMLE::toTwt(GTwtDoc* pDoc)
{
	GTwtNode* pNode = baseTwtNode(pDoc, "GNaiveMLE");
	pNode->addField(pDoc, "relation", m_pRelation->toTwt(pDoc));
	pNode->addField(pDoc, "labelDims", pDoc->newInt(m_labelDims));
	pNode->addField(pDoc, "limbo", pDoc->newDouble(m_dLimboValue));
	pNode->addField(pDoc, "ess", pDoc->newDouble(m_dEquivalentSampleSize));
	GTwtNode* pPredictions = pNode->addField(pDoc, "predictions", pDoc->newList(m_nValues));
	int i;
	for(i = 0; i < m_nValues; i++)
		pPredictions->setItem(i, m_pPredictions[i].toTwt(pDoc));
	return pNode;
}

// virtual
int GNaiveMLE::featureDims()
{
	if(m_labelDims < 1)
		ThrowError("not yet trained");
	return m_pRelation->size() - m_labelDims;
}

// virtual
int GNaiveMLE::labelDims()
{
	if(m_labelDims < 1)
		ThrowError("not yet trained");
	return m_labelDims;
}

// virtual
void GNaiveMLE::clear()
{
	delete[] m_pPredictions;
	m_pPredictions = NULL;
}

// virtual
void GNaiveMLE::train(GData* pData, int labelDims)
{
	clear();
	m_pRelation = pData->relation();
	if(!m_pRelation->areNominal(0, m_pRelation->size()))
		ThrowError("GNaiveMLE does not support continuous attributes. You should discretize first to convert real values to nominals.");
	m_nValues = 0;
	int featureDims = m_pRelation->size() - labelDims;
	for(int i = 0; i < featureDims; i++)
		m_nValues += m_pRelation->valueCount(i);
	m_pPredictions = new GCategoricalDistribution[m_nValues];
	m_labelDims = labelDims;

	// Clear the predictions
	int nOutputValues = pData->relation()->valueCount(featureDims);
	int j, k, valIn, valOut, valCount;
	GCategoricalDistribution* pCat;
	double* pValues;
	int nPrediction = 0;
	for(j = 0; j < featureDims; j++)
	{
		valCount = m_pRelation->valueCount(j);
		for(k = 0; k < valCount; k++)
		{
			pCat = &m_pPredictions[nPrediction++];
			pValues = pCat->values(nOutputValues);
			for(int i = 0; i < nOutputValues; i++)
				pValues[i] = m_dEquivalentSampleSize;
		}
	}

	// Accumulate the predictions
	double* pVec;
	for(size_t i = 0; i < pData->rows(); i++)
	{
		nPrediction = 0;
		pVec = pData->row(i);
		for(j = 0; j < featureDims; j++)
		{
			valCount = m_pRelation->valueCount(j);
			valIn = (int)pVec[j];
			valOut = (int)pVec[featureDims];
			if(valIn >= 0 && valIn < valCount && valOut >= 0 && valOut < nOutputValues)
			{
				pCat = &m_pPredictions[nPrediction + valIn];
				pValues = pCat->values(nOutputValues);
				pValues[valOut]++;
			}
			nPrediction += valCount;
		}
	}

	// Normalize the predictions
	nPrediction = 0;
	for(j = 0; j < featureDims; j++)
	{
		valCount = m_pRelation->valueCount(j);
		for(k = 0; k < valCount; k++)
		{
			pCat = &m_pPredictions[nPrediction++];
			pCat->normalize();
		}
	}
}

// virtual
void GNaiveMLE::predictDistribution(const double* pIn, GPrediction* pOut)
{
	int featureDims = m_pRelation->size() - m_labelDims;
	int nOutputValues = m_pRelation->valueCount(featureDims);
	GCategoricalDistribution* pCatOut = pOut->makeCategorical();
	double* pValuesOut = pCatOut->values(nOutputValues);
	GVec::setAll(pValuesOut, 0.0, nOutputValues);
	GCategoricalDistribution* pCatIn;
	int i, j, val, valCount;
	int nPrediction = 0;
	double x;
	for(i = 0; i < featureDims; i++)
	{
		valCount = m_pRelation->valueCount(i);
		val = (int)pIn[i];
		if(val >= 0 && val < valCount)
		{
			pCatIn = &m_pPredictions[nPrediction + val];
			for(j = 0; j < nOutputValues; j++)
			{
				x = pCatIn->likelihood(j);
				pValuesOut[j] += (1.0 - m_dLimboValue) * x + m_dLimboValue * log(x);
			}
		}
		nPrediction += valCount;
	}
	GAssert(nPrediction == m_nValues);
	pCatOut->normalizeFromLogSpace();
}

// virtual
void GNaiveMLE::predict(const double* pIn, double* pOut)
{
#ifdef WIN32
	GPrediction* out = new GPrediction[m_labelDims];
	ArrayHolder<GPrediction> hOut(out);
#else
	GPrediction out[m_labelDims];
#endif
	predictDistribution(pIn, out);
	GPrediction::predictionArrayToVector(m_labelDims, out, pOut);
}

} // namespace GClasses
