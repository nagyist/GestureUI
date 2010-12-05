/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#include "GTransform.h"
#include "GTwt.h"
#include "GVec.h"
#include "GRand.h"
#include "GManifold.h"
#include "GCluster.h"
#include "GString.h"
#include "GNeuralNet.h"
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <math.h>
#include <string>

namespace GClasses {

using std::string;
using std::vector;

GTransform::GTransform()
{
}

GTransform::GTransform(GTwtNode* pNode)
{
}

// virtual
GTransform::~GTransform()
{
}

// virtual
GTwtNode* GTransform::baseTwtNode(GTwtDoc* pDoc, const char* szClassName)
{
	GTwtNode* pNode = pDoc->newObj();
	pNode->addField(pDoc, "class", pDoc->newString(szClassName));
	return pNode;
}

// ---------------------------------------------------------------

GTransformChainer::GTransformChainer(GTransform* pFirst, GTransform* pSecond)
: GTransform()
{
	m_pFirst = pFirst;
	m_pSecond = pSecond;
}

// virtual
GTransformChainer::~GTransformChainer()
{
	delete(m_pFirst);
	delete(m_pSecond);
}

// virtual
GData* GTransformChainer::doit(GData* pIn)
{
	GData* pDataTmp = m_pFirst->doit(pIn);
	Holder<GData> hDataTmp(pDataTmp);
	return m_pSecond->doit(pDataTmp);
}

// ---------------------------------------------------------------

// virtual
GData* GIncrementalTransform::doit(GData* pIn)
{
	train(pIn);
	return transformBatch(pIn);
}

// virtual
GData* GIncrementalTransform::transformBatch(GData* pIn)
{
	if(!m_pRelationBefore.get())
		ThrowError("neither train nor enableIncrementalTraining has been called");
	size_t nRows = pIn->rows();
	GData* pOut = new GData(after());
	Holder<GData> hOut(pOut);
	pOut->newRows(nRows);
	for(size_t i = 0; i < nRows; i++)
		transform(pIn->row(i), pOut->row(i));
	return hOut.release();
}

// ---------------------------------------------------------------

// virtual
GData* GTwoWayIncrementalTransform::untransformBatch(GData* pIn)
{
	if(!m_pRelationBefore.get())
		ThrowError("neither train nor enableIncrementalTraining has been called");
	size_t nRows = pIn->rows();
	GData* pOut = new GData(before());
	Holder<GData> hOut(pOut);
	for(size_t i = 0; i < nRows; i++)
		untransform(pIn->row(i), pOut->row(i));
	return hOut.release();
}

// ---------------------------------------------------------------

GFilter::GFilter(GSupervisedLearner* pModeler, bool own)
: GIncrementalLearner(), m_beforeLabelDims(0), m_afterLabelDims(0), m_pFeatureTransform(NULL), m_pLabelTransform(NULL), m_pModeler(pModeler), m_ownFeatureTransform(false), m_ownLabelTransform(false), m_ownModeler(own), m_pScratchVector(NULL), m_pScratchVector2(NULL), m_pScratchLabels(NULL)
{
}

GFilter::GFilter(GTwtNode* pNode, GRand* pRand, GLearnerLoader* pLoader)
: GIncrementalLearner(pNode)
{
	m_pModeler = pLoader->loadModeler(pNode->field("modeler"), pRand);
	GTwtNode* pFeature = pNode->fieldIfExists("feature");
	if(pFeature)
		m_pFeatureTransform = pLoader->loadIncrementalTransform(pFeature, pRand);
	else
		m_pFeatureTransform = NULL;
	GTwtNode* pLabel = pNode->fieldIfExists("label");
	if(pLabel)
		m_pLabelTransform = pLoader->loadTwoWayIncrementalTransform(pLabel, pRand);
	else
		m_pLabelTransform = NULL;
	m_ownFeatureTransform = true;
	m_ownLabelTransform = true;
	m_ownModeler = true;
	m_beforeLabelDims = (int)pNode->field("beforeLabels")->asInt();
	m_afterLabelDims = (int)pNode->field("afterLabels")->asInt();
	m_pBeforeRelation = GRelation::fromTwt(pNode->field("beforeRelation"));
	m_pAfterRelation = GRelation::fromTwt(pNode->field("afterRelation"));
	m_pScratchVector = NULL;
	m_pScratchLabels = NULL;
	makeScratchBuffers();
}

// virtual
GFilter::~GFilter()
{
	if(m_ownFeatureTransform)
		delete(m_pFeatureTransform);
	if(m_ownLabelTransform)
		delete(m_pLabelTransform);
	if(m_ownModeler)
		delete(m_pModeler);
	delete[] m_pScratchVector;
	delete[] m_pScratchLabels;
}

// virtual
GTwtNode* GFilter::toTwt(GTwtDoc* pDoc)
{
	GTwtNode* pNode = baseTwtNode(pDoc, "GFilter");
	if(m_pFeatureTransform)
		pNode->addField(pDoc, "feature", m_pFeatureTransform->toTwt(pDoc));
	if(m_pLabelTransform)
		pNode->addField(pDoc, "label", m_pLabelTransform->toTwt(pDoc));
	pNode->addField(pDoc, "modeler", m_pModeler->toTwt(pDoc));
	pNode->addField(pDoc, "beforeLabels", pDoc->newInt(m_beforeLabelDims));
	pNode->addField(pDoc, "afterLabels", pDoc->newInt(m_afterLabelDims));
	pNode->addField(pDoc, "beforeRelation", m_pBeforeRelation->toTwt(pDoc));
	pNode->addField(pDoc, "afterRelation", m_pAfterRelation->toTwt(pDoc));
	return pNode;
}

// virtual
void GFilter::clear()
{
	m_pModeler->clear();
}

void GFilter::setFeatureTransform(GIncrementalTransform* pFeatureTransform, bool takeOwnership)
{
	if(m_ownFeatureTransform)
		delete(m_pFeatureTransform);
	m_pFeatureTransform = pFeatureTransform;
	m_ownFeatureTransform = takeOwnership;
}

void GFilter::setLabelTransform(GTwoWayIncrementalTransform* pLabelTransform, bool takeOwnership)
{
	if(m_ownLabelTransform)
		delete(m_pLabelTransform);
	m_pLabelTransform = pLabelTransform;
	m_ownLabelTransform = takeOwnership;
}

// virtual
int GFilter::featureDims()
{
	if(m_beforeLabelDims < 1)
		ThrowError("not yet trained");
	return m_pBeforeRelation->size() - m_beforeLabelDims;
}

// virtual
int GFilter::labelDims()
{
	if(m_beforeLabelDims < 1)
		ThrowError("not yet trained");
	return m_beforeLabelDims;
}

void GFilter::makeScratchBuffers()
{
	int scratch1 = 0;
	int scratch2 = 0;
	if(m_pFeatureTransform)
		scratch1 = m_pAfterRelation->size() - m_afterLabelDims;
	if(m_pLabelTransform)
	{
		scratch1 = MAX(scratch1, m_beforeLabelDims);
		scratch2 = m_afterLabelDims;
	}
	if(scratch1 + scratch2 > 0)
	{
		delete[] m_pScratchVector;
		m_pScratchVector = new double[scratch1 + scratch2];
		m_pScratchVector2 = m_pScratchVector + scratch1;
		if(scratch2 > 0)
			m_pScratchLabels = new GPrediction[scratch2];
		else
			m_pScratchLabels = NULL;
	}
	else
	{
		m_pScratchVector = NULL;
		m_pScratchVector2 = NULL;
		m_pScratchLabels = NULL;
	}
}

// virtual
void GFilter::predictDistribution(const double* pIn, GPrediction* pOut)
{
	if(m_pFeatureTransform)
	{
		m_pFeatureTransform->transform(pIn, m_pScratchVector);
		pIn = m_pScratchVector;
	}
	if(m_pLabelTransform)
	{
		m_pModeler->predictDistribution(pIn, m_pScratchLabels);
		GPrediction::predictionArrayToVector(m_afterLabelDims, m_pScratchLabels, m_pScratchVector2);
		m_pLabelTransform->untransform(m_pScratchVector2, m_pScratchVector);
		GPrediction::vectorToPredictionArray(m_pBeforeRelation.get(), m_beforeLabelDims, m_pScratchVector, pOut);
	}
	else
		m_pModeler->predictDistribution(pIn, pOut);
}

// virtual
void GFilter::predict(const double* pIn, double* pOut)
{
	if(m_pFeatureTransform)
	{
		m_pFeatureTransform->transform(pIn, m_pScratchVector);
		pIn = m_pScratchVector;
	}
	if(m_pLabelTransform)
	{
		m_pModeler->predict(pIn, m_pScratchVector2);
		m_pLabelTransform->untransform(m_pScratchVector2, pOut);
	}
	else
		m_pModeler->predict(pIn, pOut);
}

// virtual
void GFilter::train(GData* pIn, int labelDims)
{
	m_beforeLabelDims = labelDims;
	m_afterLabelDims = labelDims;
	m_pBeforeRelation = pIn->relation();
	if(!m_pFeatureTransform && !m_pLabelTransform)
	{
		m_pAfterRelation = m_pBeforeRelation;
		m_pModeler->train(pIn, labelDims);
		return;
	}

	// Transform the features
	int featureDims = pIn->cols() - labelDims;
	GData* pFeatures = pIn->attrSubset(0, featureDims);
	Holder<GData> hFeatures(pFeatures);
	GData* pTransformedFeatures = pFeatures;
	Holder<GData> hTransformedFeatures(NULL);
	if(m_pFeatureTransform)
	{
		m_pFeatureTransform->train(pFeatures);
		pTransformedFeatures = m_pFeatureTransform->transformBatch(pFeatures);
		hTransformedFeatures.reset(pTransformedFeatures);
	}

	// Transform the labels
	GData* pLabels = pIn->attrSubset(featureDims, labelDims);
	Holder<GData> hLabels(pLabels);
	GData* pTransformedLabels = pLabels;
	Holder<GData> hTransformedLabels(NULL);
	if(m_pLabelTransform)
	{
		m_pLabelTransform->train(pLabels);
		pTransformedLabels = m_pLabelTransform->transformBatch(pLabels);
		hTransformedLabels.reset(pTransformedLabels);
		m_afterLabelDims = pTransformedLabels->cols();
	}

	// Train the modeler
	GData* pUnified = GData::mergeHoriz(pTransformedFeatures, pTransformedLabels);
	Holder<GData> hUnified(pUnified);
	m_pAfterRelation = pUnified->relation();
	m_pModeler->train(pUnified, m_afterLabelDims);
	makeScratchBuffers();
}

// virtual
void GFilter::enableIncrementalLearning(sp_relation& pRelation, int labelDims, double* pMins, double* pRanges)
{
	m_beforeLabelDims = labelDims;
	m_pBeforeRelation = pRelation;
	if(!m_pModeler->canTrainIncrementally())
		ThrowError("The wrapped model does not support incremental training");
	if(!m_pFeatureTransform && !m_pLabelTransform)
	{
		m_pAfterRelation = m_pBeforeRelation;
		((GIncrementalLearner*)m_pModeler)->enableIncrementalLearning(pRelation, labelDims, pMins, pRanges);
		return;
	}

	// Make the after relation and enable incremental training on the transforms
	GMixedRelation* pRelAfter = new GMixedRelation();
	m_pAfterRelation = pRelAfter;
	int featureDims = pRelation->size() - labelDims;
	if(m_pFeatureTransform)
	{
		sp_relation pRelFeaturesBefore = new GMixedRelation();
		((GMixedRelation*)pRelFeaturesBefore.get())->addAttrs(pRelation.get(), 0, featureDims);
		m_pFeatureTransform->enableIncrementalTraining(pRelFeaturesBefore, pMins, pRanges);
		pRelAfter->addAttrs(m_pFeatureTransform->after().get());
	}
	else
		pRelAfter->addAttrs(pRelation.get(), 0, pRelation->size() - labelDims);
	if(m_pLabelTransform)
	{
		sp_relation pRelLabelsBefore = new GMixedRelation();
		((GMixedRelation*)pRelLabelsBefore.get())->addAttrs(pRelation.get(), featureDims, labelDims);
		m_pLabelTransform->enableIncrementalTraining(pRelLabelsBefore, pMins + featureDims, pRanges + featureDims);
		pRelAfter->addAttrs(m_pLabelTransform->after().get());
	}
	else
		pRelAfter->addAttrs(pRelation.get(), pRelation->size() - labelDims, labelDims);
	if(m_pLabelTransform)
		m_afterLabelDims = m_pLabelTransform->after()->size();
	else
		m_afterLabelDims = labelDims;

	// Allocate scratch buffers
	makeScratchBuffers();

	// Make the after mins and ranges
	int afterDims = pRelAfter->size();
	int afterFeatureDims = afterDims - m_afterLabelDims;
	double* pAfterMins = new double[afterDims * 2];
	ArrayHolder<double> hAfterMins(pAfterMins);
	double* pAfterRanges = pAfterMins + afterDims;
	if(m_pFeatureTransform)
	{
		GVec::copy(pAfterMins, m_pFeatureTransform->afterMins(), afterFeatureDims);
		GVec::copy(pAfterRanges, m_pFeatureTransform->afterRanges(), afterFeatureDims);
	}
	else
	{
		GVec::copy(pAfterMins, pMins, afterFeatureDims);
		GVec::copy(pAfterRanges, pRanges, afterFeatureDims);
	}
	if(m_pLabelTransform)
	{
		GVec::copy(pAfterMins + afterFeatureDims, m_pLabelTransform->afterMins(), m_afterLabelDims);
		GVec::copy(pAfterRanges + afterFeatureDims, m_pLabelTransform->afterRanges(), m_afterLabelDims);
	}
	else
	{
		GVec::copy(pAfterMins + afterFeatureDims, pMins + featureDims, m_afterLabelDims);
		GVec::copy(pAfterRanges + afterFeatureDims, pRanges + featureDims, m_afterLabelDims);
	}

	// Enable incremental learning on the modeler
	((GIncrementalLearner*)m_pModeler)->enableIncrementalLearning(m_pAfterRelation, m_afterLabelDims, pAfterMins, pAfterRanges);
}

// virtual
void GFilter::trainIncremental(const double* pIn, const double* pOut)
{
	if(m_pFeatureTransform)
	{
		m_pFeatureTransform->transform(pIn, m_pScratchVector);
		if(m_pLabelTransform)
		{
			m_pLabelTransform->transform(pOut, m_pScratchVector2);
			((GIncrementalLearner*)m_pModeler)->trainIncremental(m_pScratchVector, m_pScratchVector2);
		}
		else
			((GIncrementalLearner*)m_pModeler)->trainIncremental(m_pScratchVector, pOut);
	}
	else
	{
		if(m_pLabelTransform)
		{
			m_pLabelTransform->transform(pOut, m_pScratchVector2);
			((GIncrementalLearner*)m_pModeler)->trainIncremental(pIn, m_pScratchVector2);
		}
		else
			((GIncrementalLearner*)m_pModeler)->trainIncremental(pIn, pOut);
	}
}

// virtual
void GFilter::trainSparse(GSparseMatrix* pData, int labelDims)
{
	ThrowError("Sorry, GFilter does not support the trainSparse method due to ambiguities.");
}

// ---------------------------------------------------------------

GPCA::GPCA(int targetDims, GRand* pRand)
: GIncrementalTransform(), m_targetDims(targetDims), m_pBasisVectors(NULL), m_pEigVals(NULL), m_pRand(pRand)
{
}

GPCA::GPCA(GTwtNode* pNode, GRand* pRand)
: GIncrementalTransform(pNode), m_pEigVals(NULL), m_pRand(pRand)
{
	m_targetDims = (int)pNode->field("dims")->asInt();
	m_pRelationAfter = new GUniformRelation(m_targetDims, 0);
	m_pBasisVectors = new GData(pNode->field("basis"));
	m_pRelationBefore = new GUniformRelation(m_pBasisVectors->cols(), 0);
}

// virtual
GPCA::~GPCA()
{
	delete(m_pBasisVectors);
	delete[] m_pEigVals;
}

// virtual
GTwtNode* GPCA::toTwt(GTwtDoc* pDoc)
{
	if(!m_pRelationBefore.get())
		ThrowError("train or enableIncrementalTraining must be called before toTwt");
	GTwtNode* pNode = baseTwtNode(pDoc, "GPCA");
	pNode->addField(pDoc, "dims", pDoc->newInt(m_targetDims));
	pNode->addField(pDoc, "basis", m_pBasisVectors->toTwt(pDoc));
	return pNode;
}

void GPCA::computeEigVals()
{
	delete[] m_pEigVals;
	m_pEigVals = new double[m_targetDims];
}

// virtual
void GPCA::train(GData* pData)
{
	m_pRelationBefore = pData->relation();
	if(!m_pRelationBefore->areContinuous(0, m_pRelationBefore->size()))
		ThrowError("GPCA doesn't support nominal values. You should filter with nominaltocat to make them real.");
	delete(m_pBasisVectors);
	m_pBasisVectors = new GData(m_pRelationBefore);
	m_pBasisVectors->newRows(m_targetDims + 1);
	m_pRelationAfter = new GUniformRelation(m_targetDims, 0);

	// Compute the mean
	int nInputDims = m_pRelationBefore->size();
	double* pMean = m_pBasisVectors->row(0);
	int i;
	for(i = 0; i < nInputDims; i++)
		pMean[i] = pData->mean(i);

	// Make a copy of the data
	GData tmpData(pData->relation(), pData->heap());
	tmpData.copy(pData);

	// Compute the principle components
	double sse = 0;
	if(m_pEigVals)
		sse = tmpData.sumSquaredDistance(pMean);
	for(i = 0; i < m_targetDims; i++)
	{
		double* pVector = m_pBasisVectors->row(i + 1);
		tmpData.principalComponentIgnoreUnknowns(pVector, nInputDims, pMean, m_pRand);
		tmpData.removeComponent(pMean, pVector, nInputDims);
		if(m_pEigVals)
		{
			double t = tmpData.sumSquaredDistance(pMean);
			m_pEigVals[i] = (sse - t) / nInputDims;
			sse = t;
		}
	}
}

// virtual
void GPCA::enableIncrementalTraining(sp_relation& pRelation, double* pMins, double* pRanges)
{
	ThrowError("Sorry, PCA does not support incremental training.");
}

// virtual
void GPCA::transform(const double* pIn, double* pOut)
{
	double* pMean = m_pBasisVectors->row(0);
	int nInputDims = m_pRelationBefore->size();
	for(int i = 0; i < m_targetDims; i++)
	{
		double* pBasisVector = m_pBasisVectors->row(i + 1);
		pOut[i] = GVec::dotProductIgnoringUnknowns(pMean, pIn, pBasisVector, nInputDims);
	}
}

void GPCA::reverse(const double* pIn, double* pOut)
{
	int nInputDims = m_pRelationBefore->size();
	GVec::copy(pOut, m_pBasisVectors->row(0), nInputDims);
	for(int i = 0; i < m_targetDims; i++)
		GVec::addScaled(pOut, pIn[i], m_pBasisVectors->row(i + 1), nInputDims);
}

/*
GPCARotateOnly::GPCARotateOnly(GArffRelation* pRelation, GData* pData)
{
	m_pRelation = pRelation;
	m_pInputData = pData;
	m_pOutputData = NULL;
}

GPCARotateOnly::~GPCARotateOnly()
{
	delete(m_pOutputData);
}

// static
GData* GPCARotateOnly::DoPCA(GArffRelation* pRelation, GData* pData)
{
	GPCA pca(pRelation, pData);
	pca.DoPCA();
	return pca.ReleaseOutputData();
}

void GPCARotateOnly::DoPCA()
{
	// Compute the eigenvectors
	GMatrix m;
	m_pInputData->ComputeCovarianceMatrix(&m, m_pRelation);
	GMatrix eigenVectors;
	eigenVectors.ComputeEigenVectors(m.GetColumnCount(), &m);
	m_pOutputData = new GData(m_pInputData->rows());
	int nRowCount = m_pInputData->rows();
	int nInputCount = m_pRelation->GetInputCount();
	int nOutputCount = m_pRelation->GetOutputCount();
	int nAttributeCount = m_pRelation->size();
	double* pInputRow;
	double* pOutputRow;
	int n, i, j, nIndex;

	// Allocate space for the output
	for(n = 0; n < nRowCount; n++)
	{
		pOutputRow = new double[nAttributeCount];
		m_pOutputData->AddRow(pOutputRow);
	}

	// Compute the output
	double* pEigenVector;
	Holder<double> hInputVector(new double[nInputCount]);
	double* pInputVector = hInputVector.Get();
	for(i = 0; i < nInputCount; i++)
	{
		nIndex = m_pRelation->GetInputIndex(i);
		pEigenVector = eigenVectors.row(i);
		for(n = 0; n < nRowCount; n++)
		{
			pInputRow = m_pInputData->row(n);
			for(j = 0; j < nInputCount; j++)
				pInputVector[j] = pInputRow[m_pRelation->GetInputIndex(j)];
			pOutputRow = m_pOutputData->row(n);
			pOutputRow[nIndex] = GVec::dotProduct(pInputVector, pEigenVector, nInputCount);
		}
	}
	for(i = 0; i < nOutputCount; i++)
	{
		for(n = 0; n < nRowCount; n++)
		{
			nIndex = m_pRelation->GetOutputIndex(i);
			pInputRow = m_pInputData->row(n);
			pOutputRow = m_pOutputData->row(n);
			pOutputRow[nIndex] = pInputRow[nIndex];
		}
	}
}

GData* GPCARotateOnly::ReleaseOutputData()
{
	GData* pData = m_pOutputData;
	m_pOutputData = NULL;
	return pData;
}
*/
GData* GPCARotateOnly::transform(int nDims, int nOutputs, GData* pData, int nComponents, GRand* pRand)
{
	// Init the basis vectors
	int nElements = nDims * nDims;
	double* pBasisVectors = new double[nElements + nDims * 4];
	ArrayHolder<double> hBasisVectors(pBasisVectors);
	double* pComponent = &pBasisVectors[nElements];
	double* pA = &pBasisVectors[nElements + nDims];
	double* pB = &pBasisVectors[nElements + 2 * nDims];
	double* pMean = &pBasisVectors[nElements + 3 * nDims];
	int j;
	for(int i = 0; i < nElements; i++)
		pBasisVectors[i] = 0;
	for(int i = 0; i < nDims; i++)
		pBasisVectors[nDims * i + i] = 1;

	// Compute the mean
	for(j = 0; j < nDims; j++)
		pMean[j] = pData->mean(j);

	// Make a copy of the data
	GData* pOutData = new GData(pData->relation());
	pOutData->copy(pData);
	Holder<GData> hOutData(pOutData);

	// Rotate the basis vectors
	double dDotProd;
	for(int i = 0; i < nComponents; i++)
	{
		// Compute the next principle component
		pOutData->principalComponent(pComponent, nDims, pMean, pRand);
		pOutData->removeComponent(pMean, pComponent, nDims);

		// Use the current axis as the first plane vector
		GVec::copy(pA, &pBasisVectors[nDims * i], nDims);

		// Use the modified Gram-Schmidt process to compute the other plane vector
		GVec::copy(pB, pComponent, nDims);
		dDotProd = GVec::dotProduct(pB, pA, nDims);
		GVec::addScaled(pB, -dDotProd, pA, nDims);
		double dMag = sqrt(GVec::squaredMagnitude(pB, nDims));
		if(dMag < 1e-6)
			break; // It's already close enough. If we normalized something that small, it would just mess up our data
		GVec::multiply(pB, 1.0 / dMag, nDims);

		// Rotate the remaining basis vectors
		double dAngle = atan2(GVec::dotProduct(pComponent, pB, nDims), dDotProd);
		for(j = i; j < nDims; j++)
		{
			GVec::rotate(&pBasisVectors[nDims * j], nDims, dAngle, pA, pB);
			GAssert(ABS(GVec::squaredMagnitude(&pBasisVectors[nDims * j], nDims) - 1.0) < 1e-4);
		}
	}

	// Align data with new basis vectors
	double* pInVector;
	double* pOutVector;
	size_t nCount = pData->rows();
	for(size_t i = 0; i < nCount; i++)
	{
		pInVector = pData->row(i);
		pOutVector = pOutData->row(i);
		for(j = 0; j < nDims; j++)
			pOutVector[j] = GVec::dotProduct(pMean, pInVector, &pBasisVectors[nDims * j], nDims);
	}

	return hOutData.release();
}

#ifndef NO_TEST_CODE
//static
void GPCARotateOnly::test()
{
	GRand prng(0);
	GHeap heap(1000);
	GData data(2, &heap);
	double* pVec;
	pVec = data.newRow();	pVec[0] = 0;	pVec[1] = 0;
	pVec = data.newRow();	pVec[0] = 10;	pVec[1] = 10;
	pVec = data.newRow();	pVec[0] = 4;	pVec[1] = 6;
	pVec = data.newRow();	pVec[0] = 6;	pVec[1] = 4;
	GData* pOut2 = GPCARotateOnly::transform(2, 0, &data, 2, &prng);
	for(size_t i = 0; i < pOut2->rows(); i++)
	{
		pVec = pOut2->row(i);
		if(ABS(ABS(pVec[0]) - 7.071067) < .001)
		{
			if(ABS(pVec[1]) > .001)
				ThrowError("wrong answer");
		}
		else if(ABS(pVec[0]) < .001)
		{
			if(ABS(ABS(pVec[1]) - 1.414214) > .001)
				ThrowError("wrong answer");
		}
		else
			ThrowError("wrong answer");
	}
	delete(pOut2);
}
#endif // !NO_TEST_CODE

// --------------------------------------------------------------------------

GNoiseGenerator::GNoiseGenerator(GRand* pRand)
: GIncrementalTransform(), m_pRand(pRand), m_mean(0), m_deviation(1)
{
}

GNoiseGenerator::GNoiseGenerator(GTwtNode* pNode, GRand* pRand)
: GIncrementalTransform(pNode), m_pRand(pRand)
{
	m_mean = pNode->field("mean")->asDouble();
	m_deviation = pNode->field("dev")->asDouble();
	m_pRelationBefore = GRelation::fromTwt(pNode->field("relation"));
	m_pRelationAfter = m_pRelationBefore;
}

GNoiseGenerator::~GNoiseGenerator()
{
}

// virtual
GTwtNode* GNoiseGenerator::toTwt(GTwtDoc* pDoc)
{
	if(!m_pRelationBefore.get())
		ThrowError("train or enableIncrementalTraining must be called before toTwt");
	GTwtNode* pNode = baseTwtNode(pDoc, "GNoiseGenerator");
	pNode->addField(pDoc, "mean", pDoc->newDouble(m_mean));
	pNode->addField(pDoc, "dev", pDoc->newDouble(m_deviation));
	pNode->addField(pDoc, "relation", m_pRelationBefore->toTwt(pDoc));
	return pNode;
}

// virtual
void GNoiseGenerator::train(GData* pData)
{
	m_pRelationBefore = pData->relation();
	m_pRelationAfter = m_pRelationBefore;
}

// virtual
void GNoiseGenerator::enableIncrementalTraining(sp_relation& pRelation, double* pMins, double* pRanges)
{
	m_pRelationBefore = pRelation;
	m_pRelationAfter = m_pRelationBefore;
	int attrs = m_pRelationBefore->size();
	delete[] m_pAfterMins;
	m_pAfterMins = new double[attrs * 2];
	m_pAfterRanges = m_pAfterMins + attrs;
	GVec::copy(m_pAfterMins, pMins, attrs);
	GVec::copy(m_pAfterRanges, pRanges, attrs);
}

// virtual
void GNoiseGenerator::transform(const double* pIn, double* pOut)
{
	int nDims = m_pRelationBefore->size();
	int i, vals;
	for(i = 0; i < nDims; i++)
	{
		vals = m_pRelationBefore->valueCount(i);
		if(vals == 0)
			pOut[i] = m_pRand->normal() * m_deviation + m_mean;
		else
			pOut[i] = (double)m_pRand->next(vals);
	}
}

// --------------------------------------------------------------------------

GPairProduct::GPairProduct(int nMaxDims)
: GIncrementalTransform(), m_maxDims(nMaxDims)
{
}

GPairProduct::GPairProduct(GTwtNode* pNode)
: GIncrementalTransform(pNode)
{
	m_maxDims = (int)pNode->field("maxDims")->asInt();
	int nAttrsOut = (int)pNode->field("attrs")->asInt();
	m_pRelationAfter = new GUniformRelation(nAttrsOut, 0);
	m_pRelationBefore = GRelation::fromTwt(pNode->field("before"));
}

GPairProduct::~GPairProduct()
{
}

// virtual
GTwtNode* GPairProduct::toTwt(GTwtDoc* pDoc)
{
	if(!m_pRelationBefore.get())
		ThrowError("train or enableIncrementalTraining must be called before toTwt");
	GTwtNode* pNode = baseTwtNode(pDoc, "GPairProduct");
	pNode->addField(pDoc, "before", m_pRelationBefore->toTwt(pDoc));
	pNode->addField(pDoc, "attrs", pDoc->newInt(m_pRelationAfter->size()));
	pNode->addField(pDoc, "maxDims", pDoc->newInt(m_maxDims));
	return pNode;
}

// virtual
void GPairProduct::train(GData* pData)
{
	m_pRelationBefore = pData->relation();
	int nAttrsIn = m_pRelationBefore->size();
	int nAttrsOut = MIN(m_maxDims, nAttrsIn * (nAttrsIn - 1) / 2);
	m_pRelationAfter = new GUniformRelation(nAttrsOut, 0);
}

// virtual
void GPairProduct::enableIncrementalTraining(sp_relation& pRelation, double* pMins, double* pRanges)
{
	m_pRelationBefore = pRelation;
	int nAttrsIn = m_pRelationBefore->size();
	int nAttrsOut = MIN(m_maxDims, nAttrsIn * (nAttrsIn - 1) / 2);
	m_pRelationAfter = new GUniformRelation(nAttrsOut, 0);
	delete[] m_pAfterMins;
	m_pAfterMins = new double[2 * nAttrsOut];
	m_pAfterRanges = m_pAfterMins + nAttrsOut;
	int pos = 0;
	for(int j = 0; j < nAttrsIn && pos < nAttrsOut; j++)
	{
		for(int i = j + 1; i < nAttrsIn && pos < nAttrsOut; i++)
		{
			m_pAfterMins[pos] = MIN(MIN(pMins[i], pMins[j]), pMins[i] * pMins[j]);
			m_pAfterRanges[pos] = MAX(MAX(pMins[i] + pRanges[i], pMins[j] + pRanges[j]), (pMins[i] + pRanges[i]) * (pMins[j] + pRanges[j])) - m_pAfterMins[pos];
			pos++;
		}
	}
}

// virtual
void GPairProduct::transform(const double* pIn, double* pOut)
{
	int i, j, nAttr;
	int nAttrsIn = m_pRelationBefore->size();
	int nAttrsOut = m_pRelationAfter->size();
	nAttr = 0;
	for(j = 0; j < nAttrsIn && nAttr < nAttrsOut; j++)
	{
		for(i = j + 1; i < nAttrsIn && nAttr < nAttrsOut; i++)
			pOut[nAttr++] = pIn[i] * pIn[j];
	}
	GAssert(nAttr == nAttrsOut);
}

// --------------------------------------------------------------------------

GAttributeSelector::GAttributeSelector(GTwtNode* pNode, GRand* pRand)
: GIncrementalTransform(pNode), m_pRand(pRand)
{
	m_labelDims = (int)pNode->field("labels")->asInt();
	m_targetFeatures = (int)pNode->field("target")->asInt();
	GTwtNode* pRanksNode = pNode->field("ranks");
	m_ranks.reserve(pRanksNode->itemCount());
	for(unsigned int i = 0; i < pRanksNode->itemCount(); i++)
		m_ranks.push_back((int)pRanksNode->item(i)->asInt());
	if(m_ranks.size() + (size_t)m_labelDims != (size_t)m_pRelationBefore->size())
		ThrowError("invalid attribute selector");
	if((size_t)m_targetFeatures > m_ranks.size())
		ThrowError("invalid attribute selector");
}

// virtual
GTwtNode* GAttributeSelector::toTwt(GTwtDoc* pDoc)
{
	if(!m_pRelationBefore.get())
		ThrowError("train or enableIncrementalTraining must be called before toTwt");
	GTwtNode* pNode = baseTwtNode(pDoc, "GAttributeSelector");
	pNode->addField(pDoc, "labels", pDoc->newInt(m_labelDims));
	pNode->addField(pDoc, "target", pDoc->newInt(m_targetFeatures));
	GTwtNode* pRanksNode = pNode->addField(pDoc, "ranks", pDoc->newList(m_ranks.size()));
	for(size_t i = 0; i < m_ranks.size(); i++)
		pRanksNode->setItem(i, pDoc->newInt(m_ranks[i]));
	return pNode;
}

void GAttributeSelector::setTargetFeatures(int n)
{
	if(n < 0 || n > m_pRelationBefore->size())
		ThrowError("out of range");
	GMixedRelation* pRelAfter;
	if(m_pRelationBefore->type() == GRelation::ARFF)
		pRelAfter = new GArffRelation();
	else
		pRelAfter = new GMixedRelation();
	for(int i = 0; i < m_targetFeatures; i++)
		pRelAfter->copyAttr(m_pRelationBefore.get(), m_ranks[i]);
	int featureDims = m_pRelationBefore->size() - m_labelDims;
	for(int i = 0; i < m_labelDims; i++)
		pRelAfter->copyAttr(m_pRelationBefore.get(), featureDims + i);
	m_pRelationAfter = pRelAfter;
}

// virtual
void GAttributeSelector::train(GData* pData)
{
	m_pRelationBefore = pData->relation();
	int curDims = pData->cols() - m_labelDims;
	m_ranks.resize(curDims);
	GData* pClone = pData->clone();
	Holder<GData> hClone(pClone);
	vector<int> indexMap;
	for(int i = 0; i < curDims; i++)
		indexMap.push_back(i);

	// Produce a ranked attributed ordering by deselecting the weakest attribute each time
	while(curDims > 1)
	{
		// Train a single-layer neural network with the normalized remaining data
		GNeuralNet nn(m_pRand);
		nn.setIterationsPerValidationCheck(20);
		nn.setMinImprovement(0.002);
		GNominalToCat cat;
		GFilter fil1(&nn, false);
		fil1.setFeatureTransform(&cat, false);
		fil1.setLabelTransform(new GNominalToCat(), true);
		GFilter fil(&fil1, false);
		fil.setFeatureTransform(new GNormalize(-2.0, 2.0), true);
		fil.setLabelTransform(new GNormalize(0.0, 1.0), true);
		fil.train(pClone, m_labelDims);
		vector<int> rmap;
		cat.reverseAttrMap(rmap);

		// Identify the weakest attribute
		GNeuralNetLayer& layer = nn.getLayer(nn.layerCount() - 1);
		int pos = 0;
		double weakest = 0;
		int weakestIndex = -1;
		for(int i = 0; i < curDims; i++)
		{
			double w = 0;
			while(pos < nn.featureDims() && rmap[pos] == i)
			{
				for(vector<GNeuron>::iterator it = layer.m_neurons.begin(); it != layer.m_neurons.end(); it++)
					w = MAX(w, ABS(it->m_weights[pos + 1]));
				pos++;
			}
			if(weakestIndex < 0 || w < weakest)
			{
				weakest = w;
				weakestIndex = i;
			}
		}

		// Deselect the weakest attribute
		m_ranks[curDims - 1] = indexMap[weakestIndex];
		indexMap.erase(indexMap.begin() + weakestIndex);
		pClone->deleteColumn(weakestIndex);
		curDims--;
		GAssert(pClone->cols() - m_labelDims == curDims);
	}
	m_ranks[0] = indexMap[0];
	setTargetFeatures(m_targetFeatures);
}

// virtual
void GAttributeSelector::enableIncrementalTraining(sp_relation& pRelation, double* pMins, double* pRanges)
{
	ThrowError("not implemented yet");
}

// virtual
void GAttributeSelector::transform(const double* pIn, double* pOut)
{
	int i;
	for(i = 0; i < m_targetFeatures; i++)
		pOut[i] = pIn[m_ranks[i]];
	int featureDims = m_pRelationBefore->size() - m_labelDims;
	for(int j = 0; j < m_labelDims; j++)
		pOut[i++] = pIn[featureDims + j];
}

#ifndef NO_TEST_CODE
//static
void GAttributeSelector::test()
{
	GRand prng(0);
	GData data(21);
	for(size_t i = 0; i < 256; i++)
	{
		double* pVec = data.newRow();
		prng.cubical(pVec, 20);
		pVec[20] = 0.2 * pVec[3] * pVec[3] * - 7.0 * pVec[3] * pVec[13] + pVec[17];
	}
	GAttributeSelector as(1, 3, &prng);
	as.train(&data);
	std::vector<int>& r = as.ranks();
	if(r[1] == r[0] || r[2] == r[0] || r[2] == r[1])
		ThrowError("bogus rankings");
	if(r[0] != 3 && r[0] != 13 && r[0] != 17)
		ThrowError("failed");
	if(r[1] != 3 && r[1] != 13 && r[1] != 17)
		ThrowError("failed");
	if(r[2] != 3 && r[2] != 13 && r[2] != 17)
		ThrowError("failed");
}
#endif // NO_TEST_CODE
// --------------------------------------------------------------------------

GNominalToCat::GNominalToCat(int nValueCap)
: GTwoWayIncrementalTransform(), m_valueCap(nValueCap)
{
}

GNominalToCat::GNominalToCat(GTwtNode* pNode)
: GTwoWayIncrementalTransform(pNode)
{
	m_valueCap = (int)pNode->field("valueCap")->asInt();
	m_pRelationBefore = GRelation::fromTwt(pNode->field("before"));
	init(m_pRelationBefore);
}

void GNominalToCat::init(sp_relation& pRelation)
{
	m_pRelationBefore = pRelation;
	if(m_pRelationBefore->type() == GRelation::ARFF)
		m_pRelationAfter = new GArffRelation();
	else
		m_pRelationAfter = new GMixedRelation();
	int nDims = 0;
	int nAttrCount = m_pRelationBefore->size();
	int i, j, nValues;
	const char* szName;
	for(i = 0; i < nAttrCount; i++)
	{
		nValues = m_pRelationBefore->valueCount(i);
		if(nValues < 3 || nValues >= m_valueCap)
		{
			nDims++;
			if(m_pRelationBefore->type() == GRelation::ARFF)
			{
				szName = ((GArffRelation*)m_pRelationBefore.get())->attrName(i);
				((GArffRelation*)m_pRelationAfter.get())->addAttribute(szName, 0, NULL);
			}
			else
				((GMixedRelation*)m_pRelationAfter.get())->addAttr(0);
		}
		else
		{
			nDims += nValues;
			if(m_pRelationBefore->type() == GRelation::ARFF)
			{
				szName = ((GArffRelation*)m_pRelationBefore.get())->attrName(i);
				string sName = szName;
				sName += "_";
				for(j = 0; j < nValues; j++)
				{
					string s;
					m_pRelationBefore->attrValue(&s, i, j);
					sName += s;
					((GArffRelation*)m_pRelationAfter.get())->addAttribute(s.c_str(), 0, NULL);
				}
			}
			else
			{
				for(j = 0; j < nValues; j++)
					((GMixedRelation*)m_pRelationAfter.get())->addAttr(0);
			}
		}
	}
}

// virtual
GNominalToCat::~GNominalToCat()
{
}

// virtual
void GNominalToCat::train(GData* pData)
{
	init(pData->relation());
}

// virtual
void GNominalToCat::enableIncrementalTraining(sp_relation& pRelation, double* pMins, double* pRanges)
{
	init(pRelation);
	int attrs = m_pRelationAfter->size();
	delete[] m_pAfterMins;
	m_pAfterMins = new double[2 * attrs];
	m_pAfterRanges = m_pAfterMins + attrs;
	int nAttrCount = m_pRelationBefore->size();
	int i, nValues;
	int pos = 0;
	for(i = 0; i < nAttrCount; i++)
	{
		nValues = m_pRelationBefore->valueCount(i);
		if(nValues < 3 || nValues >= m_valueCap)
		{
			m_pAfterMins[pos] = pMins[i];
			m_pAfterRanges[pos] = pRanges[i];
			pos++;
		}
		else
		{
			GVec::setAll(m_pAfterMins + pos, 0.0, nValues);
			GVec::setAll(m_pAfterRanges + pos, 1.0, nValues);
			pos += nValues;
		}
	}

}

// virtual
GTwtNode* GNominalToCat::toTwt(GTwtDoc* pDoc)
{
	if(!m_pRelationBefore.get())
		ThrowError("train or enableIncrementalTraining must be called before toTwt");
	GTwtNode* pNode = baseTwtNode(pDoc, "GNominalToCat");
	pNode->addField(pDoc, "valueCap", pDoc->newInt(m_valueCap));
	pNode->addField(pDoc, "before", m_pRelationBefore->toTwt(pDoc));
	return pNode;
}

// virtual
void GNominalToCat::transform(const double* pIn, double* pOut)
{
	int nInAttrCount = m_pRelationBefore->size();
	for(int i = 0; i < nInAttrCount; i++)
	{
		int nValues = m_pRelationBefore->valueCount(i);
		if(nValues < 3 || nValues >= m_valueCap)
		{
			if(nValues == 0)
				*(pOut++) = *(pIn++);
			else
			{
				if(*pIn < 0)
					*(pOut++) = 0.5;
				else
					*(pOut++) = *pIn;
				pIn++;
			}
		}
		else
		{
			if(*pIn >= 0)
			{
				GAssert(*pIn < nValues);
				GVec::setAll(pOut, 0.0, nValues);
				pOut[(int)*pIn] = 1.0;
			}
			else
				GVec::setAll(pOut, 1.0 / nValues, nValues);
			pOut += nValues;
			pIn++;
		}
	}
}

// virtual
void GNominalToCat::untransform(const double* pIn, double* pOut)
{
	int nOutAttrCount = m_pRelationBefore->size();
	for(int i = 0; i < nOutAttrCount; i++)
	{
		int nValues = m_pRelationBefore->valueCount(i);
		if(nValues < 3 || nValues >= m_valueCap)
		{
			if(nValues == 0)
				*(pOut++) = *(pIn++);
			else
				*(pOut++) = (*(pIn++) < 0.5 ? 0 : 1);
		}
		else
		{
			double max = *(pIn++);
			*pOut = 0.0;
			for(int i = 1; i < nValues; i++)
			{
				if(*pIn > max)
				{
					max = *pIn;
					*pOut = (double)i;
				}
				pIn++;
			}
		}
	}
}

void GNominalToCat::reverseAttrMap(vector<int>& rmap)
{
	rmap.clear();
	int nInAttrCount = m_pRelationBefore->size();
	for(int i = 0; i < nInAttrCount; i++)
	{
		int nValues = m_pRelationBefore->valueCount(i);
		if(nValues < 3 || nValues >= m_valueCap)
			rmap.push_back(i);
		else
		{
			for(int j = 0; j < nValues; j++)
				rmap.push_back(i);
		}
	}
}

// --------------------------------------------------------------------------

GNormalize::GNormalize(double min, double max)
: GTwoWayIncrementalTransform(), m_min(min), m_max(max), m_pMins(NULL), m_pRanges(NULL)
{
}

GNormalize::GNormalize(GTwtNode* pNode)
: GTwoWayIncrementalTransform(pNode)
{
	m_pRelationBefore = GRelation::fromTwt(pNode->field("relation"));
	m_pRelationAfter = m_pRelationBefore;
	m_min = pNode->field("min")->asDouble();
	m_max = pNode->field("max")->asDouble();
	int nAttrCount = m_pRelationBefore->size();
	m_pMins = new double[2 * nAttrCount];
	m_pRanges = &m_pMins[nAttrCount];
	GVec::fromTwt(m_pMins, nAttrCount, pNode->field("mins"));
	GVec::fromTwt(m_pRanges, nAttrCount, pNode->field("ranges"));
}

// virtual
GNormalize::~GNormalize()
{
	delete[] m_pMins;
}

// virtual
GTwtNode* GNormalize::toTwt(GTwtDoc* pDoc)
{
	if(!m_pRelationBefore.get())
		ThrowError("train or enableIncrementalTraining must be called before toTwt");
	GTwtNode* pNode = baseTwtNode(pDoc, "GNormalize");
	pNode->addField(pDoc, "relation", m_pRelationBefore->toTwt(pDoc));
	pNode->addField(pDoc, "min", pDoc->newDouble(m_min));
	pNode->addField(pDoc, "max", pDoc->newDouble(m_max));
	int nAttrCount = m_pRelationBefore->size();
	pNode->addField(pDoc, "mins", GVec::toTwt(pDoc, m_pMins, nAttrCount));
	pNode->addField(pDoc, "ranges", GVec::toTwt(pDoc, m_pRanges, nAttrCount));
	return pNode;
}

void GNormalize::setMinsAndRanges(sp_relation& pRel, const double* pMins, const double* pRanges)
{
	m_pRelationBefore = pRel;
	m_pRelationAfter = m_pRelationBefore;
	int nAttrCount = m_pRelationBefore->size();
	delete[] m_pMins;
	m_pMins = new double[2 * nAttrCount];
	m_pRanges = &m_pMins[nAttrCount];
	GVec::copy(m_pMins, pMins, nAttrCount);
	GVec::copy(m_pRanges, pRanges, nAttrCount);
}

// virtual
void GNormalize::train(GData* pData)
{
	m_pRelationBefore = pData->relation();
	m_pRelationAfter = m_pRelationBefore;
	int nAttrCount = m_pRelationBefore->size();
	delete[] m_pMins;
	m_pMins = new double[2 * nAttrCount];
	m_pRanges = &m_pMins[nAttrCount];
	for(int i = 0; i < nAttrCount; i++)
	{
		if(m_pRelationBefore->valueCount(i) == 0)
		{
			pData->minAndRangeUnbiased(i, &m_pMins[i], &m_pRanges[i]);
			if(m_pRanges[i] < 1e-12)
				m_pRanges[i] = 1.0;
		}
		else
		{
			m_pMins[i] = 0;
			m_pRanges[i] = 0;
		}
	}
}

// virtual
void GNormalize::enableIncrementalTraining(sp_relation& pRelation, double* pMins, double* pRanges)
{
	m_pRelationBefore = pRelation;
	m_pRelationAfter = m_pRelationBefore;
	int nAttrCount = m_pRelationBefore->size();
	delete[] m_pMins;
	m_pMins = new double[2 * nAttrCount];
	m_pRanges = &m_pMins[nAttrCount];
	GVec::copy(m_pMins, pMins, nAttrCount);
	GVec::copy(m_pRanges, pRanges, nAttrCount);
	for(int i = 0; i < nAttrCount; i++)
	{
		if(pRelation->valueCount(i) == 0)
		{
			if(m_pRanges[i] < 1e-12)
				m_pRanges[i] = 1.0;
		}
	}
	delete[] m_pAfterMins;
	m_pAfterMins = new double[2 * nAttrCount];
	m_pAfterRanges = m_pAfterMins + nAttrCount;
	GVec::setAll(m_pAfterMins, m_min, nAttrCount);
	GVec::setAll(m_pAfterRanges, m_max - m_min, nAttrCount);
}

// virtual
void GNormalize::transform(const double* pIn, double* pOut)
{
	int nAttrCount = m_pRelationBefore->size();
	double* pMins = m_pMins;
	double* pRanges = m_pRanges;
	for(int i = 0; i < nAttrCount; i++)
	{
		if(m_pRelationBefore->valueCount(i) == 0)
			*pOut = GData::normalize(*pIn, *pMins, *pRanges, m_min, m_max - m_min);
		else
			*pOut = *pIn;
		pOut++;
		pIn++;
		pMins++;
		pRanges++;
	}
}

// virtual
void GNormalize::untransform(const double* pIn, double* pOut)
{
	int nAttrCount = m_pRelationBefore->size();
	double* pMins = m_pMins;
	double* pRanges = m_pRanges;
	for(int i = 0; i < nAttrCount; i++)
	{
		if(m_pRelationBefore->valueCount(i) == 0)
			*pOut = GData::normalize(*pIn, m_min, m_max - m_min, *pMins, *m_pRanges);
		else
			*pOut = *pIn;
		pOut++;
		pIn++;
		pMins++;
		pRanges++;
	}
}

// --------------------------------------------------------------------------

GDiscretize::GDiscretize(int buckets)
: GTwoWayIncrementalTransform()
{
	m_bucketsIn = buckets;
	m_bucketsOut = -1;
	m_pMins = NULL;
	m_pRanges = NULL;
}

GDiscretize::GDiscretize(GTwtNode* pNode)
: GTwoWayIncrementalTransform(pNode)
{
	m_bucketsIn = (int)pNode->field("bucketsIn")->asInt();
	m_bucketsOut = (int)pNode->field("bucketsOut")->asInt();
	m_pRelationBefore = GRelation::fromTwt(pNode->field("before"));
	m_pRelationAfter = GRelation::fromTwt(pNode->field("after"));
	int nAttrCount = m_pRelationBefore->size();
	m_pMins = new double[2 * nAttrCount];
	m_pRanges = &m_pMins[nAttrCount];
	GVec::fromTwt(m_pMins, nAttrCount, pNode->field("mins"));
	GVec::fromTwt(m_pRanges, nAttrCount, pNode->field("ranges"));
}

// virtual
GDiscretize::~GDiscretize()
{
	delete[] m_pMins;
}

// virtual
GTwtNode* GDiscretize::toTwt(GTwtDoc* pDoc)
{
	if(!m_pRelationBefore.get())
		ThrowError("train or enableIncrementalTraining must be called before toTwt");
	GTwtNode* pNode = baseTwtNode(pDoc, "GDiscretize");
	pNode->addField(pDoc, "before", m_pRelationBefore->toTwt(pDoc));
	pNode->addField(pDoc, "after", m_pRelationAfter->toTwt(pDoc));
	pNode->addField(pDoc, "bucketsIn", pDoc->newInt(m_bucketsIn));
	pNode->addField(pDoc, "bucketsOut", pDoc->newInt(m_bucketsOut));
	int nAttrCount = m_pRelationBefore->size();
	pNode->addField(pDoc, "mins", GVec::toTwt(pDoc, m_pMins, nAttrCount));
	pNode->addField(pDoc, "ranges", GVec::toTwt(pDoc, m_pRanges, nAttrCount));
	return pNode;
}

// virtual
void GDiscretize::train(GData* pData)
{
	// Make the relations
	m_pRelationBefore = pData->relation();
	m_bucketsOut = m_bucketsIn;
	if(m_bucketsOut < 0)
		m_bucketsOut = MAX(2, (int)sqrt((double)pData->rows()));
	int nAttrCount = m_pRelationBefore->size();
	GMixedRelation* pRelationAfter = new GMixedRelation();
	m_pRelationAfter = pRelationAfter;
	for(int i = 0; i < nAttrCount; i++)
	{
		int nValues = m_pRelationBefore->valueCount(i);
		if(nValues > 0)
			pRelationAfter->addAttr(nValues);
		else
			pRelationAfter->addAttr(m_bucketsOut);
	}

	// Determine the boundaries
	delete[] m_pMins;
	m_pMins = new double[2 * nAttrCount];
	m_pRanges = &m_pMins[nAttrCount];
	int i, nValues;
	for(i = 0; i < nAttrCount; i++)
	{
		nValues = m_pRelationBefore->valueCount(i);
		if(nValues > 0)
		{
			m_pMins[i] = 0;
			m_pRanges[i] = 0;
		}
		else
			pData->minAndRangeUnbiased(i, &m_pMins[i], &m_pRanges[i]);
	}
}

// virtual
void GDiscretize::enableIncrementalTraining(sp_relation& pRelation, double* pMins, double* pRanges)
{
	// Make the relations
	m_pRelationBefore = pRelation;
	m_bucketsOut = m_bucketsIn;
	if(m_bucketsOut < 0)
		ThrowError("The number of discretization buckets cannot be automatically determined with incremental training.");
	int nAttrCount = m_pRelationBefore->size();
	GMixedRelation* pRelationAfter = new GMixedRelation();
	m_pRelationAfter = pRelationAfter;
	for(int i = 0; i < nAttrCount; i++)
	{
		int nValues = m_pRelationBefore->valueCount(i);
		if(nValues > 0)
			pRelationAfter->addAttr(nValues);
		else
			pRelationAfter->addAttr(m_bucketsOut);
	}

	// Determine the boundaries
	delete[] m_pMins;
	m_pMins = new double[2 * nAttrCount];
	m_pRanges = &m_pMins[nAttrCount];
	int i, nValues;
	for(i = 0; i < nAttrCount; i++)
	{
		nValues = m_pRelationBefore->valueCount(i);
		if(nValues > 0)
		{
			m_pMins[i] = 0;
			m_pRanges[i] = 0;
		}
		else
		{
			m_pMins[i] = pMins[i];
			m_pRanges[i] = pRanges[i];
		}
	}
}

// virtual
void GDiscretize::transform(const double* pIn, double* pOut)
{
	if(!m_pMins)
		ThrowError("Train was not called");
	int nAttrCount = m_pRelationBefore->size();
	int i, nValues;
	for(i = 0; i < nAttrCount; i++)
	{
		nValues = m_pRelationBefore->valueCount(i);
		if(nValues > 0)
			pOut[i] = pIn[i];
		else
			pOut[i] = MAX(0, MIN(m_bucketsOut - 1, (int)(((pIn[i] - m_pMins[i]) * m_bucketsOut) / m_pRanges[i])));
	}
}

// virtual
void GDiscretize::untransform(const double* pIn, double* pOut)
{
	if(!m_pMins)
		ThrowError("Train was not called");
	int nAttrCount = m_pRelationBefore->size();
	int i, nValues;
	for(i = 0; i < nAttrCount; i++)
	{
		nValues = m_pRelationBefore->valueCount(i);
		if(nValues > 0)
			pOut[i] = pIn[i];
		else
			pOut[i] = (((double)pIn[i] + .5) * m_pRanges[i]) / m_bucketsOut + m_pMins[i];
	}
}

} // namespace GClasses

