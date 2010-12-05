/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#include "GSystemLearner.h"
#include "GActivation.h"
#include "GPolicyLearner.h"
#include "GNeuralNet.h"
#include "GNeighborFinder.h"
#include "GHillClimber.h"
#include <stdlib.h>
#include "GVec.h"
#include "GRand.h"
#include "GManifold.h"
#include "GFile.h"
#include "GHeap.h"
#include "GImage.h"
#include "GOptimizer.h"
#include "GTwt.h"
#include "GTime.h"
#include "GEvolutionary.h"
#include "GApp.h"
#include <deque>
#include <math.h>
#include <iostream>
#include <string>
#include <vector>

namespace GClasses {

using std::string;
using std::deque;
using std::cout;
using std::cerr;
using std::vector;

GTwtNode* GSystemLearner::baseTwtNode(GTwtDoc* pDoc, const char* szClassName)
{
	GTwtNode* pNode = pDoc->newObj();
	pNode->addField(pDoc, "class", pDoc->newString(szClassName));
	return pNode;
}

// ------------------------------------------------------------------------------------------

GRecurrentModel::GRecurrentModel(GSupervisedLearner* pTransition, GSupervisedLearner* pObservation, int actionDims, int contextDims, int obsDims, GRand* pRand, std::vector<size_t>* pParamDims)
: GSystemLearner(), m_actionDims(actionDims), m_contextDims(contextDims), m_obsDims(obsDims)
{
	if(m_actionDims < 0)
		ThrowError("Invalid number of action dims");
	m_pixels = 1;
	int paramDims = pParamDims ? pParamDims->size() : 0;
	for(int i = 0; i < paramDims; i++)
		m_pixels *= (*pParamDims)[i];
	m_channels = m_obsDims / m_pixels;
	if(m_obsDims != m_channels * m_pixels)
		ThrowError("Invalid observation dims");
	m_pTransitionFunc = pTransition;
	m_pObservationFunc = pObservation;
	m_paramDims = paramDims;
	m_pRand = pRand;
	m_pParamRanges = new size_t[paramDims];
	for(size_t i = 0; i < pParamDims->size(); i++)
		m_pParamRanges[i] = (*pParamDims)[i];
	int bufDims = MAX(m_channels + m_contextDims, m_actionDims + m_contextDims + m_contextDims);
	m_pParams = new double[m_paramDims + m_contextDims + bufDims];
	m_pContext = m_pParams + m_paramDims;
	m_pBuf = m_pContext + m_contextDims;
	GVec::setAll(m_pContext, 0.0, m_contextDims);
	m_transitionDelta = true;
	m_useIsomap = false;
	m_trainingSeconds = 60 * 60; // one hour
	m_validationInterval = 0;
	m_pValidationData = NULL;
	m_multiplier = 1.0;
}

GRecurrentModel::GRecurrentModel(GTwtNode* pNode, GRand* pRand)
: GSystemLearner(pNode),
	m_pRand(pRand)
{
	// Load the models
	GLearnerLoader ll;
	m_pTransitionFunc = ll.loadModeler(pNode->field("trans"), pRand);
	m_pObservationFunc = ll.loadModeler(pNode->field("obs"), pRand);

	// Load the param ranges
	GTwtNode* pContext = pNode->field("context");
	m_contextDims = pContext->itemCount();
	if(m_contextDims != m_pTransitionFunc->labelDims())
		ThrowError("invalid model");
	GTwtNode* pParamRanges = pNode->field("params");
	m_paramDims = pParamRanges->itemCount();
	if(m_paramDims != (size_t)(m_pObservationFunc->featureDims() - m_contextDims))
		ThrowError("invalid model");
	m_pParamRanges = new size_t[m_paramDims];
	for(size_t i = 0; i < m_paramDims; i++)
		m_pParamRanges[i] = (size_t)pParamRanges->item(i)->asInt();

	// Infer other stuff
	m_actionDims = m_pTransitionFunc->featureDims() - m_contextDims;
	if(m_actionDims < 0)
		ThrowError("invalid model");
	m_channels = m_pObservationFunc->labelDims();
	m_pixels = 1;
	for(size_t i = 0; i < m_paramDims; i++)
		m_pixels *= m_pParamRanges[i];
	m_obsDims = m_channels * m_pixels;
	m_transitionDelta = pNode->field("delta")->asBool();
	m_trainingSeconds = pNode->field("trainSecs")->asDouble();

	// Load the context
	int bufDims = m_actionDims + m_contextDims + m_contextDims;
	m_pParams = new double[m_paramDims + m_contextDims + bufDims];
	m_pContext = m_pParams + m_paramDims;
	m_pBuf = m_pContext + m_contextDims;
	for(int i = 0; i < m_contextDims; i++)
		m_pContext[i] = pContext->item(i)->asDouble();
	m_validationInterval = 0;
	m_pValidationData = NULL;
	m_multiplier = 1.0;
}

// virtual
GRecurrentModel::~GRecurrentModel()
{
	delete(m_pObservationFunc);
	delete(m_pTransitionFunc);
	delete[] m_pParamRanges;
	delete[] m_pParams;
}

// virtual
GTwtNode* GRecurrentModel::toTwt(GTwtDoc* pDoc)
{
	GTwtNode* pNode = baseTwtNode(pDoc, "GRecurrentModel");
	pNode->addField(pDoc, "trans", m_pTransitionFunc->toTwt(pDoc));
	pNode->addField(pDoc, "obs", m_pObservationFunc->toTwt(pDoc));
	pNode->addField(pDoc, "delta", pDoc->newBool(m_transitionDelta));
	pNode->addField(pDoc, "isomap", pDoc->newBool(m_useIsomap));
	pNode->addField(pDoc, "trainSecs", pDoc->newDouble(m_trainingSeconds));
	GTwtNode* pParamRanges = pNode->addField(pDoc, "params", pDoc->newList(m_paramDims));
	for(size_t i = 0; i < m_paramDims; i++)
		pParamRanges->setItem(i, pDoc->newInt(m_pParamRanges[i]));
	GTwtNode* pContext = pNode->addField(pDoc, "context", pDoc->newList(m_contextDims));
	for(int i = 0; i < m_contextDims; i++)
		pContext->setItem(i, pDoc->newDouble(m_pContext[i]));
	return pNode;
}

void GRecurrentModel::validateDuringTraining(double timeInterval, vector<GData*>* pValidationData)
{
	m_validationInterval = timeInterval;
	m_pValidationData = pValidationData;
}

// static
void GRecurrentModel::blurImageVector(const double* pIn, double* pOut, int wid, int hgt, int chan, double valueRange, int radius, int reps)
{
	GImage image;
	GVec::toImage(pIn, &image, wid, hgt, chan, valueRange);
	image.blurQuick(reps, radius);
	GVec::fromImage(&image, pOut, wid, hgt, chan, valueRange);
}

GData* GRecurrentModel::mosesEstimateState(GData* pActions, GData* pObservations)
{
	// Estimate the state
	int neighbors = 256;
	GPCA pca(12, m_pRand);
	pca.train(pObservations);
	GData* pReducedObs = pca.transformBatch(pObservations);
	Holder<GData> hReducedObs(pReducedObs);
	GDynamicSystemNeighborFinder nf(pReducedObs, pActions, false, neighbors, m_pRand);
//	GTemporalNeighborFinder nf(pReducedObs, neighbors, m_pRand);

	//GNeighborFinderCacheWrapper nf2(&nf, false);
	//nf2.fillCache();
	//nf2.cutShortcuts(8);

	GManifoldLearner* pML;
	if(m_useIsomap)
	{
		pML = new GIsomap(neighbors, m_contextDims, m_pRand);
		((GIsomap*)pML)->setNeighborFinder(&nf);
	}
	else
	{
		pML = new GBreadthFirstUnfolding(1/*reps*/, neighbors, m_contextDims, m_pRand);
		((GBreadthFirstUnfolding*)pML)->setNeighborFinder(&nf);
	}
	Holder<GManifoldLearner> hML(pML);
	GData* pEstState = pML->doit(pObservations);

	// Center the first state estimate at the origin
	double* pRow = pEstState->row(0);
	for(size_t i = 1; i < pEstState->rows(); i++)
		GVec::subtract(pEstState->row(i), pRow, m_contextDims);
	GVec::setAll(pRow, 0.0, m_contextDims);

	return pEstState;
}

void GRecurrentModel::trainTransitionFunction(GData* pActions, GData* pEstState)
{
	// Make data for training the transition function
	GMixedRelation* pRelation = new GMixedRelation();
	sp_relation pRel = pRelation;
	pRelation->addAttrs(pActions->relation().get());
	pRelation->addAttrs(m_contextDims + m_contextDims, 0);
	GData* pTransitionData = new GData(pRel);
	Holder<GData> hTransitionData(pTransitionData);
	pTransitionData->newRows(pActions->rows() - 1);
	for(size_t i = 0; i < pActions->rows() - 1; i++)
	{
		double* pOut = pTransitionData->row(i);
		GVec::copy(pOut, pActions->row(i), m_actionDims);
		GVec::copy(pOut + m_actionDims, pEstState->row(i), m_contextDims);
		GVec::copy(pOut + m_actionDims + m_contextDims, pEstState->row(i + 1), m_contextDims);
		if(m_transitionDelta)
			GVec::subtract(pOut + m_actionDims + m_contextDims, pEstState->row(i), m_contextDims);
	}

	// Train the transition function
	m_pTransitionFunc->train(pTransitionData, m_contextDims);
}

GData* GRecurrentModel::joshuaEstimateState(GData* pActions, GData* pObservations)
{
	// Convert actions to a categorical distribution
	GNominalToCat cat;
	GData* pActionsCat = cat.doit(pActions);
	Holder<GData> hActionsCat(pActionsCat);

	// Use PCA to reduce the data
	int prinComps = 64;//pObservations->countPrincipalComponents(0.01, m_pRand);
	GPCA* pPCA = NULL;
	Holder<GPCA> hPCA(pPCA);
	if(GFile::doesFileExist("joshua_pca.twt"))
	{
		GTwtDoc doc;
		doc.load("joshua_pca.twt");
		pPCA = new GPCA(doc.root(), m_pRand);
		hPCA.reset(pPCA);
	}
	else
	{
		pPCA = new GPCA(prinComps, m_pRand);
		hPCA.reset(pPCA);
		GData clone(pObservations->relation());
		GReleaseDataHolder hClone(&clone);
		size_t step = MAX((size_t)1, pObservations->rows() / 4096);
		for(size_t i = 0; i < pObservations->rows(); i += step)
			clone.takeRow(pObservations->row(i));
		pPCA->train(&clone);
		GTwtDoc doc;
		doc.setRoot(pPCA->toTwt(&doc));
		doc.save("joshua_pca.twt");
	}
	GData* pReducedObs = pPCA->transformBatch(pObservations);
	Holder<GData> hReducedObs(pReducedObs);
//cout << "% eigvals: ";
//double* pEigVals = pca.eigVals();
//GVec::print(cout, 8, pEigVals, 12);
//cout << "\n";

	// Make h
	GFilter* pH = NULL;
	int totalDims = m_paramDims + pActionsCat->cols() + prinComps + m_channels;
	{
		GNeuralNet* pNN = new GNeuralNet(m_pRand);
		pNN->addLayer(15);
		pNN->addLayer(50);

		sp_relation pRel = new GUniformRelation(totalDims, 0);
		double* pMins = new double[2 * (totalDims)];
		ArrayHolder<double> hMins(pMins);
		double* pRanges = pMins + totalDims;
		double* pM = pMins;
		double* pR = pRanges;
		for(size_t i = 0; i < m_paramDims; i++)
		{
			*(pM++) = 0.0;
			*(pR++) = 1.0;
		}
		for(int i = 0; i < pActionsCat->cols(); i++)
		{
			*(pM++) = 0.0;
			*(pR++) = 1.0;
		}
		for(int i = 0; i < prinComps; i++)
			pReducedObs->minAndRangeUnbiased(i, pM++, pR++);
		for(int i = 0; i < m_channels; i++)
		{
			*(pM++) = 0.0;
			*(pR++) = 255.0;
		}

		pH = new GFilter(pNN, true);
		pH->setFeatureTransform(new GNormalize(0.0, 1.0), true);
		pH->setLabelTransform(new GNormalize(0.0, 1.0), true);
		pH->enableIncrementalLearning(pRel, m_channels, pMins, pRanges);
	}

	// Train h
	int samples = 30;
	int trainSeconds = 4 * 60 * 60;
	GTEMPBUF(double, buf1, 2 * totalDims);
	GVec::setAll(buf1 + m_paramDims, 0.0, pActionsCat->cols());
	double* buf2 = buf1 + totalDims;
	double dStart = GTime::seconds();
	GCoordVectorIterator it(m_paramDims, m_pParamRanges);
	while(true)
	{
		double dTime = GTime::seconds() - dStart;
		if(dTime >= trainSeconds)
			break;
		size_t row = (size_t)m_pRand->next(pReducedObs->rows() - 1);
		GVec::copy(buf1 + m_paramDims + pActionsCat->cols(), pReducedObs->row(row), prinComps);
		GVec::copy(buf2 + m_paramDims, pActionsCat->row(row), pActionsCat->cols());
		GVec::copy(buf2 + m_paramDims + pActionsCat->cols(), pReducedObs->row(row), prinComps);
		double* pObsCur = pObservations->row(row);
		double* pObsNext = pObservations->row(row + 1);
		for(int i = 0; i < samples; i++)
		{
			// Pick random parameters
			it.setRandom(m_pRand);
			it.currentNormalized(buf1);
			it.currentNormalized(buf2);
			size_t ofs = m_channels * it.currentIndex();
			pH->trainIncremental(buf2, pObsNext + ofs);
			pH->trainIncremental(buf1, pObsCur + ofs);
			pH->trainIncremental(buf2, pObsNext + ofs);
		}
	}

	// Create a visualization of the mask
	if(m_paramDims != 2)
		ThrowError("Expected image-based observations");
	size_t frameCount = MIN((size_t)20, pObservations->rows());
	GImage image;
	image.setSize(m_pParamRanges[0] * 3, m_pParamRanges[1] * frameCount);
	double* pred1 = buf1 + m_paramDims + pActionsCat->cols() + prinComps;
	double* pred2 = buf2 + m_paramDims + pActionsCat->cols() + prinComps;
	double deltaScale = 1.0;
	for(size_t i = 0; i < frameCount; i++)
	{
		// Make the observed delta image
		int xx = 0;
		int yy = m_pParamRanges[1] * i;
		double* pFrameCur = pObservations->row(i);
		double* pFrameNext = pObservations->row(i + 1);
		for(size_t y = 0; y < m_pParamRanges[1]; y++)
		{
			for(size_t x = 0; x < m_pParamRanges[0]; x++)
			{
				int r = ClipChan(128 + (int)(deltaScale * (*pFrameNext++ - *pFrameCur++)));
				int g = ClipChan(128 + (int)(deltaScale * (*pFrameNext++ - *pFrameCur++)));
				int b = ClipChan(128 + (int)(deltaScale * (*pFrameNext++ - *pFrameCur++)));
				image.setPixel(xx + x, yy + y, gARGB(0xff, r, g, b));
			}
		}

		// Make the predicted masked delta image
		xx = m_pParamRanges[0];
		GVec::copy(buf1 + m_paramDims + pActionsCat->cols(), pReducedObs->row(i), prinComps);
		GVec::copy(buf2 + m_paramDims + pActionsCat->cols(), pReducedObs->row(i), prinComps);
		GCoordVectorIterator pi(m_paramDims, m_pParamRanges);
		for(size_t y = 0; y < m_pParamRanges[1]; y++)
		{
			for(size_t x = 0; x < m_pParamRanges[0]; x++)
			{
				pi.currentNormalized(buf1);
				pi.currentNormalized(buf2);
				GVec::copy(buf2 + m_paramDims, pActionsCat->row(i), pActionsCat->cols());
				pH->predict(buf1, pred1);
				pH->predict(buf2, pred2);
				int r = ClipChan(128 + (int)(deltaScale * (pred2[0] - pred1[0])));
				int g = ClipChan(128 + (int)(deltaScale * (pred2[1] - pred1[1])));
				int b = ClipChan(128 + (int)(deltaScale * (pred2[2] - pred1[2])));
				image.setPixel(xx + x, yy + y, gARGB(0xff, r, g, b));
				pi.advance();
			}
		}

		// Make the predicted next observation
		xx = 2 * m_pParamRanges[0];
		pi.reset();
		for(size_t y = 0; y < m_pParamRanges[1]; y++)
		{
			for(size_t x = 0; x < m_pParamRanges[0]; x++)
			{
				pi.currentNormalized(buf2);
				GVec::copy(buf2 + m_paramDims, pActionsCat->row(i), pActionsCat->cols());
				pH->predict(buf2, pred2);
				int r = ClipChan((int)pred2[0]);
				int g = ClipChan((int)pred2[1]);
				int b = ClipChan((int)pred2[2]);
				image.setPixel(xx + x, yy + y, gARGB(0xff, r, g, b));
				pi.advance();
			}
		}
	}
image.savePng("mask.png");








/*
	GDynamicSystemNeighborFinder nf(pReducedObs, pActions, false, neighbors, m_pRand);

	//GNeighborFinderCacheWrapper nf2(&nf, false);
	//nf2.fillCache();
	//nf2.cutShortcuts(8);

	GManifoldLearner* pML;
	if(m_useIsomap)
	{
		pML = new GIsomap(neighbors, m_contextDims, m_pRand);
		((GIsomap*)pML)->setNeighborFinder(&nf);
	}
	else
	{
		pML = new GBreadthFirstUnfolding(1, // reps
			neighbors, m_contextDims, m_pRand);
		((GBreadthFirstUnfolding*)pML)->setNeighborFinder(&nf);
	}
	Holder<GManifoldLearner> hML(pML);
	GData* pEstState = pML->transform(pObservations);

	// Center the first state estimate at the origin
	double* pRow = pEstState->row(0);
	for(size_t i = 1; i < pEstState->rows(); i++)
		GVec::subtract(pEstState->row(i), pRow, m_contextDims);
	GVec::setAll(pRow, 0.0, m_contextDims);

	return pEstState;
*/
	return NULL;
}

void GRecurrentModel::trainObservationFunction(GData* pEstState, GData* pObservations)
{
	GData obs(m_paramDims + m_contextDims + m_channels);
	for(size_t i = 0; i < pObservations->rows(); i++)
	{
		GCoordVectorIterator pi(m_paramDims, m_pParamRanges);
		double* pObsIn = pObservations->row(i);
		while(true)
		{
			//if(m_pRand->uint(10) == 0) // sub-sample the pixels
			{
				double* pVec = obs.newRow();
				pi.currentNormalized(pVec);
				GVec::copy(pVec + m_paramDims, pEstState->row(i), m_contextDims);
				GVec::copy(pVec + m_paramDims + m_contextDims, pObsIn, m_channels);
			}
			pObsIn += m_channels;
			if(!pi.advance())
				break;
		}
	}
	m_pObservationFunc->train(&obs, m_channels);
}

void GRecurrentModel::trainObservationFunctionIteratively(double dStart, GData* pEstState, GData* pObservations)
{
	// Enable incremental learning
	if(!m_pObservationFunc->canTrainIncrementally())
		ThrowError("This is not an incremental learner");
	GIncrementalLearner* pObsLearner = (GIncrementalLearner*)m_pObservationFunc;
	{
		sp_relation pRel = new GUniformRelation(m_paramDims + m_contextDims + m_channels, 0);
		double* pMins = new double[2 * (m_paramDims + m_contextDims + m_channels)];
		ArrayHolder<double> hMins(pMins);
		double* pRanges = pMins + m_paramDims + m_contextDims + m_channels;
		double* pM = pMins;
		double* pR = pRanges;
		for(size_t i = 0; i < m_paramDims; i++)
		{
			*(pM++) = 0.0;
			*(pR++) = 1.0;
		}
		for(int i = 0; i < m_contextDims; i++)
			pEstState->minAndRangeUnbiased(i, pM++, pR++);
		for(int i = 0; i < m_channels; i++)
		{
			double m, r;
			pObservations->minAndRangeUnbiased(0, pM, pR);
			for(int j = 1; j < m_pixels; j++)
			{
				pObservations->minAndRangeUnbiased(m_channels * j, &m, &r);
				if(m < *pM)
					*pM = m;
				if(*pM + *pR < m + r)
					*pR = m + r - *pM;
			}
			pM++;
			pR++;
		}
		pObsLearner->enableIncrementalLearning(pRel, m_channels, pMins, pRanges);
	}

	// Do training
	int samples = 30;
	int timeSlice = -1;
	GCoordVectorIterator it(m_paramDims, m_pParamRanges);
	while(true)
	{
		double dTime = GTime::seconds() - dStart;
		if(m_pValidationData)
		{
			int slice = (int)floor(dTime / m_validationInterval);
			if(slice > timeSlice)
			{
				while(timeSlice + 1 < slice)
					onObtainValidationScore(++timeSlice, UNKNOWN_REAL_VALUE, UNKNOWN_REAL_VALUE);
				timeSlice = slice;
				double mse = validate(*m_pValidationData, false, false, 1.0);
				onObtainValidationScore(timeSlice, dTime, mse);
			}
		}
		if(dTime >= m_trainingSeconds)
			break;

		size_t row = (size_t)m_pRand->next(pEstState->rows());
		GVec::copy(m_pContext, pEstState->row(row), m_contextDims);
		double* pObs = pObservations->row(row);
		for(int i = 0; i < samples; i++)
		{
			it.setRandom(m_pRand);
			it.currentNormalized(m_pParams);
			size_t ofs = it.currentIndex();
			pObsLearner->trainIncremental(m_pParams, pObs + m_channels * ofs);
		}
	}
}

void GRecurrentModel::trainMoses(GData* pActions, GData* pObservations)
{
	// Consistency checks
	if(pActions->rows() != pObservations->rows())
		ThrowError("Expected the same number of rows");
	if(pActions->cols() != m_actionDims)
		ThrowError("Expected ", gformat(m_actionDims), " action dims, got ", gformat(pActions->cols()));
	if(pObservations->cols() != m_obsDims)
		ThrowError("Expected ", gformat(m_obsDims), " action dims, got ", gformat(pObservations->cols()));

	// Estimate state
	double dStart = GTime::seconds();
	GData* pEstState = mosesEstimateState(pActions, pObservations);
	Holder<GData> hEstState(pEstState);
	onFinishedComputingStateEstimate(pEstState);

	// Train the transition function
	trainTransitionFunction(pActions, pEstState);

	// Train the observation function
	if(m_pObservationFunc->canTrainIncrementally())
		trainObservationFunctionIteratively(dStart, pEstState, pObservations);
	else
		trainObservationFunction(pEstState, pObservations);
}

class GRecurrentModelContextCalibrationTargetFunction : public GTargetFunction
{
protected:
	GRecurrentModel* m_pModel;
	const double* m_pTarget;
	size_t m_stepSize;

public:
	GRecurrentModelContextCalibrationTargetFunction(GRecurrentModel* pModel, const double* pTarget, size_t stepSize = 1)
	: GTargetFunction(pModel->contextDims()), m_pModel(pModel), m_pTarget(pTarget), m_stepSize(stepSize)
	{
	}

	virtual ~GRecurrentModelContextCalibrationTargetFunction()
	{
	}

	virtual bool isStable() { return true; }
	virtual bool isConstrained() { return false; }

	virtual void initVector(double* pVector)
	{
		GVec::copy(pVector, m_pModel->context(), m_pModel->contextDims());
	}

	virtual double computeError(const double* pVector)
	{
		GVec::copy(m_pModel->context(), pVector, m_pModel->contextDims());
		int channels = m_pModel->observationFunc()->labelDims();
		GTEMPBUF(double, prediction, channels);
		double err = 0;
		const double* pTar = m_pTarget;
		GCoordVectorIterator pi(m_pModel->paramDims(), m_pModel->paramRanges());
		while(true)
		{
			pi.currentNormalized(m_pModel->params());
			m_pModel->observationFunc()->predict(m_pModel->params(), prediction);
			err += GVec::squaredDistance(pTar, prediction, channels);

			for(size_t i = 0; i < m_stepSize; i++)
			{
				if(!pi.advance())
					return err;
				pTar++;
			}
		}
	}
};

void GRecurrentModel::trainAaron(GData* pActions, GData* pObservations)
{
	// Consistency checks
	if(pActions->rows() != pObservations->rows())
		ThrowError("Expected the same number of rows");
	if(pActions->cols() != m_actionDims)
		ThrowError("Expected ", gformat(m_actionDims), " action dims, got ", gformat(pActions->cols()));
	if(pObservations->cols() != m_obsDims)
		ThrowError("Expected ", gformat(m_obsDims), " action dims, got ", gformat(pObservations->cols()));

	// Init the context
	GData contextData(m_contextDims);
	contextData.newRows(pActions->rows());
	for(size_t i = 0; i < contextData.rows(); i++)
	{
		double* pRow = contextData.row(i);
		m_pRand->cubical(pRow, m_contextDims);
	}

	// Enable incremental learning
	m_transitionDelta = false;
	GNeuralNet* pTransFunc = (GNeuralNet*)m_pTransitionFunc;
	GNeuralNet* pObsFunc = (GNeuralNet*)m_pObservationFunc;
	prepareForOptimization(pActions, pObservations);

	// Do bias training
	GCoordVectorIterator cvi(m_paramDims, m_pParamRanges);
	size_t samples = 30;
cerr << "starting bias training...\n"; cerr.flush();
	for(size_t i = 0; i < 500000; i++)
	{
		size_t t = (size_t)m_pRand->next(pActions->rows());
		double* pObs = pObservations->row(t);
		GVec::copy(m_pContext, contextData.row(t), m_contextDims);
		for(size_t j = 0; j < samples; j++)
		{
			cvi.setRandom(m_pRand);
			cvi.currentNormalized(m_pParams);
			pObsFunc->trainIncremental(m_pParams, pObs + m_channels * cvi.currentIndex());
		}
	}
cerr << "done with bias training...\n"; cerr.flush();

	// Train
	double dStart = GTime::seconds();
	int timeSlice = -1;
	GTEMPBUF(double, pPix, m_channels + m_contextDims);
	double* pTransContext = pPix + m_channels;
	GVec::setAll(contextData.row(0), 0.0, m_contextDims);
	while(true)
	{
		// Validate
		double dTime = GTime::seconds() - dStart;
		if(m_pValidationData)
		{
			int slice = (int)floor(dTime / m_validationInterval);
			if(slice > timeSlice)
			{
				while(timeSlice + 1 < slice)
					onObtainValidationScore(++timeSlice, UNKNOWN_REAL_VALUE, UNKNOWN_REAL_VALUE);
				timeSlice = slice;
				double mse = validate(*m_pValidationData, false, false, 255.0);
				onObtainValidationScore(timeSlice, dTime, mse);
			}
		}
		if(dTime >= m_trainingSeconds)
			break;

		// Do some training
		for(size_t t = 0; t < pActions->rows(); t++)
		{
			double* pObs = pObservations->row(t);
			if(t > 0)
			{
				// Compute the context according to the transition function
				GVec::copy(m_pBuf, pActions->row(t - 1), m_actionDims);
				GVec::copy(m_pBuf + m_actionDims, contextData.row(t - 1), m_contextDims);
				pTransFunc->predict(m_pBuf, pTransContext);
			}
			GVec::copy(m_pContext, contextData.row(t), m_contextDims);
			for(size_t i = 0; i < samples; i++)
			{
				cvi.setRandom(m_pRand);
				cvi.currentNormalized(m_pParams);
				pTransFunc->trainIncremental(m_pParams, pObs + m_channels * cvi.currentIndex());
				if(t > 0)
				{
					// Update the context
					GBackProp* pBP = pObsFunc->backProp();
					GBackPropLayer& bpLayer = pBP->layer(0);
					GNeuralNetLayer& nnLayer = pObsFunc->getLayer(0);
					for(int j = 0; j < m_contextDims; j++)
					{
						for(size_t k = 0; k < nnLayer.m_neurons.size(); k++)
							m_pContext[j] += 0.1 * bpLayer.m_neurons[k].m_error * nnLayer.m_neurons[k].m_weights[1 + j];
					}
				}
			}
			if(t > 0)
			{
				// Average the context estimates
//				GVec::add(contextData.row(t), pTransContext, m_contextDims);
//				GVec::multiply(contextData.row(t), 0.5, m_contextDims);

				// Train the transition function
				GVec::copy(m_pBuf, pActions->row(t - 1), m_actionDims);
				GVec::copy(m_pBuf + m_actionDims, contextData.row(t - 1), m_contextDims);
				pTransFunc->trainIncremental(m_pBuf, contextData.row(t));
			}


/*
				// Search for a context that satisfies the observation function
				double bestErr = 1e300;
				for(size_t i = 0; i < 7; i++)
				{
					// Pick a starting point
					if(i == 0)
						GVec::copy(m_pContext, contextData.row(t), m_contextDims);
					else if(i == 1)
						GVec::copy(m_pContext, pTransContext, m_contextDims);
					else if(i == 2)
						GVec::copy(m_pContext, contextData.row(t - 1), m_contextDims);
					else
						GVec::copy(m_pContext, contextData.row(m_pRand->next(contextData.rows())), m_contextDims);

					// Refine the context
					for(size_t i = 0; i < 100; i++)
					{
						cvi.setRandom(m_pRand);
						cvi.currentNormalized(m_pParams);
						pObsFunc->singleRowFeatureGradient(m_pParams, pObs + m_channels * cvi.currentIndex(), m_pBuf, m_paramDims);
						GVec::addScaled(m_pContext, -0.01, m_pBuf, m_contextDims);
					}

					// Evaluate the context and keep the best one
					double err = 0;
					cvi.reset();
					while(true)
					{
						cvi.currentNormalized(m_pParams);
						double* pTar = pObs + m_channels * cvi.currentIndex();
						pObsFunc->predict(m_pParams, pPix);
						err += GVec::squaredDistance(pTar, pPix, m_channels);
						if(!cvi.advance(7))
							break;
					}
					if(err < bestErr)
					{
						bestErr = err;
						GVec::copy(contextData.row(t), m_pContext, m_contextDims);
					}
				}

				// Average best observation context with the transition context
//				GVec::add(contextData.row(t), pTransContext, m_contextDims);
//				GVec::multiply(contextData.row(t), 0.5, m_contextDims);

				// Train the transition function
				GVec::copy(m_pBuf, pActions->row(t - 1), m_actionDims);
				GVec::copy(m_pBuf + m_actionDims, contextData.row(t - 1), m_contextDims);
				pTransFunc->trainIncremental(m_pBuf, contextData.row(t));
			}

			// Train the observation function
			GVec::copy(m_pContext, contextData.row(t), m_contextDims);
			for(size_t i = 0; i < samples; i++)
			{
				cvi.setRandom(m_pRand);
				cvi.currentNormalized(m_pParams);
				double* pTar = pObs + m_channels * cvi.currentIndex();
				pObsFunc->trainIncremental(m_pParams, pTar);
			}
*/
		}
	}

contextData.saveArff("aaron_context.arff");
}

void GRecurrentModel::trainJoshua(GData* pActions, GData* pObservations)
{
	// Consistency checks
	if(pActions->rows() != pObservations->rows())
		ThrowError("Expected the same number of rows");
	if(pActions->cols() != m_actionDims)
		ThrowError("Expected ", gformat(m_actionDims), " action dims, got ", gformat(pActions->cols()));
	if(pObservations->cols() != m_obsDims)
		ThrowError("Expected ", gformat(m_obsDims), " action dims, got ", gformat(pObservations->cols()));

	// Estimate state
//	double dStart = GTime::seconds();
	GData* pEstState = joshuaEstimateState(pActions, pObservations);
	Holder<GData> hEstState(pEstState);
/*
	onFinishedComputingStateEstimate(pEstState);

	// Train the transition function
	trainTransitionFunction(pActions, pEstState);

	// Train the observation function
	if(m_pObservationFunc->canTrainIncrementally())
		trainObservationFunctionIteratively(dStart, pEstState, pObservations);
	else
		trainObservationFunction(pEstState, pObservations);
*/
}

class GRecurrentModelTargetFunction : public GTargetFunction
{
protected:
	GRecurrentModel* m_pModel;
	size_t m_paramDims;
	size_t* m_pParamRanges;
	GRand* m_pRand;
	GData* m_pObservations;
	GData* m_pActions;
	GNeuralNet* m_pNNTrans;
	GNeuralNet* m_pNNObs;
	size_t m_transWeightCount;
	size_t m_obsWeightCount;
	double* m_pParamArray;
	GCoordVectorIterator m_it;

public:
	GRecurrentModelTargetFunction(GRecurrentModel* pModel, size_t paramDims, size_t* paramRanges, GRand* pRand, GData* pObservations, GData* pActions)
	: GTargetFunction(countWeights(pModel)), m_pModel(pModel), m_paramDims(paramDims), m_pParamRanges(paramRanges), m_pRand(pRand), m_it(paramDims, paramRanges)
	{
		m_pObservations = new GData(pObservations->relation());
		m_pActions = new GData(pActions->relation());
		for(size_t i = 0; i < 40 && i < pActions->rows(); i++)
		{
			m_pObservations->takeRow(pObservations->row(i));
			m_pActions->takeRow(pActions->row(i));
		}
		m_pParamArray = new double[30 * m_paramDims];
		pickRandomParams();
	}

	virtual ~GRecurrentModelTargetFunction()
	{
		m_pObservations->releaseAllRows();
		m_pActions->releaseAllRows();
		delete(m_pObservations);
		delete(m_pActions);
		delete[] m_pParamArray;
	}

	void pickRandomParams()
	{
		double* p = m_pParamArray;
		for(size_t i = 0; i < 30; i++)
		{
			m_it.setRandom(m_pRand);
			m_it.currentNormalized(p);
			p += m_paramDims;
		}
	}

	size_t countWeights(GRecurrentModel* pModel)
	{
		GIncrementalLearner* pTrans = (GIncrementalLearner*)pModel->transitionFunc();
		while(true)
		{
			if(!pTrans->canTrainIncrementally())
				ThrowError("Expected a filter or a neural net");
			if(!pTrans->isFilter())
				break;
			pTrans = (GIncrementalLearner*)((GFilter*)pTrans)->modeler();
		}
		m_pNNTrans = (GNeuralNet*)pTrans;
		m_transWeightCount = m_pNNTrans->countWeights();

		GIncrementalLearner* pObs = (GIncrementalLearner*)pModel->observationFunc();
		while(true)
		{
			if(!pObs->canTrainIncrementally())
				ThrowError("Expected a filter or a neural net");
			if(!pObs->isFilter())
				break;
			pObs = (GIncrementalLearner*)((GFilter*)pObs)->modeler();
		}
		m_pNNObs = (GNeuralNet*)pObs;
		m_obsWeightCount = m_pNNObs->countWeights();
		return m_transWeightCount + m_obsWeightCount;
	}

	virtual bool isStable() { return false; }
	virtual bool isConstrained() { return false; }

	virtual void initVector(double* pVector)
	{
		m_pNNTrans->weights(pVector);
		m_pNNObs->weights(pVector + m_transWeightCount);
	}

	virtual double computeError(const double* pVector)
	{
		m_pNNTrans->setWeights(pVector);
		m_pNNObs->setWeights(pVector + m_transWeightCount);
		return m_pModel->quickValidate(m_pActions, m_pObservations, 30, m_pParamArray, true);
	}
};

void GRecurrentModel::prepareForOptimization(GData* pActions, GData* pObservations)
{
	// Enable incremental learning
	if(!m_pTransitionFunc->canTrainIncrementally())
		ThrowError("Expected an incremental learner");
	{
		GMixedRelation* pMixedRel = new GMixedRelation();
		sp_relation pRel = pMixedRel;
		pMixedRel->addAttrs(pActions->relation().get());
		pMixedRel->addAttrs(2 * m_contextDims, 0);
		double* pMins = new double[2 * (pActions->cols() + m_contextDims + m_contextDims)];
		ArrayHolder<double> hMins(pMins);
		double* pRanges = pMins + pActions->cols() + m_contextDims + m_contextDims;
		GVec::setAll(pMins, 0.0, pActions->cols() + m_contextDims + m_contextDims);
		GVec::setAll(pRanges, 1.0, pActions->cols() + m_contextDims + m_contextDims);
		((GIncrementalLearner*)m_pTransitionFunc)->enableIncrementalLearning(pRel, m_contextDims, pMins, pRanges);
	}

	if(!m_pObservationFunc->canTrainIncrementally())
		ThrowError("Expected an incremental learner");
	GIncrementalLearner* pObsLearner = (GIncrementalLearner*)m_pObservationFunc;
	{
		sp_relation pRel = new GUniformRelation(m_paramDims + m_contextDims + m_channels, 0);
		double* pMins = new double[2 * (m_paramDims + m_contextDims + m_channels)];
		ArrayHolder<double> hMins(pMins);
		double* pRanges = pMins + m_paramDims + m_contextDims + m_channels;
		double* pM = pMins;
		double* pR = pRanges;
		for(size_t i = 0; i < m_paramDims; i++)
		{
			*pM++ = 0.0;
			*pR++ = 1.0;
		}
		for(int i = 0; i < m_contextDims; i++)
		{
			*pM++ = 0.0;
			*pR++ = 1.0;
		}
		for(int i = 0; i < m_channels; i++)
		{
			double m, r;
			pObservations->minAndRangeUnbiased(0, pM, pR);
			for(int j = 1; j < m_pixels; j++)
			{
				pObservations->minAndRangeUnbiased(m_channels * j, &m, &r);
				if(m < *pM)
					*pM = m;
				if(*pM + *pR < m + r)
					*pR = m + r - *pM;
			}
			pM++;
			pR++;
		}
		pObsLearner->enableIncrementalLearning(pRel, m_channels, pMins, pRanges);
	}
}

void GRecurrentModel::trainEvolutionary(GData* pActions, GData* pObservations)
{
	// Consistency checks
	if(pActions->rows() != pObservations->rows())
		ThrowError("Expected the same number of rows");
	if(pActions->cols() != m_actionDims)
		ThrowError("Expected ", gformat(m_actionDims), " action dims, got ", gformat(pActions->cols()));
	if(pObservations->cols() != m_obsDims)
		ThrowError("Expected ", gformat(m_obsDims), " action dims, got ", gformat(pObservations->cols()));

	prepareForOptimization(pActions, pObservations);
	GRecurrentModelTargetFunction tar(this, m_paramDims, m_pParamRanges, m_pRand, pObservations, pActions);
	GEvolutionaryOptimizer op(&tar, 30, m_pRand, 0.9);
	double dStart = GTime::seconds();
	int timeSlice = -1;
	while(true)
	{
		// Validate
		double dTime = GTime::seconds() - dStart;
		if(m_pValidationData)
		{
			int slice = (int)floor(dTime / m_validationInterval);
			if(slice > timeSlice)
			{
				while(timeSlice + 1 < slice)
					onObtainValidationScore(++timeSlice, UNKNOWN_REAL_VALUE, UNKNOWN_REAL_VALUE);
				timeSlice = slice;
				double mse = validate(*m_pValidationData, false, false, 1.0);
				onObtainValidationScore(timeSlice, dTime, mse);
			}
		}
		if(dTime >= m_trainingSeconds)
			break;

		// Do some training
		for(int i = 0; i < 30; i++)
		{
			tar.pickRandomParams();
			op.iterate();
		}
	}
}

void GRecurrentModel::trainHillClimber(GData* pActions, GData* pObservations, double dev, double decay, double seconds, bool climb, bool anneal)
{
	// Consistency checks
	if(pActions->rows() != pObservations->rows())
		ThrowError("Expected the same number of rows");
	if(pActions->cols() != m_actionDims)
		ThrowError("Expected ", gformat(m_actionDims), " action dims, got ", gformat(pActions->cols()));
	if(pObservations->cols() != m_obsDims)
		ThrowError("Expected ", gformat(m_obsDims), " action dims, got ", gformat(pObservations->cols()));

	prepareForOptimization(pActions, pObservations);
	GRecurrentModelTargetFunction tar(this, m_paramDims, m_pParamRanges, m_pRand, pObservations, pActions);
	GHillClimber op(&tar);
	double dStart = GTime::seconds();
	int timeSlice = -1;
	while(true)
	{
		// Validate
		double dTime = GTime::seconds() - dStart;
		if(m_pValidationData)
		{
			int slice = (int)floor(dTime / m_validationInterval);
			if(slice > timeSlice)
			{
				while(timeSlice + 1 < slice)
					onObtainValidationScore(++timeSlice, UNKNOWN_REAL_VALUE, UNKNOWN_REAL_VALUE);
				timeSlice = slice;
				double mse = validate(*m_pValidationData, false, false, 1.0);
				onObtainValidationScore(timeSlice, dTime, mse);
			}
		}
		if(dTime >= m_trainingSeconds)
			break;

		// Do some training
		if(climb)
		{
			tar.pickRandomParams();
			op.iterate();
		}
		if(anneal)
		{
			for(int i = 0; i < 30; i++)
			{
				tar.pickRandomParams();
				op.anneal(dev * pow(decay, dTime / seconds), m_pRand);
			}
		}
	}
}

size_t GRecurrentModel::trainBackPropThroughTime(GData* pActions, GData* pObservations, int depth, int itersPerSeqLen)
{
	// Consistency checks
	if(pActions->rows() != pObservations->rows())
		ThrowError("Expected the same number of rows");
	if(pActions->cols() != m_actionDims)
		ThrowError("Expected ", gformat(m_actionDims), " action dims, got ", gformat(pActions->cols()));
	if(!pActions->relation()->areContinuous(0, pActions->cols()))
		ThrowError("Expected continuous actions");
	if(pObservations->cols() != m_obsDims)
		ThrowError("Expected ", gformat(m_obsDims), " action dims, got ", gformat(pObservations->cols()));

	m_transitionDelta = false;
	GNeuralNet* pTransFunc = (GNeuralNet*)m_pTransitionFunc;
	GNeuralNet* pObsFunc = (GNeuralNet*)m_pObservationFunc;
	prepareForOptimization(pActions, pObservations);

	// Make a copy of the transition net for each unfolding-through-time
	sp_relation pRel = new GUniformRelation(pTransFunc->featureDims() + pTransFunc->labelDims(), 0);
	vector<GNeuralNet*> transNets;
	VectorOfPointersHolder<GNeuralNet> hTransNets(transNets);
	for(int i = 0; i < depth; i++)
	{
		GNeuralNet* pNet = new GNeuralNet(m_pRand);
		pNet->copyStructure(pTransFunc);
		transNets.push_back(pNet);
	}

	size_t obsLayers = pObsFunc->layerCount();
	size_t transLayers = pTransFunc->layerCount();
	size_t transWeightCount = pTransFunc->countWeights();
	double* pTransWeights = new double[2 * transWeightCount];
	ArrayHolder<double> hTransWeights(pTransWeights);
	double* pTransWeightsAccum = pTransWeights + transWeightCount;
	double* pAct = m_pBuf;
	double* pContextIn = pAct + m_actionDims;
	double* pContextOut = pContextIn + m_contextDims;
	GData contexts(m_contextDims);
	contexts.newRows(depth + 1);
	GTEMPBUF(double, pPix, m_channels);
	int seqIters = 0;
	size_t seqLen = 1;
	double dStart = GTime::seconds();
	int timeSlice = -1;
	GCoordVectorIterator it(m_paramDims, m_pParamRanges);
	while(true)
	{
		// Validate
		double dTime = GTime::seconds() - dStart;
		if(m_pValidationData)
		{
			int slice = (int)floor(dTime / m_validationInterval);
			if(slice > timeSlice)
			{
				while(timeSlice + 1 < slice)
					onObtainValidationScore(++timeSlice, UNKNOWN_REAL_VALUE, UNKNOWN_REAL_VALUE);
				timeSlice = slice;
				double mse = validate(*m_pValidationData, false, false, 255.0);
				onObtainValidationScore(timeSlice, dTime, mse);
			}
		}
		if(dTime >= m_trainingSeconds)
			break;

		// See if it's time to increase the effective sequence length
		if(seqIters >= itersPerSeqLen)
		{
			seqIters = 0;
			seqLen++;
		}

		// Pick a random pixel
		it.setRandom(m_pRand);
		it.currentNormalized(m_pParams);
		size_t ofs = it.currentIndex() * m_channels;

		// Set the initial context
		GVec::setAll(contexts.row(0), 0.0, m_contextDims);

		// Do back-prop-through-time across the whole sequence
		for(size_t dest = 0; dest < pActions->rows() && dest < seqLen; dest++) // dest is the index of the target observation
		{
			// Get the target pixel
			double* pTar = pObservations->row(dest) + ofs;

			// Forward prop over the unfolded transition function to compute the dest context
			pTransFunc->weights(pTransWeights);
			int passDepth = MIN(depth, (int)dest);
			size_t src = dest - passDepth;
			for(int i = 0; i < passDepth; i++)
			{
				GNeuralNet* pNN = transNets[i];
				pNN->setWeights(pTransWeights);
				GVec::copy(pAct, pActions->row(src + i), m_actionDims);
				GVec::copy(pContextIn, contexts.row(i), m_contextDims);
				pNN->predict(pAct, pContextOut);
				GVec::copy(contexts.row(i + 1), pContextOut, m_contextDims);
			}

			// Forward prop over observation function to compute target observation
			GVec::copy(m_pContext, contexts.row(passDepth), m_contextDims);
			pObsFunc->predict(m_pParams, pPix);

			// Compute the error on the output nodes
			GBackPropLayer& bpOutputLayer = pObsFunc->backProp()->layer(obsLayers - 1);
			GNeuralNetLayer& nnOutputLayer = pObsFunc->getLayer(obsLayers - 1);
			for(int i = 0; i < m_channels; i++)
				bpOutputLayer.m_neurons[i].m_error = (pTar[i] - pPix[i]) * nnOutputLayer.m_pActivationFunction->derivativeInverse(pPix[i]);

			// Compute the context gradient
			GBackProp* pBPObs = pObsFunc->backProp();

			// Back prop over observation function
			pBPObs->backpropagate();
			pBPObs->descendGradient(m_pParams, pObsFunc->learningRate(), pObsFunc->momentum());

			// Repeatedly back prop over transition function back to the src
			for(int i = passDepth - 1; i >= 0; i--)
			{
				GBackProp* pBPTrans = transNets[i]->backProp();
				if(i == passDepth - 1)
				{
					GBackPropLayer& bpTo = pBPTrans->layer(transLayers - 1);
					GBackProp::backPropLayer(&pObsFunc->getLayer(0), &transNets[i]->getLayer(transLayers - 1), &pBPObs->layer(0), &bpTo, m_paramDims);
				}
				else
					GBackProp::backPropLayer(&transNets[i + 1]->getLayer(0), &transNets[i]->getLayer(transLayers - 1), &transNets[i + 1]->backProp()->layer(0), &pBPTrans->layer(transLayers - 1), m_actionDims);
				GVec::copy(pAct, pActions->row(src + i), m_actionDims);
				GVec::copy(pContextIn, contexts.row(i), m_contextDims);
				GBackProp* pBPTemp = transNets[i]->backProp();
				pBPTemp->backpropagate();
				pBPTemp->descendGradient(pAct, transNets[i]->learningRate(), transNets[i]->momentum());
			}

			// Average the weights over the transition functions
			if(passDepth > 0)
			{
				GVec::setAll(pTransWeightsAccum, 0.0, transWeightCount);
				for(int i = 0; i < passDepth; i++)
				{
					transNets[i]->weights(pTransWeights);
					GVec::add(pTransWeightsAccum, pTransWeights, transWeightCount);
				}
				GVec::multiply(pTransWeightsAccum, 1.0 / passDepth, transWeightCount);
				pTransFunc->setWeights(pTransWeightsAccum);
			}

			// Shift the initial context
			if(depth > 1 && passDepth >= depth)
				GVec::copy(contexts.row(0), contexts.row(1), m_contextDims);
		}
		seqIters++;
	}
//cout << "% seq len = " << seqLen << "\n";
	return seqLen;
}

void GRecurrentModel::doAction(const double* pAction)
{
	GVec::copy(m_pBuf, pAction, m_actionDims);
	double* pContext = m_pBuf + m_actionDims;
	GVec::copy(pContext, m_pContext, m_contextDims);
	if(m_transitionDelta)
	{
		double* pDelta = pContext + m_contextDims;
		m_pTransitionFunc->predict(m_pBuf, pDelta);
		GVec::add(m_pContext, pDelta, m_contextDims);
	}
	else
		m_pTransitionFunc->predict(m_pBuf, m_pContext);
}

void GRecurrentModel::predict(double* pObs)
{
	GCoordVectorIterator pi(m_paramDims, m_pParamRanges);
	while(true)
	{
		pi.currentNormalized(m_pParams);
		m_pObservationFunc->predict(m_pParams, pObs);
		pObs += m_pObservationFunc->labelDims();
		if(!pi.advance())
			break;
	}
}

void GRecurrentModel::predictPixel(const double* pParams, double* pObs)
{
	GVec::copy(m_pParams, pParams, m_paramDims);
	m_pObservationFunc->predict(m_pParams, pObs);
}

// virtual
void GRecurrentModel::calibrate(const double* pTarget)
{
/*
	// Assuming the observation function is a neural net, compute the gradient with respect to
	// the context, so we can calibrate the context
	GTEMPBUF(double, gradient, m_contextDims * 3);
	double* sumGradient = gradient + m_contextDims;
	double* prevSumGradient = sumGradient + m_contextDims;
	double rate = 0.1;
	double prevErr = 1e308;
	GVec::setAll(prevSumGradient, 0.0, m_contextDims);
	GBackPropLayer& layer = m_pObservationFunc->backProp()->layer(m_pObservationFunc->layerCount() - 1);
	GBackProp* pBackProp = m_pObservationFunc->backProp();
	for(int i = 0; i < 100; i++)
	{
		// Compute the gradient
		GVec::setAll(m_pBuf, 0.0, m_paramDims);
		GVec::setAll(sumGradient, 0.0, m_contextDims);
		GVec::copy(m_pBuf + m_paramDims, m_pContext, m_contextDims);
		double* pPrediction = m_pBuf + m_paramDims + m_contextDims;
		double* pTar = pTarget;
		size_t pixels = 0;
		double err = 0;
		while(true)
		{
			m_pObservationFunc->predict(m_pBuf, pPrediction);
			err += GVec::squaredDistance(pTar, pPrediction, layer.m_neurons.size());
			for(size_t i = 0; i < layer.m_neurons.size(); i++)
				layer.m_neurons[i].m_error = pTar[i] - pPrediction[i];
			pBackProp->computeInputGradient(gradient);
			GVec::add(sumGradient, gradient, m_contextDims);
			pTar += m_pObservationFunc->labelDims();
			pixels++;

			// Increment the parameters
			if(incParams(m_pBuf))
				break;
		}
		GAssert(pixels == m_pixels);
		GVec::multiply(sumGradient, 1.0 / (pixels * layer.m_neurons.size()), m_contextDims);

		if(err < prevErr)
		{
			// Step and accelerate
			GVec::addScaled(m_pContext, -rate, sumGradient, m_contextDims);
			prevErr = err;
			std::swap(sumGradient, prevSumGradient);
		}
		else
		{
			// Undo the previous acceleration, back up, and decelerate
			rate /= 1.15;
			GVec::addScaled(m_pContext, 0.75 * rate, prevSumGradient, m_contextDims);
			rate = MAX(rate * 0.25, 1e-12);
		}
		rate *= 1.15;
	}
*/
	// Use a hill climber to calibrate the context
	GRecurrentModelContextCalibrationTargetFunction tar(this, pTarget);
	//GAnnealing hc(&tar, 0.1, 1.0, m_pRand);
	GHillClimber hc(&tar);
	for(int i = 0; i < 30; i++)
		hc.iterate();
	GVec::copy(m_pContext, hc.currentVector(), m_contextDims);
}

double GRecurrentModel::validate(vector<GData*>& validationData, bool calibrateContext, bool monotonic, double multiplier)
{
	double mse = 0;
	double* pPrediction = new double[m_obsDims];
	ArrayHolder<double> hPrediction(pPrediction);
	if(validationData.size() & 1)
		ThrowError("Expected an even number of datasets, one action set after each observation set");
	size_t count = 0;
	m_multiplier = multiplier;
	for(size_t j = 0; j < validationData.size() / 2; j++)
	{
		GData* pDataObs = validationData[2 * j];
		if(pDataObs->cols() != m_obsDims)
			ThrowError("Wrong number of dims in the obs validation data at index ", gformat(2 * j));
		GData* pDataAction = validationData[2 * j + 1];
		if(pDataAction->cols() != m_actionDims)
			ThrowError("Wrong number of dims in the action validation data at index ", gformat(2 * j + 1));
		GVec::setAll(m_pContext, 0.0, m_contextDims);
		double err = 0.0;
		for(size_t i = 0; i + 1 < pDataAction->rows(); i++)
		{
			if(calibrateContext)
				calibrate(pDataObs->row(i));
			doAction(pDataAction->row(i));
			predict(pPrediction);
			double d = GVec::squaredDistance(pDataObs->row(i + 1), pPrediction, m_obsDims) / m_obsDims;
			count++;
			if(multiplier != 1.0)
			{
				d = sqrt(d) * multiplier;
				d *= d;
			}
			if(monotonic)
				err = MAX(err, d);
			else
				err = d;
			mse += err;
		}
	}
	return mse / count;
}

double GRecurrentModel::quickValidate(GData* pDataAction, GData* pDataObs, int pixelSamples, double* paramArray, bool monotonic)
{
	if(m_paramDims != 2)
		ThrowError("Sorry, this method is currently only implemented for 2 param dims");
	GVec::setAll(m_pContext, 0.0, m_contextDims);
	double sse = 0;
	GTEMPBUF(double, pObs, m_channels + pixelSamples);
	double* err = pObs + m_channels;
	GVec::setAll(err, 0.0, pixelSamples);
	for(size_t i = 0; i + 1 < pDataAction->rows(); i++)
	{
		doAction(pDataAction->row(i));

		double* pParams = paramArray;
		double* pTarget = pDataObs->row(i + 1);
		for(int j = 0; j < pixelSamples; j++)
		{
			GVec::copy(m_pParams, pParams, m_paramDims);
			m_pObservationFunc->predict(m_pParams, pObs);
			double d = GVec::squaredDistance(pTarget + m_channels * ((int)pParams[0] + (int)m_pParamRanges[0] * (int)pParams[1]), pObs, m_channels);
			pParams += m_paramDims;
			if(monotonic)
				err[j] = MAX(err[j], d);
			else
				err[j] = d;
			sse += err[j];
		}
	}
	return sse;
}

GImage* GRecurrentModel::frames(GData* pDataAction, GData* pDataObs, bool calibrateContext, unsigned int frameWidth, int stepsPerImage, double scalePredictions)
{
	// Allocate the image
	if(pDataObs && pDataAction->rows() != pDataObs->rows())
		ThrowError("Expected same number of rows");
	if(m_paramDims != 2)
		ThrowError("Expected 2 params (x and y)");
	if(frameWidth < m_pParamRanges[0])
		ThrowError("Expected frameWidth to be at least ", gformat(m_pParamRanges[0]));
	if(frameWidth % m_pParamRanges[0] != 0)
		ThrowError("Expected frameWidth to be a multiply of ", gformat(m_pParamRanges[0]));
	if(m_channels != 3)
		ThrowError("Expected 3 channels");
	if(m_pObservationFunc->labelDims() != m_channels)
		ThrowError("Something is wrong");
	if(calibrateContext && !pDataObs)
		ThrowError("Cannot calibrate if observations are not given");
	double prediction[3];
	if(pDataObs && (size_t)pDataObs->cols() != m_pParamRanges[0] * m_pParamRanges[1] * m_channels)
		ThrowError("observation cols don't match specified parameters");
	unsigned int frameHeight = frameWidth / m_pParamRanges[0] * m_pParamRanges[1];
	unsigned int imageHeight = frameHeight * ((pDataAction->rows() - 1) / stepsPerImage);
	GImage* pImage = new GImage();
	Holder<GImage> hImage(pImage);
	pImage->setSize(frameWidth * 2, imageHeight);

	// Walk the sequence
	m_multiplier = scalePredictions;
	int scale = frameWidth / m_pParamRanges[0];
	GVec::setAll(m_pContext, 0.0, m_contextDims);
	int top = 0;
	for(size_t i = 0; i + 1 < pDataAction->rows(); i++)
	{
		if(calibrateContext)
			calibrate(pDataObs->row(i));
		if((i % stepsPerImage) == 0 && (unsigned int)top < imageHeight)
		{
			// Make the target image
			if(pDataObs)
			{
				double* pTarget = pDataObs->row(i);
				for(unsigned int y = 0; y < frameHeight; y++)
				{
					int rowStart = m_pParamRanges[0] * m_channels * (y * m_pParamRanges[1] / frameHeight);
					unsigned int* pPix = pImage->pixelRef(0, top + y);
					for(unsigned int x = 0; x < frameWidth; x++)
					{
						int col = x * m_pParamRanges[0] / frameWidth;
						int r = ClipChan((int)pTarget[rowStart + m_channels * col]);
						int g = ClipChan((int)pTarget[rowStart + m_channels * col + 1]);
						int b = ClipChan((int)pTarget[rowStart + m_channels * col + 2]);
						*(pPix++) = gARGB(0xff, r, g, b);
					}
				}
			}

			// Make the predicted image
			size_t scaledParams[2];
			scaledParams[0] = scale * m_pParamRanges[0];
			scaledParams[1] = scale * m_pParamRanges[1];
			GCoordVectorIterator it(2, scaledParams);
			int x = 0;
			int y = 0;
			while(true)
			{
				it.currentNormalized(m_pParams);
				m_pObservationFunc->predict(m_pParams, prediction);
				int r = ClipChan((int)(prediction[0] * scalePredictions));
				int g = ClipChan((int)(prediction[1] * scalePredictions));
				int b = ClipChan((int)(prediction[2] * scalePredictions));
				unsigned int pix = gARGB(0xff, r, g, b);
				pImage->setPixel(frameWidth + x, top + y, pix);
				if(++x >= (int)frameWidth)
				{
					x = 0;
					y++;
				}
				if(!it.advance())
					break;
			}
			top += frameHeight;
		}
		doAction(pDataAction->row(i));
	}
	return hImage.release();
}

// ------------------------------------------------------------------------------------------
/*

GManifoldDynamicsLearner::GManifoldDynamicsLearner(sp_relation& pRelation, int actionDims, int shortTermMemorySize, int contextDims, double minCorrelation, GAgentActionIterator* pActionIterator, GRand* pRand)
: m_pRelation(pRelation), m_senseDims(pRelation->size() - actionDims), m_actionDims(actionDims),
	m_minCorrelation(minCorrelation),
	m_contextDims(contextDims),
	m_longTermMemory((1 + contextDims) * m_senseDims + pActionIterator->actionCount() * contextDims),
	m_shortTermMemory(pRelation),
	m_shortTermMemoryPos(0),
	m_shortTermMemoryCount(0),
	m_pActionIterator(pActionIterator),
	m_pRand(pRand)
{
	m_neighborCount = 14;
	m_shortTermMemory.reserve(shortTermMemorySize);
	int actionCount = m_pActionIterator->actionCount();
	if(actionCount > 200)
		ThrowError("Sorry, only small discrete numbers of actions are supported by this algorithm");
	if(m_contextDims > m_senseDims)
		ThrowError("Sorry, context dims must be smaller than sense dims");
	m_shortTermMemory.newRows(shortTermMemorySize);
	double* pCurPat = m_shortTermMemory.row(m_shortTermMemoryPos);
	GVec::setAll(pCurPat, UNKNOWN_REAL_VALUE, m_senseDims + m_actionDims);
	m_actionVecStart = m_senseDims + contextDims * m_senseDims;
	m_patchSize = m_actionVecStart + actionCount * contextDims;
	m_pNeighborRelation = new GUniformRelation(m_senseDims);
	GDissimilarityMetric* pMetric = new GRowDistance();
	m_pNeighborFinder = new GKdTree(&m_longTermMemory, m_longTermMemory.cols() - m_senseDims, m_neighborCount, pMetric, true);
//	m_pNeighborFinder = new GBruteForceNeighborFinder(&m_longTermMemory, m_neighborCount, pMetric, true, NULL);
	m_currentPatch = -1;
	m_pContext = new double[m_contextDims + m_senseDims + MAX(m_senseDims, m_actionDims)];
	GVec::setAll(m_pContext, 0.0, m_contextDims);
	m_pBuf = m_pContext + m_contextDims;
}

// virtual
GManifoldDynamicsLearner::~GManifoldDynamicsLearner()
{
	delete[] m_pContext;
	delete(m_pNeighborFinder);
}

#ifdef _DEBUG
int g_logImageId = 0;
void GManifoldDynamicsLearner::logImage(int patch, double a, double b, const char* szMessage)
{
	double* pPatch = m_longTermMemory.row(patch);

	// Copy the mean
	GVec::copy(m_pBuf + m_senseDims, pPatch, m_senseDims);

	// Add in the context
	GVec::addScaled(m_pBuf + m_senseDims, a, pPatch + m_senseDims, m_senseDims);
	GVec::addScaled(m_pBuf + m_senseDims, b, pPatch + 2 * m_senseDims, m_senseDims);

	// Make the image
	GImage tmp;
	tmp.setSize(12, 12);
	double* pVec = m_pBuf + m_senseDims;
	for(int y = 0; y < 12; y++)
	{
		for(int x = 0; x < 12; x++)
		{
			int r = ClipChan(((int)*pVec * 256));
			pVec++;
			int g = ClipChan(((int)*pVec * 256));
			pVec++;
			int b = ClipChan(((int)*pVec * 256));
			pVec++;
			tmp.setPixel(x, y, gARGB(0xff, r, g, b));
		}
	}

	// Save the image
	string filename = "log";
	filename += gformat(g_logImageId++);
	filename += ".png";
	tmp.savePng(filename.c_str());

	// Print the message
	//cout << filename << " = " << szMessage;
}
#endif // _DEBUG

double GManifoldDynamicsLearner::computeDihedralCorrelation(size_t a, size_t b)
{
	if(a == b)
		return 1.0;
	double* pPatchA = m_longTermMemory.row(a);
	double* pPatchB = m_longTermMemory.row(b);
	m_pRand->spherical(m_pBuf, m_senseDims);
	GVec::copy(m_pBuf + m_senseDims, m_pBuf, m_senseDims);
	for(int i = 0; i < m_contextDims; i++)
	{
		GVec::subtractComponent(m_pBuf, pPatchA + (i + 1) * m_senseDims, m_senseDims);
		GVec::subtractComponent(m_pBuf + m_senseDims, pPatchB + (i + 1) * m_senseDims, m_senseDims);
	}
	return ABS(GVec::correlation(m_pBuf, m_pBuf + m_senseDims, m_senseDims));
}

void GManifoldDynamicsLearner::makeNewPatch()
{
	m_shortTermMemoryCount = MAX(0, m_shortTermMemoryCount - 8); // how often to make a new patch

	// Make a new patch
	double* pPatch = new double[m_patchSize];

	// Compute the mean of sense data in short term memory
	for(int i = 0; i < m_senseDims; i++)
		pPatch[i] = m_shortTermMemory.mean(i);

	// Compute the basis vectors of the sense data in short term memory
	GHeap heapTmp(2048);
	GData dataTmp(m_shortTermMemory.relation(), &heapTmp);
	dataTmp.copy(&m_shortTermMemory);
	for(int i = 0; i < m_contextDims; i++)
	{
		dataTmp.principalComponent(pPatch + (i + 1) * m_senseDims, m_senseDims, pPatch, m_pRand);
		dataTmp.removeComponent(pPatch, pPatch + (i + 1) * m_senseDims, m_senseDims);
	}

	// Compute the action vectors
	m_pActionIterator->reset(NULL);
	for(int action = 0; true; action++)
	{
		if(!m_pActionIterator->nextAction(m_pBuf + m_senseDims))
			break;

		// Compute the average observation delta for all samples in short-term memory due to this action
		GVec::setAll(m_pBuf, 0.0, m_senseDims);
		int count = 0;
		for(size_t j = 0; j < m_shortTermMemory.rows() - 1; j++)
		{
			double* pPat = m_shortTermMemory.row((m_shortTermMemoryPos + j) % m_shortTermMemory.rows());
			GAssert(*pPat != UNKNOWN_REAL_VALUE);
			if(GVec::squaredDistance(m_pBuf + m_senseDims, pPat + m_senseDims, m_actionDims) < 0.5)
			{
				double* pPatNext = m_shortTermMemory.row((m_shortTermMemoryPos + j + 1) % m_shortTermMemory.rows());
				GVec::add(m_pBuf, pPatNext, m_senseDims);
				GVec::subtract(m_pBuf, pPat, m_senseDims);
				count++;
			}
		}
		if(count > 0)
			GVec::multiply(m_pBuf, 1.0 / count, m_senseDims);
		else
			cout << "uh oh, action " << action << " has no representation in patch " << m_longTermMemory.rows() << "\n";

		// Dot-product deltas with basis vectors
		double* pActionVec = pPatch + m_actionVecStart + action * m_contextDims;
		for(int dim = 0; dim < m_contextDims; dim++)
		{
			pActionVec[dim] = GVec::dotProduct(m_pBuf, pPatch + (dim + 1) * m_senseDims, m_senseDims);
			GAssert(pActionVec[dim] > -1e20 && pActionVec[dim] < 1e20); // unreasonable value
		}
	}

	// todo: merge patches if they are very close?

	// Add the patch to the neighbor finder
	m_pNeighborFinder->addVector(pPatch);
//LogImage(m_longTermMemory.rows() - 1, 0, 0, "New patch");
}

void GManifoldDynamicsLearner::recomputeContext(const double* pSenses)
{
	// Find the best patch
	GTEMPBUF(size_t, neighbors, m_neighborCount);
	GTEMPBUF(double, squaredDist, m_neighborCount);
	m_pNeighborFinder->neighbors(neighbors, squaredDist, pSenses);
	m_pNeighborFinder->sortNeighbors(neighbors, squaredDist);
	for(int i = 0; i < m_neighborCount; i++)
	{
		if(neighbors[i] < 0)
			break;
		double corr = computeDihedralCorrelation(m_currentPatch, neighbors[i]);
		if(corr >= m_minCorrelation)
		{
			m_currentPatch = neighbors[i];
			break;
		}
	}

	// Recompute context
	double* pPatch = m_longTermMemory.row(m_currentPatch);
	for(int dim = 0; dim < m_contextDims; dim++)
		m_pContext[dim] = GVec::dotProduct(pPatch, pSenses, pPatch + (dim + 1) * m_senseDims, m_senseDims);
}

// virtual
void GManifoldDynamicsLearner::doAction(const double* pActions)
{
	// Learn
	if(m_contextStack.empty())
	{
		double* pCurPat = m_shortTermMemory.row(m_shortTermMemoryPos);
		if(++m_shortTermMemoryPos >= (int)m_shortTermMemory.rows())
			m_shortTermMemoryPos = 0;
		if(pCurPat[0] == UNKNOWN_REAL_VALUE)
		{
			// Didn't store actual observations, so we shouldn't do any learning
			m_shortTermMemoryCount = 0;
		}
		else
		{
			// Store the action
			GVec::copy(pCurPat + m_senseDims, pActions, m_actionDims);
			m_shortTermMemoryCount = MIN((int)m_shortTermMemory.rows(), m_shortTermMemoryCount + 1);

			// Create a new patch
			if(m_shortTermMemoryCount >= (int)m_shortTermMemory.rows())
				makeNewPatch();

			// Clear the new memory slot
			pCurPat = m_shortTermMemory.row(m_shortTermMemoryPos);
			GVec::setAll(pCurPat, UNKNOWN_REAL_VALUE, m_senseDims + m_actionDims);
		}
	}

	// Simulate the action
	if(m_currentPatch < m_longTermMemory.rows())
	{
		// Determine the action index
		m_pActionIterator->reset(NULL);
		int action;
		for(action = 0; true; action++)
		{
			if(!m_pActionIterator->nextAction(m_pBuf))
				ThrowError("action not found");
			if(GVec::squaredDistance(m_pBuf, pActions, m_actionDims) < 0.5)
				break;
		}

		// Adjust the context
		double* pPatch = m_longTermMemory.row(m_currentPatch);
		double* pActionVec = pPatch + m_actionVecStart + action * m_contextDims;
		for(int dim = 0; dim < m_contextDims; dim++)
			m_pContext[dim] += pActionVec[dim];

		// If we're just predicting, transition to a better patch
		if(m_shortTermMemoryCount == 0)
		{
			predict(m_pBuf);
			recomputeContext(m_pBuf);
		}
	}
}

// virtual
void GManifoldDynamicsLearner::predict(double* pSenses)
{
	// Find the current patch
	if(m_currentPatch >= m_longTermMemory.rows())
	{
		if(m_longTermMemory.rows() > 0)
			m_currentPatch = 0;
		else
		{
			GVec::setAll(pSenses, 0.0, m_senseDims);
			return;
		}
	}
	double* pPatch = m_longTermMemory.row(m_currentPatch);

	// Copy the mean
	GVec::copy(pSenses, pPatch, m_senseDims);

	// Add in the context
	for(int dim = 0; dim < m_contextDims; dim++)
	{
		double* pBasis = pPatch + (dim + 1) * m_senseDims;
		GVec::addScaled(pSenses, m_pContext[dim], pBasis, m_senseDims);
	}
}

// virtual
void GManifoldDynamicsLearner::calibrate(const double* pSenses)
{
	if(!m_contextStack.empty())
		ThrowError("You're not supposed to call this method when the context stack is not empty");

	// Remember the observations
	double* pCurPat = m_shortTermMemory.row(m_shortTermMemoryPos);
	GVec::copy(pCurPat, pSenses, m_senseDims);

	// Calibrate the context
	if(m_currentPatch < m_longTermMemory.rows())
		recomputeContext(pSenses);
}

void GManifoldDynamicsLearner::pushContext()
{
	m_contextStack.push_back((double)m_currentPatch);
	for(int i = 0; i < m_contextDims; i++)
		m_contextStack.push_back(m_pContext[i]);
}

void GManifoldDynamicsLearner::popContext()
{
	for(int i = m_contextDims - 1; i >= 0; i--)
	{
		m_pContext[i] = *(m_contextStack.end() - 1);
		m_contextStack.pop_back();
	}
	m_currentPatch = (size_t)*(m_contextStack.end() - 1);
	m_contextStack.pop_back();
}


// ---------------------------------------------------------------------


GTemporalInstanceLearner::GTemporalInstanceLearner(sp_relation& pRelation, int actionDims, int shortTermMemorySize, GRand* pRand)
: m_pRelation(pRelation), m_senseDims(pRelation->size() - actionDims), m_actionDims(actionDims),
	m_shortTermMemory(pRelation),
	m_shortTermMemoryPos(0),
	m_shortTermMemorySize(0),
	m_pRand(pRand)
{
	m_actionCount = pRelation->valueCount(pRelation->size() - 1);
	if(actionDims != 1 || m_actionCount == 0)
		ThrowError("Expected exactly 1 nominal action dimension");
	m_balance = 0.5;
	m_pActionIterator = new GDiscreteActionIterator(m_actionCount);
	m_shortTermMemory.newRows(shortTermMemorySize);
	m_longTermMemories = new GData*[m_actionCount];
	for(int i = 0; i < m_actionCount; i++)
		m_longTermMemories[i] = new GData((m_actionCount + 2) * m_senseDims);
	m_pActionEffects = new double[m_actionCount * m_senseDims];
}

// virtual
GTemporalInstanceLearner::~GTemporalInstanceLearner()
{
	delete(m_pActionIterator);
	for(int i = 0; i < m_actionCount; i++)
		delete(m_longTermMemories[i]);
	delete[] m_longTermMemories;
	delete[] m_pActionEffects;
}


bool GTemporalInstanceLearner::computeAverageActionEffects(double* pVec, int offset)
{
	// Compute the average action effects
	bool full = true;
	for(int j = 0; j < m_actionCount; j++)
	{
		GVec::setAll(pVec, 0.0, m_senseDims);
		int count = 0;
		for(size_t i = 0; i < m_shortTermMemory.rows() - 2; i++)
		{
			double* pPat = m_shortTermMemory.row((m_shortTermMemoryPos + 1 + offset + i) % m_shortTermMemory.rows());
			if((int)pPat[m_senseDims] == j)
			{
				double* pPat2 = m_shortTermMemory.row((m_shortTermMemoryPos + 2 + offset + i) % m_shortTermMemory.rows());
				GVec::add(pVec, pPat2, m_senseDims);
				GVec::subtract(pVec, pPat, m_senseDims);
				count++;
			}
		}
		if(count > 0)
			GVec::multiply(pVec, 1.0 / count, m_senseDims);
		else
			full = false;
		pVec += m_senseDims;
	}
	return full;
}

void GTemporalInstanceLearner::makeNewInstance()
{
	int prevPos = (int)((m_shortTermMemoryPos + m_shortTermMemory.rows() - 1) % m_shortTermMemory.rows());
	double* pPrevPat = m_shortTermMemory.row(prevPos);
	double* pCurPat = m_shortTermMemory.row(m_shortTermMemoryPos);
	int action = (int)pPrevPat[m_senseDims];
	GAssert(action >= 0 && action < m_actionCount);

	// Compute the average action effects
	double* pNewInstance = m_longTermMemories[action]->newRow();
	if(!computeAverageActionEffects(pNewInstance, 0))
	{
		cout << "not all actions are represented--nixing instance\n";
		m_longTermMemories[action]->deleteRow(m_longTermMemories[action]->rows() - 1);
		return;
	}

	// Store the position and delta
	double* pVec = pNewInstance + m_actionCount * m_senseDims;
	GVec::copy(pVec, pPrevPat, m_senseDims);
	pVec += m_senseDims;
	GVec::copy(pVec, pCurPat, m_senseDims);
	GVec::subtract(pVec, pPrevPat, m_senseDims);
}

size_t GTemporalInstanceLearner::findBestMatch(int action)
{
	computeAverageActionEffects(m_pActionEffects, 1);
	double* pCurPat = m_shortTermMemory.row(m_shortTermMemoryPos);

	size_t bestIndex = m_longTermMemories[action]->rows();
	double bestErr = 1e308;
	GAssert(action >= 0 && action < m_actionCount);
	for(size_t i = 0; i < m_longTermMemories[action]->rows(); i++)
	{
		double* pCand = m_longTermMemories[action]->row(i);

		// Compute the error
		double err = 0;
		err += (1.0 - m_balance) * GVec::squaredDistance(m_pActionEffects, pCand, m_actionCount * m_senseDims) / m_actionCount;
		err += m_balance * GVec::squaredDistance(pCurPat, pCand + m_actionCount * m_senseDims, m_senseDims);

		if(err < bestErr)
		{
			bestIndex = i;
			bestErr = err;
		}
	}
	return bestIndex;
}

// virtual
void GTemporalInstanceLearner::doAction(const double* pActions)
{
	// Store the action
	int nextPos = (int)((m_shortTermMemoryPos + 1) % m_shortTermMemory.rows());
	double* pCurPat = m_shortTermMemory.row(m_shortTermMemoryPos);
	double* pNextPat = m_shortTermMemory.row(nextPos);
	pCurPat[m_senseDims] = pActions[0];

	// Advance prediction
	int action = (int)pActions[0];
	GVec::copy(pNextPat, pCurPat, m_senseDims);
	if(m_shortTermMemorySize >= m_shortTermMemory.rows())
	{
		size_t match = findBestMatch(action);
		if(match < m_longTermMemories[action]->rows())
		{
			double* pMatch = m_longTermMemories[action]->row(match);
			GVec::add(pNextPat, pMatch + (m_actionCount + 1) * m_senseDims, m_senseDims);
		}
	}

	// Advance the memory
	m_shortTermMemoryPos = nextPos;
	m_shortTermMemorySize = MIN(m_shortTermMemorySize + 1, m_shortTermMemory.rows());
}

// virtual
void GTemporalInstanceLearner::predict(double* pSenses)
{
	double* pCurPat = m_shortTermMemory.row(m_shortTermMemoryPos);
	GVec::copy(pSenses, pCurPat, m_senseDims);
}

// virtual
void GTemporalInstanceLearner::calibrate(const double* pSenses)
{
	// Remember the observations
	double* pCurPat = m_shortTermMemory.row(m_shortTermMemoryPos);
	GVec::copy(pCurPat, pSenses, m_senseDims);

	// Calibrate the context
	if(m_shortTermMemorySize + 1 >= m_shortTermMemory.rows())
		makeNewInstance();
}
*/
} // namespace GClasses

