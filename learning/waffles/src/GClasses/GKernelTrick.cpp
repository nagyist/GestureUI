/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#include "GKernelTrick.h"
#include "GHillClimber.h"
#include "GDistribution.h"
#include "GMath.h"

using namespace GClasses;

GKernelMachine::GKernelMachine()
: GSupervisedLearner(), m_featureDims(0), m_labelDims(0), m_pWeights(NULL), m_pSupportVectors(NULL)
{
	m_pKernel = NULL;
	m_pBuf = NULL;
}

GKernelMachine::GKernelMachine(GTwtNode* pLearner)
: GSupervisedLearner(pLearner)
{
	ThrowError("Sorry, not implemented yet");
}

// virtual
GKernelMachine::~GKernelMachine()
{
	delete(m_pKernel);
	delete[] m_pWeights;
	delete(m_pSupportVectors);
	delete[] m_pBuf;
}

// virtual
int GKernelMachine::featureDims()
{
	if(m_labelDims < 1)
		ThrowError("not yet trained");
	return m_featureDims;
}

// virtual
int GKernelMachine::labelDims()
{
	if(m_labelDims < 1)
		ThrowError("not yet trained");
	return m_labelDims;
}

double* GKernelMachine::swapWeights(double* pNewWeights)
{
	double* pOld = m_pWeights;
	m_pWeights = pNewWeights;
	return pOld;
}

// static
GKernel* GKernelMachine::makeKernel(int dims)
{
	//return new GKernelIdentity(dims);
	//return new GKernelPolynomial(dims, 1, 7);

	GKernel* pK1 = new GKernelPolynomial(dims, 0, 3);
	GKernel* pK2 = new GKernelPolynomial(dims, 1, 7);
	GKernel* pK3 = new GKernelAdd(pK1, pK2);
	GKernel* pK4 = new GKernelNormalize(pK3);

	GKernel* pK5 = new GKernelGaussianRBF(dims, 0.01);
	GKernel* pK6 = new GKernelGaussianRBF(dims, 0.1);
	GKernel* pK7 = new GKernelAdd(pK5, pK6);
	GKernel* pK8 = new GKernelNormalize(pK7);

	GKernel* pK9 = new GKernelGaussianRBF(dims, 1.0);
	GKernel* pK10 = new GKernelGaussianRBF(dims, 10.0);
	GKernel* pK11 = new GKernelMultiply(pK9, pK10);
	GKernel* pK12 = new GKernelNormalize(pK11);

	GKernel* pK13 = new GKernelAdd(pK8, pK12);
	GKernel* pK14 = new GKernelAdd(pK4, pK13);
	return pK14;
}

// virtual
GTwtNode* GKernelMachine::toTwt(GTwtDoc* pDoc)
{
	ThrowError("Sorry, not implemented yet");
	return NULL;
}

class GKernelMachineOptimizer : public GTargetFunction
{
protected:
	int m_inputs, m_outputs;
	GKernelMachine* m_pMachine;
	GData* m_pTrainingData;

public:
	GKernelMachineOptimizer(int weights, int inputs, int outputs, GKernelMachine* pMachine, GData* pTrainingData)
	: GTargetFunction(weights), m_inputs(inputs), m_outputs(outputs), m_pMachine(pMachine), m_pTrainingData(pTrainingData)
	{
	}

	virtual ~GKernelMachineOptimizer()
	{
	}

	virtual bool isStable() { return true; }
	virtual bool isConstrained() { return false; }

	virtual void initVector(double* pVector)
	{
		GVec::setAll(pVector, 0.001, relation()->size());
	}

	virtual double computeError(const double* pVector)
	{
		double* pOld = m_pMachine->swapWeights((double*)pVector);
#ifdef WIN32
		GPrediction* out = new GPrediction[m_outputs];
		ArrayHolder<GPrediction> hOut(out);
#else
		GPrediction out[m_outputs];
#endif
		double err = 0;
		for(size_t i = 0; i < m_pTrainingData->rows(); i++)
		{
			double* pPat = m_pTrainingData->row(i);
			m_pMachine->predictDistribution(pPat, out);
			for(int j = 0; j < m_outputs; j++)
			{
				double d = pPat[m_inputs + j] - out[j].asNormal()->mean();
				err += (d * d);
			}
		}
		m_pMachine->swapWeights(pOld);
		return err;
	}
};

// virtual
void GKernelMachine::train(GData* pData, int labelDims)
{
	m_labelDims = labelDims;
	delete[] m_pBuf;
	m_pBuf = new double[labelDims];
	int featureDims = pData->cols() - labelDims;
	m_featureDims = featureDims;
	delete m_pKernel;
	m_pKernel = makeKernel(featureDims);
	if(!pData->relation()->areContinuous(0, pData->relation()->size()))
		ThrowError("Sorry, only continuous attributes are supported");

	// I'm feeling too lazy to compute the support vectors
	// today, so let's just pick 16 patterns
	delete[] m_pWeights;
	m_pWeights = NULL;
	int weightCount = (int)(m_labelDims * pData->rows());
	m_pSupportVectors = new GData(featureDims);
	m_pSupportVectors->newRows(16);
	for(int i = 0; i < 16; i++)
		m_pSupportVectors->copyRow(pData->row((17947 * i) % pData->rows()));

	// I don't feel like coding up a QP solver, so let's just
	// use a hill climber
	GKernelMachineOptimizer critic(weightCount, featureDims, m_labelDims, this, pData);
	GHillClimber search(&critic);
	search.searchUntil(50, 20, 0.005);
	m_pWeights = new double[weightCount];
	GVec::copy(m_pWeights, search.currentVector(), weightCount);
}

// virtual
void GKernelMachine::predictDistribution(const double* pIn, GPrediction* pOut)
{
	if(!m_pWeights)
		ThrowError("Not trained yet");
	GVec::setAll(m_pBuf, 0.0, m_labelDims);
	double* pWeight = m_pWeights;
	for(size_t i = 0; i < m_pSupportVectors->rows(); i++)
	{
		double val = m_pKernel->apply(pIn, m_pSupportVectors->row(i));
		double* pNet = m_pBuf;
		for(int j = 0; j < m_labelDims; j++)
		{
			*pNet += (*pWeight) * val;
			pNet++;
			pWeight++;
		}
	}
	double* pNet = m_pBuf;
	for(int j = 0; j < m_labelDims; j++)
	{
		GNormalDistribution* pNorm = pOut[j].makeNormal();
		pNorm->setMeanAndVariance(GMath::logistic(*pNet), 1.0);
		pNet++;
	}
}

// virtual
void GKernelMachine::predict(const double* pIn, double* pOut)
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


