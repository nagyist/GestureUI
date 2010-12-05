/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#include "GNaiveInstance.h"
#include "GVec.h"
#include "GTwt.h"
#include "GDistribution.h"
#include "GRand.h"
#include "GTransform.h"
#include <map>

using std::multimap;
using std::make_pair;

namespace GClasses {

class GNaiveInstanceAttr
{
protected:
	multimap<double,const double*> m_instances;

public:
	GNaiveInstanceAttr() {}

	GNaiveInstanceAttr(GTwtNode* pAttr, size_t labelDims, GHeap* pHeap)
	{
		size_t count = pAttr->itemCount() / (1 + labelDims);
		if(count * (1 + labelDims) != pAttr->itemCount())
			ThrowError("invalid list size");
		size_t pos = 0;
		for(size_t i = 0; i < count; i++)
		{
			double d = pAttr->item(pos++)->asDouble();
			double* pLabel = (double*)pHeap->allocAligned(sizeof(double) * labelDims);
			m_instances.insert(make_pair(d, pLabel));
			for(size_t j = 0; j < labelDims; j++)
				*(pLabel++) = pAttr->item(pos++)->asDouble();
		}
	}

	virtual ~GNaiveInstanceAttr()
	{
	}

	multimap<double,const double*>& instances() { return m_instances; }

	GTwtNode* toTwt(GTwtDoc* pDoc, size_t labelDims)
	{
		GTwtNode* pList = pDoc->newList((1 + labelDims) * m_instances.size());
		size_t pos = 0;
		for(multimap<double,const double*>::iterator it = m_instances.begin(); it != m_instances.end(); it++)
		{
			pList->setItem(pos++, pDoc->newDouble(it->first));
			for(size_t i = 0; i < labelDims; i++)
				pList->setItem(pos++, pDoc->newDouble(it->second[i]));
		}
		return pList;
	}

	void addInstance(double dInput, const double* pOutputs)
	{
		m_instances.insert(make_pair(dInput, pOutputs));
	}
};

// -----------------------------------------------------------

GNaiveInstance::GNaiveInstance(int nNeighbors)
: GIncrementalLearner(), m_pHeap(NULL)
{
	m_nNeighbors = nNeighbors;
	m_pAttrs = NULL;
	m_labelDims = 0;
	m_featureDims = 0;
	m_pValueSums = NULL;
}

GNaiveInstance::GNaiveInstance(GTwtNode* pNode)
: GIncrementalLearner(pNode), m_pHeap(NULL)
{
	m_pAttrs = NULL;
	m_pValueSums = NULL;
	m_nNeighbors = (int)pNode->field("neighbors")->asInt();
	sp_relation pRel = GRelation::fromTwt(pNode->field("relation"));
	int labelDims = (int)pNode->field("labelDims")->asInt();
	m_featureDims = 0;
	m_labelDims = 0;
	enableIncrementalLearning(pRel, labelDims, NULL, NULL);
	GTwtNode* pAttrs = pNode->field("attrs");
	if((int)pAttrs->itemCount() != m_featureDims)
		ThrowError("Expected ", gformat(m_featureDims), " attrs, got ", gformat(pAttrs->itemCount()), " attrs");
	m_pHeap = new GHeap(1024);
	for(int i = 0; i < m_featureDims; i++)
	{
		delete(m_pAttrs[i]);
		m_pAttrs[i] = new GNaiveInstanceAttr(pAttrs->item(i), m_labelDims, m_pHeap);
	}
}

// virtual
GNaiveInstance::~GNaiveInstance()
{
	clear();
}

void GNaiveInstance::clear()
{
	for(int i = 0; i < m_featureDims; i++)
		delete(m_pAttrs[i]);
	delete[] m_pAttrs;
	m_pAttrs = NULL;
	delete[] m_pValueSums;
	m_pValueSums = NULL;
	m_labelDims = 0;
	m_featureDims = 0;
	m_pRelation.reset();
	delete(m_pHeap);
	m_pHeap = NULL;
}

// virtual
GTwtNode* GNaiveInstance::toTwt(GTwtDoc* pDoc)
{
	GTwtNode* pNode = baseTwtNode(pDoc, "GNaiveInstance");
	pNode->addField(pDoc, "relation", m_pRelation->toTwt(pDoc));
	pNode->addField(pDoc, "labelDims", pDoc->newInt(m_labelDims));
	pNode->addField(pDoc, "neighbors", pDoc->newInt(m_nNeighbors));
	GTwtNode* pAttrs = pNode->addField(pDoc, "attrs", pDoc->newList(m_featureDims));
	for(int i = 0; i < m_featureDims; i++)
		pAttrs->setItem(i, m_pAttrs[i]->toTwt(pDoc, m_labelDims));
	return pNode;
}

// virtual
int GNaiveInstance::featureDims()
{
	if(m_labelDims < 1)
		ThrowError("not yet trained");
	return m_featureDims;
}

// virtual
int GNaiveInstance::labelDims()
{
	if(m_labelDims < 1)
		ThrowError("not yet trained");
	return m_labelDims;
}

// virtual
void GNaiveInstance::enableIncrementalLearning(sp_relation& pRelation, int labelDims, double* pMins, double* pRanges)
{
	if(!pRelation->areContinuous(0, pRelation->size()))
		ThrowError("Only continuous attributes are supported.");
	clear();
	m_pRelation = pRelation;
	m_labelDims = labelDims;
	m_featureDims = pRelation->size() - labelDims;
	m_pAttrs = new GNaiveInstanceAttr*[m_featureDims];
	for(int i = 0; i < m_featureDims; i++)
		m_pAttrs[i] = new GNaiveInstanceAttr();
	m_pValueSums = new double[4 * m_labelDims + m_featureDims];
	m_pWeightSums = &m_pValueSums[m_labelDims];
	m_pSumBuffer = &m_pValueSums[2 * m_labelDims];
	m_pSumOfSquares = &m_pValueSums[3 * m_labelDims];
}

// virtual
void GNaiveInstance::trainIncremental(const double* pIn, const double* pOut)
{
	if(!m_pHeap)
		m_pHeap = new GHeap(1024);
	double* pOutputs = (double*)m_pHeap->allocAligned(sizeof(double) * m_labelDims);
	GVec::copy(pOutputs, pOut, m_labelDims);
	for(int i = 0; i < m_featureDims; i++)
	{
		if(*pIn != UNKNOWN_REAL_VALUE)
			m_pAttrs[i]->addInstance(*(pIn++), pOutputs);
	}
}

// virtual
void GNaiveInstance::train(GData* pData, int labelDims)
{
	enableIncrementalLearning(pData->relation(), labelDims, NULL, NULL);
	double* pPat;
	for(size_t i = 0; i < pData->rows(); i++)
	{
		pPat = pData->row(i);
		trainIncremental(pPat, pPat + m_featureDims);
	}
}

// virtual
void GNaiveInstance::trainSparse(GSparseMatrix* pData, int labelDims)
{
	ThrowError("Sorry, trainSparse is not implemented yet in GNaiveInstance");
}

void GNaiveInstance::evalInput(int nInputDim, double dInput)
{
	// Init the accumulators
	GVec::setAll(m_pSumBuffer, 0.0, m_labelDims);
	GVec::setAll(m_pSumOfSquares, 0.0, m_labelDims);

	// Find the nodes on either side of dInput
	GNaiveInstanceAttr* pAttr = m_pAttrs[nInputDim];
	multimap<double,const double*>& instances = pAttr->instances();
	multimap<double,const double*>::iterator itLeft = instances.lower_bound(dInput);
	multimap<double,const double*>::iterator itRight = itLeft;
	bool leftValid = true;
	if(itLeft == instances.end())
	{
		if(instances.size() > 0)
			itLeft--;
		else
			leftValid = false;
	}
	else
		itRight++;

	// Compute the mean and variance of the values for the k-nearest neighbors
	int nNeighbors = 0;
	bool goRight;
	while(true)
	{
		// Pick the closer of the two nodes
		if(!leftValid)
		{
			if(itRight == instances.end())
				break;
			goRight = true;
		}
		else if(itRight == instances.end())
			goRight = false;
		else if(dInput - itLeft->first < itRight->first - dInput)
			goRight = false;
		else
			goRight = true;

		// Accumulate values
		const double* pOutputVec = goRight ? itRight->second : itLeft->second;
		GVec::add(m_pSumBuffer, pOutputVec, m_labelDims);
		for(int j = 0; j < m_labelDims; j++)
			m_pSumOfSquares[j] += (pOutputVec[j] * pOutputVec[j]);

		// See if we're done
		if(++nNeighbors >= m_nNeighbors)
			break;

		// Advance
		if(goRight)
			itRight++;
		else
		{
			if(itLeft == instances.begin())
				leftValid = false;
			else
				itLeft--;
		}
	}
	GVec::multiply(m_pSumBuffer, 1.0 / nNeighbors, m_labelDims);
	GVec::multiply(m_pSumOfSquares, 1.0 / nNeighbors, m_labelDims);

	// Accumulate the predictions across all dimensions
	int dims = 0;
	double weight;
	for(int i = 0; i < m_labelDims; i++)
	{
		weight = 1.0 / MAX(m_pSumOfSquares[i] - (m_pSumBuffer[i] * m_pSumBuffer[i]), 1e-5);
		m_pWeightSums[dims] += weight;
		m_pValueSums[dims] += weight * m_pSumBuffer[dims];
		dims++;
	}
}

// virtual
void GNaiveInstance::predictDistribution(const double* pIn, GPrediction* pOut)
{
	GVec::setAll(m_pWeightSums, 0.0, m_labelDims);
	GVec::setAll(m_pValueSums, 0.0, m_labelDims);
	int i;
	for(i = 0; i < m_featureDims; i++)
		evalInput(i, pIn[i]);
	for(i = 0; i < m_labelDims; i++)
	{
		GNormalDistribution* pNorm = pOut[i].makeNormal();
		pNorm->setMeanAndVariance(m_pValueSums[i] / m_pWeightSums[i], 1.0 / m_pWeightSums[i]);
	}
}

// virtual
void GNaiveInstance::predict(const double* pIn, double* pOut)
{
	GVec::setAll(m_pWeightSums, 0.0, m_labelDims);
	GVec::setAll(m_pValueSums, 0.0, m_labelDims);
	int i;
	for(i = 0; i < m_featureDims; i++)
		evalInput(i, pIn[i]);
	for(i = 0; i < m_labelDims; i++)
		pOut[i] = m_pValueSums[i] / m_pWeightSums[i];
}

#ifndef NO_TEST_CODE
//static
void GNaiveInstance::test()
{
	GRand prng(0);
	GNaiveInstance ni(8);
	GFilter tl(&ni, false);
	tl.setLabelTransform(new GNominalToCat(), true);
	tl.basicTest(0.71, &prng, 0.02);
}
#endif

} // namespace GClasses
