/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#include "GDecisionTree.h"
#include "GMacros.h"
#include <stdlib.h>
#include "GVec.h"
#include "GPolynomial.h"
#include "GHillClimber.h"
#include "GDistribution.h"
#include "GRand.h"
#include "GTwt.h"
#include "GTransform.h"
#include <string>
#include <iostream>

using namespace GClasses;
using std::cout;
using std::string;
using std::ostream;
using std::vector;

namespace GClasses {

class GDecisionTreeNode
{
public:
	GDecisionTreeNode()
	{
	}

	virtual ~GDecisionTreeNode()
	{
	}

	virtual bool IsLeaf() = 0;
	virtual int GetBranchSize() = 0;
	virtual GDecisionTreeNode* DeepCopy(int nOutputCount, GDecisionTreeNode* pInterestingNode, GDecisionTreeNode** ppOutInterestingCopy) = 0;
	virtual void print(GRelation* pRelation, ostream& stream, int labelDims, int depth, const char* parentValue) = 0;
	virtual void CountValues(int nOutput, int* pnCounts) = 0;
	virtual double FindSumOutputValue(int nOutput) = 0;
	static GDecisionTreeNode* fromTwt(GTwtNode* pNode);
	virtual GTwtNode* toTwt(GTwtDoc* pDoc, int outputCount) = 0;
};

class GDecisionTreeInteriorNode : public GDecisionTreeNode
{
friend class GDecisionTree;
protected:
	int m_nAttribute;
	double m_dPivot;
	int m_nChildren;
	GDecisionTreeNode** m_ppChildren;

public:
	GDecisionTreeInteriorNode(int nAttribute, double dPivot, int children) : GDecisionTreeNode()
	{
		m_nAttribute = nAttribute;
		m_dPivot = dPivot;
		m_nChildren = children;
		m_ppChildren = new GDecisionTreeNode*[children];
		memset(m_ppChildren, '\0', sizeof(GDecisionTreeNode*) * children);
	}

	GDecisionTreeInteriorNode(GTwtNode* pNode) : GDecisionTreeNode()
	{
		m_nAttribute = (int)pNode->field("attr")->asInt();
		m_dPivot = pNode->field("pivot")->asDouble();
		GTwtNode* pChildren = pNode->field("children");
		m_nChildren = (int)pChildren->itemCount();
		m_ppChildren = new GDecisionTreeNode*[m_nChildren];
		int i;
		for(i = 0; i < m_nChildren; i++)
			m_ppChildren[i] = GDecisionTreeNode::fromTwt(pChildren->item(i));
	}

	virtual ~GDecisionTreeInteriorNode()
	{
		if(m_ppChildren)
		{
			int n;
			for(n = 0; n < m_nChildren; n++)
				delete(m_ppChildren[n]);
			delete[] m_ppChildren;
		}
	}

	virtual GTwtNode* toTwt(GTwtDoc* pDoc, int outputCount)
	{
		GTwtNode* pNode = pDoc->newObj();
		pNode->addField(pDoc, "attr", pDoc->newInt(m_nAttribute));
		pNode->addField(pDoc, "pivot", pDoc->newDouble(m_dPivot));
		GTwtNode* pChildren = pDoc->newList(m_nChildren);
		pNode->addField(pDoc, "children", pChildren);
		int i;
		for(i = 0; i < m_nChildren; i++)
			pChildren->setItem(i, m_ppChildren[i]->toTwt(pDoc, outputCount));
		return pNode;
	}

	virtual bool IsLeaf() { return false; }

	virtual int GetBranchSize()
	{
		int size = 1;
		int i;
		for(i = 0; i < m_nChildren; i++)
			size += m_ppChildren[i]->GetBranchSize();
		return size;
	}

	virtual GDecisionTreeNode* DeepCopy(int nOutputCount, GDecisionTreeNode* pInterestingNode, GDecisionTreeNode** ppOutInterestingCopy)
	{
		GDecisionTreeInteriorNode* pNewNode = new GDecisionTreeInteriorNode(m_nAttribute, m_dPivot, m_nChildren);
		for(int n = 0; n < m_nChildren; n++)
			pNewNode->m_ppChildren[n] = m_ppChildren[n]->DeepCopy(nOutputCount, pInterestingNode, ppOutInterestingCopy);
		if(this == pInterestingNode)
			*ppOutInterestingCopy = pNewNode;
		return pNewNode;
	}

	virtual void print(GRelation* pRelation, ostream& stream, int labelDims, int depth, const char* parentValue)
	{
		for(int n = 0; n < depth; n++)
			stream << "  ";
		if(parentValue)
			stream << parentValue << " -> ";
		if(pRelation->valueCount(m_nAttribute) == 0)
		{
			string s;
			pRelation->attrValue(&s, m_nAttribute, m_dPivot);
			if(pRelation->type() == GRelation::ARFF)
				stream << "Is " << ((GArffRelation*)pRelation)->attrName(m_nAttribute) << " < " << s.c_str() << "?\n";
			else
				stream << "Is attr " << m_nAttribute << " < " << s.c_str() << "?\n";
			if(m_nChildren != 2)
				ThrowError("expected this node to have two child nodes");
			m_ppChildren[0]->print(pRelation, stream, labelDims, depth + 1, "Yes");
			m_ppChildren[1]->print(pRelation, stream, labelDims, depth + 1, "No");
		}
		else
		{
			if(pRelation->type() == GRelation::ARFF)
				stream << "What is the value of " << ((GArffRelation*)pRelation)->attrName(m_nAttribute) << "?\n";
			else
				stream << "What is the value of attr " << m_nAttribute << "?\n";
			for(int n = 0; n < m_nChildren; n++)
			{
				string s;
				pRelation->attrValue(&s, m_nAttribute, (double)n);
				m_ppChildren[n]->print(pRelation, stream, labelDims, depth + 1, s.c_str());
			}
		}
	}

	// Recursive function that counts the number of times a particular
	// value is found in a particular output in this branch of the tree
	virtual void CountValues(int nOutput, int* pnCounts)
	{
		int n;
		for(n = 0; n < m_nChildren; n++)
			m_ppChildren[n]->CountValues(nOutput, pnCounts);
	}

	virtual double FindSumOutputValue(int nOutput)
	{
		double dSum = 0;
		int n;
		for(n = 0; n < m_nChildren; n++)
			dSum += m_ppChildren[n]->FindSumOutputValue(nOutput);
		return dSum;
	}
};

class GDecisionTreeLeafNode : public GDecisionTreeNode
{
public:
	double* m_pOutputValues;
	size_t m_nSampleSize;

public:
	GDecisionTreeLeafNode(double* pOutputValues, size_t nSampleSize) : GDecisionTreeNode()
	{
		m_pOutputValues = pOutputValues;
		m_nSampleSize = nSampleSize;
	}

	GDecisionTreeLeafNode(GTwtNode* pNode) : GDecisionTreeNode()
	{
		m_nSampleSize = (size_t)pNode->field("size")->asInt();
		GTwtNode* pOut = pNode->field("out");
		size_t count = pOut->itemCount();
		m_pOutputValues = new double[count];
		for(size_t i = 0; i < count; i++)
			m_pOutputValues[i] = pOut->item(i)->asDouble();
	}

	virtual ~GDecisionTreeLeafNode()
	{
		delete[] m_pOutputValues;
	}

	virtual GTwtNode* toTwt(GTwtDoc* pDoc, int outputCount)
	{
		GTwtNode* pNode = pDoc->newObj();
		pNode->addField(pDoc, "size", pDoc->newInt(m_nSampleSize));
		GTwtNode* pOut = pDoc->newList(outputCount);
		pNode->addField(pDoc, "out", pOut);
		for(int i = 0; i < outputCount; i++)
			pOut->setItem(i, pDoc->newDouble(m_pOutputValues[i]));
		return pNode;
	}

	virtual bool IsLeaf() { return true; }

	virtual int GetBranchSize()
	{
		return 1;
	}

	virtual GDecisionTreeNode* DeepCopy(int nOutputCount, GDecisionTreeNode* pInterestingNode, GDecisionTreeNode** ppOutInterestingCopy)
	{
		double* pOutputValues = new double[nOutputCount];
		GVec::copy(pOutputValues, m_pOutputValues, nOutputCount);
		GDecisionTreeLeafNode* pNewNode = new GDecisionTreeLeafNode(pOutputValues, m_nSampleSize);
		if(this == pInterestingNode)
			*ppOutInterestingCopy = pNewNode;
		return pNewNode;
	}

	virtual void print(GRelation* pRelation, ostream& stream, int labelDims, int depth, const char* parentValue)
	{
		for(int n = 0; n < depth; n++)
			stream << "  ";
		if(parentValue)
			stream << parentValue << " -> ";
		int featureDims = pRelation->size() - labelDims;
		for(int n = 0; n < labelDims; n++)
		{
			int labelAttr = featureDims + n;
			if(n > 0)
				stream << ", ";
			string s;
			pRelation->attrValue(&s, labelAttr, m_pOutputValues[n]);
			if(pRelation->type() == GRelation::ARFF)
				stream << ((GArffRelation*)pRelation)->attrName(labelAttr) << "=" << s.c_str();
			else
				stream << s.c_str();
		}
		stream << "\n";
	}

	virtual void CountValues(int nOutput, int* pnCounts)
	{
		int nVal = (int)m_pOutputValues[nOutput];
		pnCounts[nVal] += m_nSampleSize;
	}

	virtual double FindSumOutputValue(int nOutput)
	{
		return m_pOutputValues[nOutput] * m_nSampleSize;
	}
};

}

// static
GDecisionTreeNode* GDecisionTreeNode::fromTwt(GTwtNode* pNode)
{
	if(pNode->fieldIfExists("children"))
		return new GDecisionTreeInteriorNode(pNode);
	else
		return new GDecisionTreeLeafNode(pNode);
}

// -----------------------------------------------------------------

GDecisionTree::GDecisionTree(GRand* pRand)
: GSupervisedLearner(), m_leafThresh(1), m_featureDims(0), m_labelDims(0)
{
	m_pRoot = NULL;
	m_eAlg = GDecisionTree::MINIMIZE_ENTROPY;
	m_pRand = pRand;
}

GDecisionTree::GDecisionTree(GDecisionTree* pThat, GDecisionTreeNode* pInterestingNode, GDecisionTreeNode** ppOutInterestingCopy)
: GSupervisedLearner(), m_leafThresh(1)
{
	if(pThat->m_labelDims == 0)
		ThrowError("not trained");
	m_pRoot = pThat->m_pRoot->DeepCopy(pThat->m_labelDims, pInterestingNode, ppOutInterestingCopy);
}

GDecisionTree::GDecisionTree(GTwtNode* pNode, GRand* pRand)
: GSupervisedLearner(pNode), m_pRand(pRand), m_leafThresh(1)
{
	m_pRelation = GRelation::fromTwt(pNode->field("relation"));
	m_labelDims = (int)pNode->field("labelDims")->asInt();
	m_featureDims = m_pRelation->size() - m_labelDims;
	m_eAlg = (DivisionAlgorithm)pNode->field("alg")->asInt();
	m_pRoot = GDecisionTreeNode::fromTwt(pNode->field("root"));
}

// virtual
GDecisionTree::~GDecisionTree()
{
	clear();
}

// virtual
GTwtNode* GDecisionTree::toTwt(GTwtDoc* pDoc)
{
	if(!m_pRoot)
		ThrowError("not trained yet");
	GTwtNode* pNode = baseTwtNode(pDoc, "GDecisionTree");
	pNode->addField(pDoc, "relation", m_pRelation->toTwt(pDoc));
	pNode->addField(pDoc, "labelDims", pDoc->newInt(m_labelDims));
	pNode->addField(pDoc, "alg", pDoc->newInt(m_eAlg));
	pNode->addField(pDoc, "root", m_pRoot->toTwt(pDoc, m_labelDims));
	return pNode;
}

// virtual
int GDecisionTree::featureDims()
{
	if(m_labelDims < 1)
		ThrowError("not yet trained");
	return m_featureDims;
}

// virtual
int GDecisionTree::labelDims()
{
	if(m_labelDims < 1)
		ThrowError("not yet trained");
	return m_labelDims;
}

int GDecisionTree::treeSize()
{
	return m_pRoot->GetBranchSize();
}

void GDecisionTree::print(ostream& stream, GArffRelation* pRelation)
{
	if(!m_pRoot)
		ThrowError("not trained yet");
	GRelation* pRel = pRelation;
	if(!pRel)
		pRel = m_pRelation.get();
	m_pRoot->print(pRel, stream, m_labelDims, 0, NULL);
}

void GDecisionTree::train(GData* pData, int labelDims)
{
	m_pRelation = pData->relation();
	clear();
	m_labelDims = labelDims;
	m_featureDims = pData->cols() - m_labelDims;

	// Make a list of available features
	vector<size_t> attrPool;
	attrPool.reserve(m_featureDims);
	for(size_t i = 0; i < (size_t)m_featureDims; i++)
		attrPool.push_back(i);

	// We must make a copy of the data because buildNode will mess with it
	// by calling RandomlyReplaceMissingData
	GData tmp(pData->relation(), pData->heap());
	tmp.copy(pData);

	m_pRoot = buildBranch(&tmp, attrPool, 0/*depth*/, 4/*tolerance*/);
}

double GDecisionTree_measureRealSplitInfo(GData* pData, GData* pOther, int labelDims, int attr, double pivot)
{
	GAssert(pOther->rows() == 0);
	size_t rowCount = pData->rows();
	pData->splitByPivot(pOther, attr, pivot);
	double d;
	if(pData->rows() >= 1 && pOther->rows() >= 1)
		d = (pData->measureLabelInfo(labelDims) * pData->rows() + pOther->measureLabelInfo(labelDims) * pOther->rows()) / rowCount;
	else
		d = 1e308;
	pData->mergeVert(pOther);
	return d;
}

double GDecisionTree_pickPivotToReduceInfo(GData* pData, double* pPivot, int labelDims, int attr, GRand* pRand)
{
	size_t nRows = pData->rows();
	double bestPivot = UNKNOWN_REAL_VALUE;
	double bestInfo = 1e100;
	double* pRow1;
	double* pRow2;
	GData dataTmp(pData->relation(), pData->heap());
	dataTmp.reserve(nRows);
	size_t attempts = MIN(pData->rows() - 1, (pData->rows() * pData->cols() > 100000 ? (size_t)1 : (size_t)8));
	for(size_t n = 0; n < attempts; n++)
	{
		pRow1 = pData->row((size_t)pRand->next(nRows));
		pRow2 = pData->row((size_t)pRand->next(nRows));
		double pivot = 0.5 * (pRow1[attr] + pRow2[attr]);
		double info = GDecisionTree_measureRealSplitInfo(pData, &dataTmp, labelDims, attr, pivot);
		if(info < bestInfo)
		{
			bestInfo = info;
			bestPivot = pivot;
		}
	}
	*pPivot = bestPivot;
	return bestInfo;
}

double GDecisionTree_measureNominalSplitInfo(GData* pData, int nAttribute, int labelDims)
{
	size_t nRowCount = pData->rows() - pData->countValue(nAttribute, UNKNOWN_DISCRETE_VALUE);
	GData dataClass(pData->relation(), pData->heap());
	int values = pData->relation()->valueCount(nAttribute);
	double dInfo = 0;
	for(int n = 0; n < values; n++)
	{
		pData->splitByDiscreteValue(&dataClass, nAttribute, n);
		dInfo += ((double)dataClass.rows() / nRowCount) * dataClass.measureLabelInfo(labelDims);
		pData->mergeVert(&dataClass);
	}
	return dInfo;
}

size_t GDecisionTree::pickDivision(GData* pData, double* pPivot, vector<size_t>& attrPool, int nDepth)
{
	if(m_eAlg == MINIMIZE_ENTROPY)
	{
		// Pick the best attribute to divide on
		GAssert(pData->rows() > 0); // Can't work without data
		double bestInfo = 1e100;
		double pivot = 0.0;
		double bestPivot = 0;
		size_t index = 0;
		size_t bestIndex = attrPool.size();
		for(vector<size_t>::iterator it = attrPool.begin(); it != attrPool.end(); it++)
		{
			double info;
			if(pData->relation()->valueCount(*it) == 0)
				info = GDecisionTree_pickPivotToReduceInfo(pData, &pivot, m_labelDims, *it, m_pRand);
			else
				info = GDecisionTree_measureNominalSplitInfo(pData, *it, m_labelDims);
			if(info < bestInfo)
			{
				bestInfo = info;
				bestIndex = index;
				bestPivot = pivot;
			}
			index++;
		}
		*pPivot = bestPivot;
		return bestIndex;
	}
	else if(m_eAlg == RANDOM)
	{
		// Pick the best of m_randomDraws random attributes from the attribute pool
		GAssert(pData->rows() > 0); // Can't work without data
		double bestInfo = 1e100;
		double bestPivot = 0;
		size_t bestIndex = attrPool.size();
		size_t patience = MAX((size_t)6, m_randomDraws * 2);
		GData dataTmp(pData->relation(), pData->heap());
		dataTmp.reserve(pData->rows());
		for(size_t i = 0; i < m_randomDraws && patience > 0; i++)
		{
			size_t index = (size_t)m_pRand->next(attrPool.size());
			size_t attr = attrPool[index];
			double pivot = 0.0;
			double info;
			if(pData->relation()->valueCount(attr) == 0)
			{
				// Pick a random pivot. (Note that this is not a uniform distribution. This
				// distribution is biased in favor of pivots that will divide the data well.)
				double a = pData->row((size_t)m_pRand->next(pData->rows()))[attr];
				double b = pData->row((size_t)m_pRand->next(pData->rows()))[attr];
				pivot = 0.5 * (a + b);
				if(m_randomDraws > 1)
					info = GDecisionTree_measureRealSplitInfo(pData, &dataTmp, m_labelDims, attr, pivot);
				else
					info = 0.0;
			}
			else
			{
				if(m_randomDraws > 1)
					info = GDecisionTree_measureNominalSplitInfo(pData, attr, m_labelDims);
				else
					info = 0.0;
			}
			if(info < bestInfo)
			{
				bestInfo = info;
				bestIndex = index;
				bestPivot = pivot;
			}
		}
		if(bestIndex < attrPool.size() && !pData->isAttrHomogenous(attrPool[bestIndex]))
		{
			*pPivot = bestPivot;
			return bestIndex;
		}

		// We failed to find a useful attribute with random draws. (This may happen if there is a large
		// ratio of homogenous attributes.) Now, we need to be a little more systematic about finding a good
		// attribute. (This is not specified in the random forest algorithm, but it makes a big difference
		// with some problems.)
		size_t k = (size_t)m_pRand->next(attrPool.size());
		for(size_t i = 0; i < attrPool.size(); i++)
		{
			size_t index = (i + k) % attrPool.size();
			size_t attr = attrPool[index];
			if(pData->relation()->valueCount(attr) == 0)
			{
				// Find the min
				double m = pData->row(0)[attr];
				for(size_t j = 1; j < pData->rows(); j++)
				{
					double d = pData->row(j)[attr];
					if(d != UNKNOWN_REAL_VALUE)
						m = MIN(m, d);
				}

				// Randomly pick one of the non-min values
				size_t candidates = 0;
				for(size_t j = 0; j < pData->rows(); j++)
				{
					double d = pData->row(j)[attr];
					if(d != UNKNOWN_REAL_VALUE && d > m)
					{
						if(m_pRand->next(++candidates) == 0)
							*pPivot = d;
					}
				}
				if(candidates == 0)
					continue; // This attribute is worthless
			}
			else
			{
				if(pData->isAttrHomogenous(attr))
					continue; // This attribute is worthless
			}
			return index;
		}
	}
	else
		GAssert(false); // unknown division algorithm
	return attrPool.size();
}

double* GDecisionTreeNode_labelVec(GData* pData, int labelDims)
{
	double* pVec = new double[labelDims];
	pData->baselineVector(labelDims, pVec);
	return pVec;
}

double* GDecisionTreeNode_labelVec(GData* pData, GDataArray* pParts, int labelDims)
{
	double* pVec = new double[labelDims];
	GData* pB = pParts->sets()[pParts->largestSet()];
	if(pData->rows() > pB->rows())
		pData->baselineVector(labelDims, pVec);
	else
		pB->baselineVector(labelDims, pVec);
	return pVec;
}

// This constructs the decision tree in a recursive depth-first manner
GDecisionTreeNode* GDecisionTree::buildBranch(GData* pData, vector<size_t>& attrPool, int nDepth, int tolerance)
{
	// Make a leaf if we're out of tolerance or the output is homogenous or there are no attributes left
	if(tolerance <= 0 || pData->rows() <= m_leafThresh || attrPool.size() == 0 || pData->areLabelsHomogenous(m_labelDims))
		return new GDecisionTreeLeafNode(GDecisionTreeNode_labelVec(pData, m_labelDims), pData->rows());

	// Pick the division
	double pivot = 0.0;
	size_t bestIndex = pickDivision(pData, &pivot, attrPool, nDepth);

	// Make a leaf if there are no good divisions
	if(bestIndex >= attrPool.size())
		return new GDecisionTreeLeafNode(GDecisionTreeNode_labelVec(pData, m_labelDims), pData->rows());
	size_t attr = attrPool[bestIndex];

	// Make sure there aren't any missing values in the decision attribute
	pData->randomlyReplaceMissingValues(attr, m_pRand);

	// Split the data
	GDataArray parts(pData->relation());
	int nonEmptyBranchCount = 0;
	if(pData->relation()->valueCount(attr) == 0)
	{
		// Split on a continuous attribute
		GData* pOther = parts.newSet(0);
		pData->splitByPivot(pOther, attr, pivot);
		nonEmptyBranchCount += (pData->rows() > 0 ? 1 : 0) + (pOther->rows() > 0 ? 1 : 0);
	}
	else
	{
		// Split on a nominal attribute
		int valueCount = pData->relation()->valueCount(attr);
		for(int i = 1; i < valueCount; i++)
		{
			GData* pOther = parts.newSet(0);
			pData->splitByDiscreteValue(pOther, attr, i);
			if(pOther->rows() > 0)
				nonEmptyBranchCount++;
		}
		if(pData->rows() > 0)
			nonEmptyBranchCount++;

		// Remove this attribute from the pool of available attributes
		std::swap(attrPool[bestIndex], attrPool[attrPool.size() - 1]);
		attrPool.erase(attrPool.end() - 1);
	}

	// If we didn't actually separate anything
	if(nonEmptyBranchCount < 2)
	{
		size_t setCount = parts.sets().size();
		for(size_t i = 0; i < setCount; i++)
			pData->mergeVert(parts.sets()[i]);
		if(m_eAlg == MINIMIZE_ENTROPY)
			return new GDecisionTreeLeafNode(GDecisionTreeNode_labelVec(pData, m_labelDims), pData->rows());
		else
		{
			// Try another division
			GDecisionTreeNode* pNode = buildBranch(pData, attrPool, nDepth, tolerance - 1);
			attrPool.push_back(attr);
			return pNode;
		}
	}

	// Make an interior node
	GDecisionTreeInteriorNode* pNode = new GDecisionTreeInteriorNode(attr, pivot, parts.sets().size() + 1);
	Holder<GDecisionTreeInteriorNode> hNode(pNode);
	if(pData->rows() > 0)
		pNode->m_ppChildren[0] = buildBranch(pData, attrPool, nDepth + 1, tolerance);
	else
		pNode->m_ppChildren[0] = new GDecisionTreeLeafNode(GDecisionTreeNode_labelVec(pData, &parts, m_labelDims), 0);
	for(size_t i = 0; i < parts.sets().size(); i++)
	{
		if(parts.sets()[i]->rows() > 0)
			pNode->m_ppChildren[i + 1] = buildBranch(parts.sets()[i], attrPool, nDepth + 1, tolerance);
		else
			pNode->m_ppChildren[i + 1] = new GDecisionTreeLeafNode(GDecisionTreeNode_labelVec(pData, &parts, m_labelDims), 0);
	}
	attrPool.push_back(attr);
	return hNode.release();
}

GDecisionTreeLeafNode* GDecisionTree::findLeaf(const double* pIn, int* pDepth)
{
	if(!m_pRoot)
		ThrowError("Not trained yet");
	GDecisionTreeNode* pNode = m_pRoot;
	int nVal;
	int nDepth = 1;
	while(!pNode->IsLeaf())
	{
		GDecisionTreeInteriorNode* pInterior = (GDecisionTreeInteriorNode*)pNode;
		if(m_pRelation->valueCount(pInterior->m_nAttribute) == 0)
		{
			if(pIn[pInterior->m_nAttribute] == UNKNOWN_REAL_VALUE)
				pNode = pInterior->m_ppChildren[m_pRand->next(2)]; // todo: we could pick in proportion to the number of training patterns that took each side
			else if(pIn[pInterior->m_nAttribute] < pInterior->m_dPivot)
				pNode = pInterior->m_ppChildren[0];
			else
				pNode = pInterior->m_ppChildren[1];
		}
		else
		{
			nVal = (int)pIn[pInterior->m_nAttribute];
			if(nVal < 0)
			{
				GAssert(nVal == UNKNOWN_DISCRETE_VALUE); // out of range
				nVal = (int)m_pRand->next(m_pRelation->valueCount(pInterior->m_nAttribute));
			}
			GAssert(nVal < m_pRelation->valueCount(pInterior->m_nAttribute)); // value out of range
			pNode = pInterior->m_ppChildren[nVal];
		}
		nDepth++;
	}
	*pDepth = nDepth;
	return (GDecisionTreeLeafNode*)pNode;
}


void GDecisionTree::predict(const double* pIn, double* pOut)
{
	int depth;
	GDecisionTreeLeafNode* pLeaf = findLeaf(pIn, &depth);
	GVec::copy(pOut, pLeaf->m_pOutputValues, m_labelDims);
}

void GDecisionTree::predictDistribution(const double* pIn, GPrediction* pOut)
{
	// Copy the output values into the row
	int depth;
	GDecisionTreeLeafNode* pLeaf = findLeaf(pIn, &depth);
	int n, nValues;
	for(n = 0; n < m_labelDims; n++)
	{
		nValues = m_pRelation->valueCount(m_featureDims + n);
		if(nValues == 0)
			pOut[n].makeNormal()->setMeanAndVariance(pLeaf->m_pOutputValues[n], depth);
		else
			pOut[n].makeCategorical()->setSpike(nValues, (int)pLeaf->m_pOutputValues[n], depth);
	}
}

// virtual
void GDecisionTree::clear()
{
	delete(m_pRoot);
	m_pRoot = NULL;
}

#ifndef NO_TEST_CODE
// static
void GDecisionTree::test()
{
	GRand prng(0);
	GDecisionTree tree(&prng);
	tree.basicTest(0.64, &prng);
}
#endif

// ----------------------------------------------------------------------

namespace GClasses {
class GMeanMarginsTreeNode
{
public:
	GMeanMarginsTreeNode()
	{
	}

	virtual ~GMeanMarginsTreeNode()
	{
	}

	virtual bool IsLeaf() = 0;
	virtual GTwtNode* toTwt(GTwtDoc* pDoc, int nInputs, int nOutputs) = 0;

	static GMeanMarginsTreeNode* fromTwt(GTwtNode* pNode);
};


class GMeanMarginsTreeInteriorNode : public GMeanMarginsTreeNode
{
protected:
	double* m_pCenter;
	double* m_pNormal;
	GMeanMarginsTreeNode* m_pLeft;
	GMeanMarginsTreeNode* m_pRight;

public:
	GMeanMarginsTreeInteriorNode(double* pCenter, double* pNormal)
	: GMeanMarginsTreeNode()
	{
		m_pCenter = pCenter;
		m_pNormal = pNormal;
		m_pLeft = NULL;
		m_pRight = NULL;
	}

	GMeanMarginsTreeInteriorNode(GTwtNode* pNode)
	: GMeanMarginsTreeNode()
	{
		GTwtNode* pCenter = pNode->field("center");
		int dims = (int)pCenter->itemCount();
		m_pCenter = new double[dims];
		GVec::fromTwt(m_pCenter, dims, pCenter);
		m_pNormal = new double[dims];
		GVec::fromTwt(m_pNormal, dims, pNode->field("normal"));
		m_pLeft = GMeanMarginsTreeNode::fromTwt(pNode->field("left"));
		m_pRight = GMeanMarginsTreeNode::fromTwt(pNode->field("right"));
	}

	virtual ~GMeanMarginsTreeInteriorNode()
	{
		delete[] m_pCenter;
		delete[] m_pNormal;
		delete(m_pLeft);
		delete(m_pRight);
	}

	virtual GTwtNode* toTwt(GTwtDoc* pDoc, int nInputs, int nOutputs)
	{
		GTwtNode* pNode = pDoc->newObj();
		pNode->addField(pDoc, "center", GVec::toTwt(pDoc, m_pCenter, nInputs));
		pNode->addField(pDoc, "normal", GVec::toTwt(pDoc, m_pNormal, nInputs));
		pNode->addField(pDoc, "left", m_pLeft->toTwt(pDoc, nInputs, nOutputs));
		pNode->addField(pDoc, "right", m_pRight->toTwt(pDoc, nInputs, nOutputs));
		return pNode;
	}

	virtual bool IsLeaf()
	{
		return false;
	}

	bool Test(double* pInputVector, int nInputs)
	{
		return GVec::dotProductIgnoringUnknowns(m_pCenter, pInputVector, m_pNormal, nInputs) >= 0;
	}

	void SetLeft(GMeanMarginsTreeNode* pNode)
	{
		m_pLeft = pNode;
	}

	void SetRight(GMeanMarginsTreeNode* pNode)
	{
		m_pRight = pNode;
	}

	GMeanMarginsTreeNode* GetRight()
	{
		return m_pRight;
	}

	GMeanMarginsTreeNode* GetLeft()
	{
		return m_pLeft;
	}
};

class GMeanMarginsTreeLeafNode : public GMeanMarginsTreeNode
{
protected:
	double* m_pOutputs;

public:
	GMeanMarginsTreeLeafNode(int nOutputCount, double* pOutputs)
	: GMeanMarginsTreeNode()
	{
		m_pOutputs = new double[nOutputCount];
		GVec::copy(m_pOutputs, pOutputs, nOutputCount);
	}

	GMeanMarginsTreeLeafNode(GTwtNode* pNode)
	: GMeanMarginsTreeNode()
	{
		int dims = (int)pNode->itemCount();
		m_pOutputs = new double[dims];
		GVec::fromTwt(m_pOutputs, dims, pNode);
	}

	virtual ~GMeanMarginsTreeLeafNode()
	{
		delete[] m_pOutputs;
	}

	virtual GTwtNode* toTwt(GTwtDoc* pDoc, int nInputs, int nOutputs)
	{
		return GVec::toTwt(pDoc, m_pOutputs, nOutputs);
	}

	virtual bool IsLeaf()
	{
		return true;
	}

	double* GetOutputs()
	{
		return m_pOutputs;
	}
};
}

// static
GMeanMarginsTreeNode* GMeanMarginsTreeNode::fromTwt(GTwtNode* pNode)
{
	if(pNode->type() == GTwtNode::type_list)
		return new GMeanMarginsTreeLeafNode(pNode);
	else
	{
		// todo: if it's a GMVNTreeInteriorNode { ... } else
		return new GMeanMarginsTreeInteriorNode(pNode);
	}
}

// ---------------------------------------------------------------

GMeanMarginsTree::GMeanMarginsTree(GRand* pRand)
: GSupervisedLearner(), m_labelDims(0), m_pRand(pRand)
{
	m_pRoot = NULL;
	m_pEvalVector = NULL;
}

GMeanMarginsTree::GMeanMarginsTree(GTwtNode* pNode, GRand* pRand)
: GSupervisedLearner(pNode), m_pRand(pRand)
{
	m_pRoot = GMeanMarginsTreeNode::fromTwt(pNode->field("root"));
	m_pRelation = GRelation::fromTwt(pNode->field("relation"));
	m_labelDims = (int)pNode->field("labelDims")->asInt();
	m_nInputVectorSize = (int)pNode->field("in")->asInt();
	m_nOutputVectorSize = (int)pNode->field("out")->asInt();
	m_pEvalVector = new double[m_nInputVectorSize];
}

GMeanMarginsTree::~GMeanMarginsTree()
{
	delete(m_pRoot);
	delete[] m_pEvalVector;
}

// virtual
GTwtNode* GMeanMarginsTree::toTwt(GTwtDoc* pDoc)
{
	GTwtNode* pNode = baseTwtNode(pDoc, "GMeanMarginsTree");
	pNode->addField(pDoc, "relation", m_pRelation->toTwt(pDoc));
	pNode->addField(pDoc, "labelDims", pDoc->newInt(m_labelDims));
	pNode->addField(pDoc, "in", pDoc->newInt(m_nInputVectorSize));
	pNode->addField(pDoc, "out", pDoc->newInt(m_nOutputVectorSize));
	pNode->addField(pDoc, "root", m_pRoot->toTwt(pDoc, m_nInputVectorSize, m_nOutputVectorSize));
	return pNode;
}

// virtual
int GMeanMarginsTree::featureDims()
{
	if(m_labelDims < 1)
		ThrowError("not yet trained");
	return m_pRelation->size() - m_labelDims;
}

// virtual
int GMeanMarginsTree::labelDims()
{
	if(m_labelDims < 1)
		ThrowError("not yet trained");
	return m_labelDims;
}

// virtual
void GMeanMarginsTree::train(GData* pData, int labelDims)
{
	clear();
	m_pRelation = pData->relation();
	m_labelDims = labelDims;
	int featureDims = pData->cols() - m_labelDims;
	m_nInputVectorSize = pData->relation()->countRealSpaceDims(0, featureDims);
	m_nOutputVectorSize = pData->relation()->countRealSpaceDims(featureDims, m_labelDims);
	m_pEvalVector = new double[m_nInputVectorSize];

	int nCount = (int)pData->rows();
	GHeap heap(1000);
	GData internalData(m_nInputVectorSize + m_nOutputVectorSize, &heap);
	internalData.reserve(nCount);
	int i;
	double* pVectorIn;
	double* pVectorOut;
	for(i = 0; i < nCount; i++)
	{
		pVectorIn = pData->row(i);
		pVectorOut = internalData.newRow();
		m_pRelation->toRealSpace(pVectorIn, pVectorOut, 0, featureDims + m_labelDims);
		if(GVec::doesContainUnknowns(pVectorOut + m_nInputVectorSize, m_nOutputVectorSize))
			ThrowError("GMeanMarginsTree doesn't support unknown real output values");
	}
	double* pBuf = new double[2 * m_nOutputVectorSize + 3 * m_nInputVectorSize];
	ArrayHolder<double> hBuf(pBuf);
	m_pRoot = buildNode(&internalData, pBuf);
}

// virtual
void GMeanMarginsTree::predictDistribution(const double* pIn, GPrediction* pOut)
{
	int featureDims = m_pRelation->size() - m_labelDims;
	m_pRelation->toRealSpace(pIn, m_pEvalVector, 0, featureDims);
	GMeanMarginsTreeNode* pNode = m_pRoot;
	int nDepth = 1;
	while(!pNode->IsLeaf())
	{
		if(((GMeanMarginsTreeInteriorNode*)pNode)->Test(m_pEvalVector, m_nInputVectorSize))
			pNode = ((GMeanMarginsTreeInteriorNode*)pNode)->GetRight();
		else
			pNode = ((GMeanMarginsTreeInteriorNode*)pNode)->GetLeft();
		nDepth++;
	}
	m_pRelation->fromRealSpace(((GMeanMarginsTreeLeafNode*)pNode)->GetOutputs(), pOut, featureDims, m_labelDims);
}

// virtual
void GMeanMarginsTree::predict(const double* pIn, double* pOut)
{
	int featureDims = m_pRelation->size() - m_labelDims;
	m_pRelation->toRealSpace(pIn, m_pEvalVector, 0, featureDims);
	GMeanMarginsTreeNode* pNode = m_pRoot;
	int nDepth = 1;
	while(!pNode->IsLeaf())
	{
		if(((GMeanMarginsTreeInteriorNode*)pNode)->Test(m_pEvalVector, m_nInputVectorSize))
			pNode = ((GMeanMarginsTreeInteriorNode*)pNode)->GetRight();
		else
			pNode = ((GMeanMarginsTreeInteriorNode*)pNode)->GetLeft();
		nDepth++;
	}
	m_pRelation->fromRealSpace(((GMeanMarginsTreeLeafNode*)pNode)->GetOutputs(), pOut, featureDims, m_labelDims, m_pRand);
}

class GDataShifter
{
protected:
	GData* m_pData;
	int m_featureDims;

public:
	GDataShifter(GData* pData, int nInputCount)
	{
		m_pData = pData;
		m_featureDims = nInputCount;

		// Shift out the input data
		size_t nCount = m_pData->rows();
		for(size_t i = 0; i < nCount; i++)
			m_pData->replaceRow(i, m_pData->row(i) + m_featureDims);
	}

	~GDataShifter()
	{
		// Shift the input data back in
		size_t nCount = m_pData->rows();
		for(size_t i = 0; i < nCount; i++)
			m_pData->replaceRow(i, m_pData->row(i) - m_featureDims);
	}
};

GMeanMarginsTreeNode* GMeanMarginsTree::buildNode(GData* pData, double* pBuf)
{
	// Check for a leaf node
	size_t nCount = pData->rows();
	if(nCount < 2)
	{
		GAssert(nCount > 0); // no data left
		return new GMeanMarginsTreeLeafNode(m_nOutputVectorSize, pData->row(0) + m_nInputVectorSize);
	}

	double* pOutputMean = pBuf;
	double* pPrincipleComponent = pBuf + m_nOutputVectorSize;
	{
		GDataShifter shifter(pData, m_nInputVectorSize);

		// Compute the output mean and principle component
		for(int i = 0; i < m_nOutputVectorSize; i++)
			pOutputMean[i] = pData->mean(i);
		pData->principalComponentIgnoreUnknowns(pPrincipleComponent, m_nOutputVectorSize, pOutputMean, m_pRand);
	}

	// Find the input mean of each cluster
	double* pInputMean1 = pPrincipleComponent + m_nOutputVectorSize;
	double* pInputMean2 = pInputMean1 + m_nInputVectorSize;
	GAssert(sizeof(double) == 2 * sizeof(int)); // doubles have unexpected size
	int* pCounts1 = (int*)(pInputMean2 + m_nInputVectorSize);
	int* pCounts2 = pCounts1 + m_nInputVectorSize;
	GVec::setAll(pInputMean1, 0.0, m_nInputVectorSize);
	GVec::setAll(pInputMean2, 0.0, m_nInputVectorSize);
	memset(pCounts1, '\0', sizeof(int) * m_nInputVectorSize);
	memset(pCounts2, '\0', sizeof(int) * m_nInputVectorSize);
	double* pVector;
	for(size_t i = 0; i < nCount; i++)
	{
		pVector = pData->row(i);
		if(GVec::dotProductIgnoringUnknowns(pOutputMean, pVector + m_nInputVectorSize, pPrincipleComponent, m_nOutputVectorSize) >= 0)
		{
			for(int j = 0; j < m_nInputVectorSize; j++)
			{
				if(pVector[j] != UNKNOWN_REAL_VALUE)
				{
					pInputMean2[j] += pVector[j];
					pCounts2[j]++;
				}
			}
		}
		else
		{
			for(int j = 0; j < m_nInputVectorSize; j++)
			{
				if(pVector[j] != UNKNOWN_REAL_VALUE)
				{
					pInputMean1[j] += pVector[j];
					pCounts1[j]++;
				}
			}
		}
	}
	for(int j = 0; j < m_nInputVectorSize; j++)
	{
		if(pCounts1[j] == 0 || pCounts2[j] == 0)
			return new GMeanMarginsTreeLeafNode(m_nOutputVectorSize, pOutputMean);
		pInputMean1[j] /= pCounts1[j];
		pInputMean2[j] /= pCounts2[j];
	}

	// Compute the input center
	double* pCenter = new double[m_nInputVectorSize];
	GVec::copy(pCenter, pInputMean1, m_nInputVectorSize);
	GVec::add(pCenter, pInputMean2, m_nInputVectorSize);
	GVec::multiply(pCenter, .5, m_nInputVectorSize);

	// Compute the input normal
	double* pNormal = new double[m_nInputVectorSize];
	GVec::copy(pNormal, pInputMean2, m_nInputVectorSize);
	GVec::subtract(pNormal, pInputMean1, m_nInputVectorSize);
	if(GVec::squaredMagnitude(pNormal, m_nInputVectorSize) == 0)
		m_pRand->spherical(pNormal, m_nInputVectorSize);
	else
		GVec::normalize(pNormal, m_nInputVectorSize);

	// Make the interior node
	GMeanMarginsTreeInteriorNode* pNode = new GMeanMarginsTreeInteriorNode(pCenter, pNormal);

	// Divide the data
	GData other(pData->relation(), pData->heap());
	for(size_t i = pData->rows() - 1; i < pData->rows(); i--)
	{
		pVector = pData->row(i);
		if(pNode->Test(pVector, m_nInputVectorSize))
			other.takeRow(pData->releaseRow(i));
	}

	// If we couldn't separate anything, just return a leaf node
	if(other.rows() == 0 || pData->rows() == 0)
	{
		delete(pNode);
		return new GMeanMarginsTreeLeafNode(m_nOutputVectorSize, pOutputMean);
	}

	// Build the child nodes
	pNode->SetLeft(buildNode(pData, pBuf));
	pNode->SetRight(buildNode(&other, pBuf));

	return pNode;
}

// virtual
void GMeanMarginsTree::clear()
{
	delete(m_pRoot);
	m_pRoot = NULL;
	delete[] m_pEvalVector;
	m_pEvalVector = NULL;
}

#ifndef NO_TEST_CODE
// static
void GMeanMarginsTree::test()
{
	GRand prng(0);
	GMeanMarginsTree mm(&prng);
	GFilter tl(&mm, false);
	tl.setFeatureTransform(new GNominalToCat(), true);
	mm.basicTest(0.65, &prng);
}
#endif
