/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#ifndef __GDECISIONTREE_H__
#define __GDECISIONTREE_H__

#include "GLearner.h"
#include <vector>

namespace GClasses {

class GDecisionTreeNode;
class GRegressionTreeNode;
class GRand;
class GMeanMarginsTreeNode;
class GDecisionTreeLeafNode;


/// This is an efficient learning algorithm. It divides
/// on the attributes that reduce entropy the most, or alternatively
/// can make random divisions.
class GDecisionTree : public GSupervisedLearner
{
public:
	enum DivisionAlgorithm
	{
		MINIMIZE_ENTROPY,
		RANDOM,
	};

protected:
	sp_relation m_pRelation;
	GDecisionTreeNode* m_pRoot;
	DivisionAlgorithm m_eAlg;
	GRand* m_pRand;
	size_t m_leafThresh;
	size_t m_randomDraws;
	int m_featureDims, m_labelDims;

public:
	GDecisionTree(GRand* pRand);

	/// Makes a deep copy of another decision tree.  Also, if pInterestingNode
	/// is non-NULL, then ppOutInterestingNode will return the node that is
	/// a copy of pInterestingNode
	GDecisionTree(GDecisionTree* pThat, GDecisionTreeNode* pInterestingNode, GDecisionTreeNode** ppOutInterestingCopy);

	/// Loads the model from a text file in ".twt" format
	GDecisionTree(GTwtNode* pNode, GRand* pRand);

	virtual ~GDecisionTree();

#ifndef NO_TEST_CODE
	/// Performs unit tests for this class. Throws an exception if there is a failure.
	static void test();
#endif

	/// Saves to a text-based format
	virtual GTwtNode* toTwt(GTwtDoc* pDoc);

	/// See the comment for GSupervisedLearner::featureDims
	virtual int featureDims();

	/// See the comment for GSupervisedLearner::labelDims
	virtual int labelDims();

	/// Specifies for this decision tree to use random divisions (instead of
	/// divisions that reduce entropy). Random divisions make the algorithm
	/// train somewhat faster, and also increase model variance, so it is better
	/// suited for ensembles, but random divisions also make the decision tree
	/// vulnerable to problems with irrelevant attributes.
	void useRandomDivisions(size_t randomDraws = 1)
	{
		m_eAlg = RANDOM;
		m_randomDraws = randomDraws;
	}

	/// Sets the leaf threshold. When the number of samples is <= this value,
	/// it will no longer try to divide the data, but will create a leaf node.
	/// The default value is 1. For noisy data, a larger value may be advantageous.
	void setLeafThresh(size_t n) { m_leafThresh = n; }

	/// Frees the model
	virtual void clear();

	/// Trains this decision tree
	virtual void train(GData* pData, int labelDims);

	/// See the comment for GSupervisedLearner::predict
	void predict(const double* pIn, double* pOut);

	/// See the comment for GSupervisedLearner::predictDistribution
	virtual void predictDistribution(const double* pIn, GPrediction* pOut);

	/// Returns the number of nodes in this tree
	int treeSize();

	/// Prints an ascii representation of the decision tree to the specified stream.
	/// pRelation is an optional relation that can be supplied in order to provide
	/// better meta-data to make the print-out richer.
	void print(std::ostream& stream, GArffRelation* pRelation = NULL);

protected:
	/// Finds the leaf node that corresponds with the specified feature vector
	GDecisionTreeLeafNode* findLeaf(const double* pIn, int* pDepth);

	/// A recursive helper method used to construct the decision tree
	GDecisionTreeNode* buildBranch(GData* pData, std::vector<size_t>& attrPool, int nDepth, int tolerance);

	/// InfoGain is defined as the difference in entropy in the data
	/// before and after dividing it based on the specified attribute. For
	/// continuous attributes it uses the difference between the original
	/// variance and the sum of the variances of the two parts after
	/// dividing at the point the maximizes this value.
	double measureInfoGain(GData* pData, int nAttribute, double* pPivot);

	size_t pickDivision(GData* pData, double* pPivot, std::vector<size_t>& attrPool, int nDepth);
};



/// A GMeanMarginsTree is similar a DecisionTree, except it divides
/// as follows:
/// It finds the mean and principle component of the output vectors.
/// It divides all the vectors into two groups, one that has a
/// positive dot-product with the principle component (after subtracting
/// the mean) and one that has a negative dot-product with the
/// principle component (after subtracting the mean). Next it finds the
/// average input vector for each of the two groups. Then it finds
/// the mean and principle component of those two vectors. The dividing
/// criteria for this node is to subtract the mean and then see whether
/// the dot-product with the principle component is positive or negative
class GMeanMarginsTree : public GSupervisedLearner
{
protected:
	int m_labelDims;
	GMeanMarginsTreeNode* m_pRoot;
	int m_nInputVectorSize;
	int m_nOutputVectorSize;
	double* m_pEvalVector;
	GRand* m_pRand;
	sp_relation m_pRelation;

public:
	/// nOutputs specifies the number of output dimensions
	GMeanMarginsTree(GRand* pRand);

	/// Load from a text-based format
	GMeanMarginsTree(GTwtNode* pNode, GRand* pRand);

	virtual ~GMeanMarginsTree();

#ifndef NO_TEST_CODE
	static void test();
#endif

	/// Save to a text-based format
	virtual GTwtNode* toTwt(GTwtDoc* pDoc);

	/// See the comment for GSupervisedLearner::featureDims
	virtual int featureDims();

	/// See the comment for GSupervisedLearner::labelDims
	virtual int labelDims();

	/// Inductively build the tree
	virtual void train(GData* pData, int labelDims);

	/// See the comment for GSupervisedLearner::predict
	void predict(const double* pIn, double* pOut);

	/// See the comment for GSupervisedLearner::predictDistribution
	virtual void predictDistribution(const double* pIn, GPrediction* pOut);

	virtual void clear();

protected:
	GMeanMarginsTreeNode* buildNode(GData* pData, double* pBuf);
};


} // namespace GClasses

#endif // __GDECISIONTREE_H__
