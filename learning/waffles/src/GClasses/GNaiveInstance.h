/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#ifndef __GNAIVEINSTANCE_H__
#define __GNAIVEINSTANCE_H__

#include "GLearner.h"

namespace GClasses {

class GNaiveInstanceAttr;
class GHeap;


/// This is an instance-based learner. Instead of finding the k-nearest
/// neighbors of a feature vector, it finds the k-nearst neighbors in each
/// dimension. That is, it finds n*k neighbors, considering each dimension
/// independently. It then combines the label from all of these neighbors
/// to make a prediction. Finding neighbors in this way makes it more robust to
/// high-dimensional datasets. It tends to perform worse than k-nn in low-dimensional space, and better
/// than k-nn in high-dimensional space. (It may be thought of as a cross
/// between a k-nn instance learner and a Naive Bayes learner. It only
/// supports continuous features and labels (so it is common to wrap it
/// in a Categorize filter which will convert nominal features to a categorical
/// distribution of continuous values).
class GNaiveInstance : public GIncrementalLearner
{
protected:
	sp_relation m_pRelation;
	int m_labelDims, m_featureDims;
	int m_nNeighbors;
	GNaiveInstanceAttr** m_pAttrs;
	double* m_pValueSums;
	double* m_pWeightSums;
	double* m_pSumBuffer;
	double* m_pSumOfSquares;
	GHeap* m_pHeap;

public:
	/// nNeighbors is the number of neighbors (in each dimension)
	/// that will contribute to the output value.
	GNaiveInstance(int nNeighbors);

	/// Deserializing constructor
	GNaiveInstance(GTwtNode* pNode);
	virtual ~GNaiveInstance();

#ifndef NO_TEST_CODE
	/// Performs unit tests for this class. Throws an exception if there is a failure.
	static void test();
#endif

	/// Save to a text-based format
	virtual GTwtNode* toTwt(GTwtDoc* pDoc);

	/// See the comment for GSupervisedLearner::featureDims
	virtual int featureDims();

	/// See the comment for GSupervisedLearner::labelDims
	virtual int labelDims();

	/// See the comment for GIncrementalLearner::enableIncrementalLearning
	virtual void enableIncrementalLearning(sp_relation& pRelation, int labelDims, double* pMins, double* pRanges);

	/// Incrementally train with a single instance
	virtual void trainIncremental(const double* pIn, const double* pOut);

	/// Train using all the samples in a collection
	virtual void train(GData* pData, int labelDims);

	/// See the comment for GIncrementalLearner::trainSparse
	virtual void trainSparse(GSparseMatrix* pData, int labelDims);

	/// See the comment for GSupervisedLearner::predictDistribution
	virtual void predictDistribution(const double* pIn, GPrediction* pOut);

	/// See the comment for GSupervisedLearner::predict
	virtual void predict(const double* pIn, double* pOut);

	/// See the comment for GSupervisedLearner::clear
	virtual void clear();

protected:
	void evalInput(int nInputDim, double dInput);
};

} // namespace GClasses

#endif // __GNAIVEINSTANCE_H__

