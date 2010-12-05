/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#ifndef __GENSEMBLE_H__
#define __GENSEMBLE_H__

#include "GLearner.h"
#include <vector>
#include <exception>

namespace GClasses {

class GRelation;
class GRand;


typedef void (*EnsembleProgressCallback)(void* pThis, int i, int n);

/// BAG stands for bootstrap aggregator. It represents an ensemble
/// of voting modelers. Each model is trained with a slightly different
/// training set, which is produced by drawing randomly from the original
/// training set with replacement until we have a new training set of
/// the same size. Each model is given equal weight in the vote.
class GBag : public GSupervisedLearner
{
protected:
	int m_labelDims;
	std::vector<GSupervisedLearner*> m_models;
	int m_nAccumulatorDims;
	double* m_pAccumulator;
	sp_relation m_pRelation;
	GRand* m_pRand;
	EnsembleProgressCallback m_pCB;
	void* m_pThis;

public:
	/// nInitialSize tells it how fast to grow the dynamic array that holds the
	/// models. It's not really important to get it right, just guess how many
	/// models will go in the ensemble.
	GBag(GRand* pRand);

	/// Load from a text-based format
	GBag(GTwtNode* pNode, GRand* pRand, GLearnerLoader* pLoader);

	virtual ~GBag();

#ifndef NO_TEST_CODE
	static void test();
#endif

	/// Saves to a text-based format
	virtual GTwtNode* toTwt(GTwtDoc* pDoc);

	/// See the comment for GSupervisedLearner::featureDims
	virtual int featureDims();

	/// See the comment for GSupervisedLearner::labelDims
	virtual int labelDims();

	/// See the comment for GSupervisedLearner::clear
	virtual void clear();

	/// Removes and deletes all the learners
	void flush();

	/// Adds a learner to the bag. This takes ownership of pLearner (so
	/// it will delete it when it's done with it)
	void addLearner(GSupervisedLearner* pLearner);

	/// Trains all the learning algorithms with the complete data set. Note
	/// that it's common to train each algorithm with a different portion
	/// of the data set. That's not what this method does. todo: write one
	/// that does that.
	virtual void train(GData* pData, int labelDims);

	/// See the comment for GSupervisedLearner::predictDistribution
	virtual void predictDistribution(const double* pIn, GPrediction* pOut);

	/// See the comment for GSupervisedLearner::predict
	virtual void predict(const double* pIn, double* pOut);

	/// If you want to be notified when another instance begins training, you can set this callback
	void setProgressCallback(EnsembleProgressCallback pCB, void* pThis)
	{
		m_pCB = pCB;
		m_pThis = pThis;
	}

protected:
	void accumulate(const double* pOut);
	void tally(int nCount, GPrediction* pOut);
	void tally(int nCount, double* pOut);
};



/// When Train is called, this performs cross-validation on the training
/// set to determine which learner is the best. It then trains that learner
/// with the entire training set.
class GBucket : public GSupervisedLearner
{
protected:
	int m_featureDims, m_labelDims;
	int m_nBestLearner;
	std::vector<GSupervisedLearner*> m_models;
	GRand* m_pRand;

public:
	/// nInitialSize tells it how fast to grow the dynamic array that holds the
	/// models. It's not really important to get it right, just guess how many
	/// models will go in the ensemble.
	GBucket(GRand* pRand);

	/// Load from a text-based format
	GBucket(GTwtNode* pNode, GRand* pRand, GLearnerLoader* pLoader);

	virtual ~GBucket();

#ifndef NO_TEST_CODE
	static void test();
#endif

	/// Saves to a text-based format
	virtual GTwtNode* toTwt(GTwtDoc* pDoc);

	/// See the comment for GSupervisedLearner::featureDims
	virtual int featureDims();

	/// See the comment for GSupervisedLearner::labelDims
	virtual int labelDims();

	/// See the comment for GSupervisedLearner::clear
	virtual void clear();

	/// Removes and deletes all the learners
	void flush();

	/// Adds a modeler to the list. This takes ownership of pLearner (so
	/// it will delete it when it's done with it)
	void addLearner(GSupervisedLearner* pLearner);

	/// Returns the modeler that did the best with the training set. It is
	/// your responsibility to delete the modeler this returns. Throws if
	/// you haven't trained yet.
	GSupervisedLearner* releaseBestModeler();

	/// Uses cross-validation to select the best learner, and then trains
	/// it with the entire training set.
	virtual void train(GData* pData, int labelDims);

	/// See the comment for GSupervisedLearner::predictDistribution
	virtual void predictDistribution(const double* pIn, GPrediction* pOut);

	/// See the comment for GSupervisedLearner::predict
	virtual void predict(const double* pIn, double* pOut);

	/// If one of the algorithms throws during training,
	/// it will catch it and call this no-op method. Overload
	/// this method if you don't want to ignore exceptions.
	virtual void onError(std::exception& e);
};


} // namespace GClasses

#endif // __GENSEMBLE_H__
