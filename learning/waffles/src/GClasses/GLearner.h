/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#ifndef __GLEARNER_H__
#define __GLEARNER_H__

#include "GData.h"

namespace GClasses {

class GTwtDoc;
class GTwtNode;
class GDistribution;
class GCategoricalDistribution;
class GNormalDistribution;
class GUniformDistribution;
class GPlotWindow;
class GUnivariateDistribution;
class GIncrementalTransform;
class GTwoWayIncrementalTransform;
class GSparseMatrix;


/// This class is used to represent the predicted distribution made by a supervised learning algorithm.
/// (It is just a shallow wrapper around GDistribution.) It is used in conjunction with calls
/// to GSupervisedLearner::predictDistribution. The predicted distributions will be either
/// categorical distributions (for nominal values) or Normal distributions (for continuous values).
class GPrediction
{
protected:
	GUnivariateDistribution* m_pDistribution;

public:
	GPrediction() : m_pDistribution(NULL)
	{
	}

	~GPrediction();

	/// Converts an array of prediction objects to a vector of most-likely values.
	static void predictionArrayToVector(int nOutputCount, GPrediction* pOutputs, double* pVector);

	/// Converts an array of values to an array of predictions. There's not really
	/// enough information for this conversion, so it simply fabricates the variance
	/// and class-probability information as needed. Only the mean (for normal distributions)
	/// and the most-likely class (for categorical distributions) is reliable after this
	/// conversion.
	static void vectorToPredictionArray(GRelation* pRelation, int nOutputCount, double* pVector, GPrediction* pOutputs);

	/// Returns true if this wraps a normal distribution, false otherwise
	bool isContinuous();

	/// Returns the mode (most likely value). For the Normal distribution, this is the same as the mean.
	double mode();

	/// If the current distribution is not a categorical distribution, then it
	/// replaces it with a new categorical distribution. Then it returns the
	/// current (categorical) distribution.
	GCategoricalDistribution* makeCategorical();

	/// If the current distribution is not a normal distribution, then it
	/// replaces it with a new normal distribution. Then it returns the
	/// current (normal) distribution.
	GNormalDistribution* makeNormal();

	/// Returns the current distribution. Throws if it is not a categorical distribution
	GCategoricalDistribution* asCategorical();

	/// Returns the current distribution. Throws if it is not a normal distribution
	GNormalDistribution* asNormal();
};




// nRep and nFold are zero-indexed
typedef void (*RepValidateCallback)(void* pThis, int nRep, int nFold, int labelDims, double* pFoldResults);


/// This is the base class of supervised learning algorithms (that may or may not
/// have an internal model allowing them to generalize rows that were not available
/// at training time). Note that the literature typically refers to supervised learning
/// algorithms that can't generalize (because they lack an internal hypothesis model)
/// as "Semi-supervised". (You cannot generalize with a semi-supervised algorithm--you have to
/// train again with the new rows.)
class GTransducer
{
public:
	GTransducer();
	GTransducer(GTwtNode* pLearner);
	virtual ~GTransducer();

	/// Returns false because semi-supervised learners have no internal
	/// model, so they can't evaluate previously unseen rows.
	virtual bool canGeneralize() { return false; }

	/// Returns false because semi-supervised learners cannot be trained incrementally.
	virtual bool canTrainIncrementally() { return false; }

	/// Trains with pDataLabeled, and then predicts labels for all the
	/// rows in pDataUnlabeled. (Works with model-free algorithms that
	/// cannot be trained.) The rows in pDataUnlabeled are expected
	/// to have space allocated for labels, but the values will be overwritten
	/// with the predicted values.
	virtual void transduce(GData* pDataLabeled, GData* pDataUnlabeled, int labelDims) = 0;

	/// Trains and tests this learner. pOutResults should have
	/// an element for each label dim.
	virtual void trainAndTest(GData* pTrainingSet, GData* pTestSet, int labelDims, double* pOutResults);

	/// Perform n-fold cross validation on pData. Uses trainAndTest
	/// for each fold. pCB is an optional callback method for reporting
	/// intermediate stats. It can be NULL if you don't want intermediate reporting.
	/// nRep is just the rep number that will be passed to the callback.
	/// pThis is just a pointer that will be passed to the callback for you
	/// to use however you want. It doesn't affect this method.
	/// The results of each fold is returned in a dataset.
	GData* crossValidate(GData* pData, int nFolds, int labelDims, RepValidateCallback pCB = NULL, int nRep = 0, void* pThis = NULL);

	/// Perform cross validation "nReps" times and return the
	/// average score. (5 reps with 2 folds is preferred over 10-fold cross
	/// validation because it yields less type 1 error.)
	/// pCB is an optional callback method for reporting intermediate stats
	/// It can be NULL if you don't want intermediate reporting.
	/// pThis is just a pointer that will be passed to the callback for you
	/// to use however you want. It doesn't affect this method.
	/// The results of each fold is returned in a dataset.
	GData* repValidate(GData* pData, int nReps, int nFolds, int labelDims, GRand* pRand, RepValidateCallback pCB = NULL, void* pThis = NULL);

	/// This performs two-fold cross-validation on a shuffled
	/// non-uniform split of the data, and returns an error value that
	/// represents the results of all labels combined. (This is
	/// for heuristic optimization, not for reporting accuracy.)
	double heuristicValidate(GData* pData, int labelDims, GRand* pRand);

protected:
	/// Child classes should use this in their implementation of toTwt
	GTwtNode* baseTwtNode(GTwtDoc* pDoc, const char* szClassName);
};


/// This is the base class of algorithms that learn with supervision and
/// have an internal hypothesis model that allows them to generalize
/// rows that were not available at training time.
class GSupervisedLearner : public GTransducer
{
public:
	GSupervisedLearner();
	GSupervisedLearner(GTwtNode* pLearner);
	virtual ~GSupervisedLearner();

	/// Save to a text-based format. Implementations of this method should
	/// use MakeBaseTwtNode
	virtual GTwtNode* toTwt(GTwtDoc* pDoc) = 0;

	/// Returns true because fully supervised learners have an internal
	/// model that allows them to generalize previously unseen rows.
	virtual bool canGeneralize() { return true; }

	/// Returns the number of feature dims. Throws an exception if the
	/// model has not yet been trained. (In the case of incremental
	/// learners, it should work if enableIncrementalLearning has been
	/// called.)
	virtual int featureDims() = 0;

	/// Returns the number of label dims. Throws an exception if the
	/// model has not yet been trained. (In the case of incremental
	/// learners, it should work if enableIncrementalLearning has been
	/// called.)
	virtual int labelDims() = 0;

	/// Trains the model with pData. (Any previous training will be
	/// discarded.) The last "labelDims" columns in the data represent the
	/// label. All of the preceding columns are the feature vectors.
	virtual void train(GData* pData, int labelDims) = 0;

	/// Evaluate pIn and compute a prediction for pOut. pOut is expected
	/// to point to an array of GPrediction objects which have already been
	/// allocated. There should be labelDims() elements in this array.
	/// Even though the prediction is a probability distribution,
	/// this distribution is not calibrated. For continuous
	/// label dimensions, the mean should be a good prediction, but the variance
	/// may be a very rough estimate, or even a bogus value, depending on
	/// the algorithm. For discrete label dimensions, the mode is the
	/// most-likely prediction, but the confidence values of
	/// other categories may be determined by weak heuristics, depending
	/// on the algorithm.
	virtual void predictDistribution(const double* pIn, GPrediction* pOut) = 0;

	/// Evaluate pIn and compute a prediction for pOut.
	virtual void predict(const double* pIn, double* pOut) = 0;

	/// Discards all training for the purpose of freeing memory.
	/// If you call this method, you must train before making any predictions.
	/// No settings or options are discarded.
	virtual void clear() = 0;

	/// Computes predictive accuracy for nominal label dimensions, and
	/// the mean squared error for continuous outputs.
	/// pOutResults should have m_labelDims dimensions.
	void accuracy(GData* pData, double* pOutResults);

	/// label specifies which output to measure. (It should be 0 if there is only one label dimension.)
	/// The measurement will be performed "nReps" times and results averaged together
	/// nPrecisionSize specifies the number of points at which the function is sampled
	/// pOutPrecision should be an array big enough to hold nPrecisionSize elements for every possible
	/// label value. (If the attribute is continuous, it should just be big enough to hold nPrecisionSize elements.)
	/// If bLocal is true, it computes the local precision instead of the global precision.
	void precisionRecall(double* pOutPrecision, int nPrecisionSize, GData* pData, int labelDims, int label, int nReps, GRand* pRand);

	/// See GSemiSupervisedLearner::transduce
	virtual void transduce(GData* pDataLabeled, GData* pDataUnlabeled, int labelDims);

	/// Trains and tests this learner
	virtual void trainAndTest(GData* pTrainingSet, GData* pTestSet, int labelDims, double* pOutResults);

#ifndef NO_TEST_CODE
	/// This is a helper method used by the unit tests of several model learners
	void basicTest(double minAccuracy, GRand* pRand, double deviation = 1e-6);
#endif
protected:
	/// This is a helper method used by precisionRecall.
	int precisionRecallNominal(GPrediction* pOutput, double* pFunc, GData* pTrain, GData* pTest, int labelDims, int label, int value);

	/// This is a helper method used by precisionRecall.
	int precisionRecallContinuous(GPrediction* pOutput, double* pFunc, GData* pTrain, GData* pTest, int labelDims, int label);
};




/// This is the base class of supervised learning algorithms that
/// can learn one row at a time.
class GIncrementalLearner : public GSupervisedLearner
{
public:
	GIncrementalLearner()
	: GSupervisedLearner()
	{
	}

	GIncrementalLearner(GTwtNode* pLearner)
	 : GSupervisedLearner(pLearner)
	{
	}

	virtual ~GIncrementalLearner()
	{
	}

	/// Only the GFilter class should return true to this method
	virtual bool isFilter() { return false; }

	/// Returns true
	virtual bool canTrainIncrementally() { return true; }

	/// You must call this method before you call trainIncremental
	/// pMins specifies the minimum values. pRanges specifies (max - min) values.
	/// (Only the min and range values for continuous attributes will actually be used,
	/// so it is safe to pass in NULL for both arrays if all attributes are nominal.)
	virtual void enableIncrementalLearning(sp_relation& pRelation, int labelDims, double* pMins, double* pRanges) = 0;

	/// Pass a single input row and the corresponding label to
	/// incrementally train this model
	virtual void trainIncremental(const double* pIn, const double* pOut) = 0;

	/// Train using a sparse matrix. (Typically, implementations of this
	/// method will iterate over the rows in pData, and for each row it
	/// will convert the sparse row to a full row, call trainIncremental,
	/// and then discard the full row.)
	virtual void trainSparse(GSparseMatrix* pData, int labelDims) = 0;
};



/// This class is for loading various learning algorithms from a Twt node. When any
/// learning algorithm is saved, it calls MakeBaseTwtNode, which creates (among other
/// things) a field named "class" which specifies the class name of the algorithm.
/// This class contains methods that will recognize any of the classes in this library
/// and load them. If it doesn't recognize a class, it will either return NULL or throw
/// and exception, depending on the flags you pass to the constructor.
/// Obviously this loader won't recognize any classes that you make. Therefore, you should
/// overload the corresponding method in this class with a new method that will first
/// recognize and load your classes, and then call these methods to handle other types.
class GLearnerLoader
{
protected:
	bool m_throwIfClassNotFound;

public:
	GLearnerLoader(bool throwIfClassNotFound = true) { m_throwIfClassNotFound = throwIfClassNotFound; }
	virtual ~GLearnerLoader() {}

	virtual GIncrementalTransform* loadIncrementalTransform(GTwtNode* pNode, GRand* pRand);
	virtual GTwoWayIncrementalTransform* loadTwoWayIncrementalTransform(GTwtNode* pNode, GRand* pRand);
	virtual GSupervisedLearner* loadModeler(GTwtNode* pNode, GRand* pRand);
	virtual GIncrementalLearner* loadIncrementalLearner(GTwtNode* pNode, GRand* pRand);
};



/// Always outputs the label mean (for continuous labels) and the most common
/// class (for nominal labels).
class GBaselineLearner : public GSupervisedLearner
{
protected:
	int m_featureDims, m_labelDims;
	sp_relation m_pRelation;
	double* m_pLabels;
	GRand* m_pRand;

public:
	GBaselineLearner(GRand* pRand);

	/// Load from a text-based format
	GBaselineLearner(GTwtNode* pNode, GRand* pRand);

	virtual ~GBaselineLearner();

#ifndef NO_TEST_CODE
	static void test();
#endif

	/// Save to a text-based format
	virtual GTwtNode* toTwt(GTwtDoc* pDoc);

	/// See the comment for GSupervisedLearner::labelDims
	virtual int featureDims();

	/// See the comment for GSupervisedLearner::labelDims
	virtual int labelDims();

	/// See the comment for GSupervisedLearner::clear
	virtual void clear();

	/// Trains this model
	virtual void train(GData* pData, int labelDims);

	/// Predicts a distribution
	virtual void predictDistribution(const double* pIn, GPrediction* pOut);

	/// Makes a prediction (which is always the same)
	virtual void predict(const double* pIn, double* pOut);
};


/// This is an implementation of the identity function. It might be
/// useful, for example, as the observation function in a GRecurrentModel
/// if you want to create a Jordan network.
class GIdentityFunction : public GSupervisedLearner
{
protected:
	int m_labelDims;
	int m_featureDims;

public:
	GIdentityFunction();

	/// Load from a text-based format
	GIdentityFunction(GTwtNode* pNode);

	virtual ~GIdentityFunction();

	/// Save to a text-based format
	virtual GTwtNode* toTwt(GTwtDoc* pDoc);

	/// See the comment for GSupervisedLearner::labelDims
	virtual int featureDims();

	/// See the comment for GSupervisedLearner::labelDims
	virtual int labelDims();

	/// See the comment for GSupervisedLearner::clear
	virtual void clear();

	/// This method is a no-op
	virtual void train(GData* pData, int labelDims);

	/// Not implemented yet
	virtual void predictDistribution(const double* pIn, GPrediction* pOut);

	/// Copies pIn to pOut. If there are more feature dims than
	/// label dims, then the extras are clipped. If there are
	/// more label dims than feature dims, then the extra values
	/// are predicted to be 0.
	virtual void predict(const double* pIn, double* pOut);
};


} // namespace GClasses

#endif // __GLEARNER_H__

