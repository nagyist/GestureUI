/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#ifndef __GTRANSFORM_H__
#define __GTRANSFORM_H__

#include "GLearner.h"

namespace GClasses {

/// This is the base class of algorithms that transform data without supervision
class GTransform
{
public:
	GTransform();
	GTransform(GTwtNode* pNode);
	virtual ~GTransform();

	/// Applies the transformation to pIn and returns the results. For
	/// transformations with an internal model (including all transforms
	/// that inherit from GIncrementalTransform), this is equivalent to calling
	/// train, and then calling transformBatch.
	virtual GData* doit(GData* pIn) = 0;

protected:
	/// Child classes should use this in their implementation of toTwt
	GTwtNode* baseTwtNode(GTwtDoc* pDoc, const char* szClassName);
};




/// This wraps two unsupervised learners to make a single unsupervised learner
/// that performs both transforms sequentially
class GTransformChainer : public GTransform
{
protected:
	GTransform* m_pFirst;
	GTransform* m_pSecond;

public:
	GTransformChainer(GTransform* pFirst, GTransform* pSecond);
	virtual ~GTransformChainer();

	/// Transform the data using both transforms
	virtual GData* doit(GData* pIn);
};





/// This is the base class of algorithms that can transform data
/// one row at a time without supervision.
class GIncrementalTransform : public GTransform
{
protected:
	sp_relation m_pRelationBefore;
	sp_relation m_pRelationAfter;
	double* m_pAfterMins;
	double* m_pAfterRanges;

public:
	GIncrementalTransform() : GTransform(), m_pAfterMins(NULL), m_pAfterRanges(NULL) {}
	GIncrementalTransform(GTwtNode* pNode) : GTransform(pNode), m_pAfterMins(NULL), m_pAfterRanges(NULL) {}
	virtual ~GIncrementalTransform() { delete[] m_pAfterMins; }

	/// Save to a text-based format
	virtual GTwtNode* toTwt(GTwtDoc* pDoc) = 0;

	/// sets m_pRelationBefore and m_pRelationAfter, and trains the transform.
	virtual void train(GData* pData) = 0;

	/// Prepares the transform to be used with incremental training
	virtual void enableIncrementalTraining(sp_relation& pRelation, double* pMins, double* pRanges) = 0;

	/// train must be called before this method is used
	sp_relation& before() { return m_pRelationBefore; }

	/// train must be called before this method is used
	sp_relation& after() { return m_pRelationAfter; }

	/// enableIncrementalTraining must be called before this method is used. It returns a
	/// vector of minimum values for the data after the transform has been applied.
	const double* afterMins() { return m_pAfterMins; }

	/// enableIncrementalTraining must be called before this method is used. It returns a
	/// vector of range values for the data after the transform has been applied.
	const double* afterRanges() { return m_pAfterRanges; }

	/// pIn is the source row. pOut is a buffer that will hold the transformed row.
	/// train must be called before this method is used
	virtual void transform(const double* pIn, double* pOut) = 0;

	/// This calls Train with pIn, then transforms pIn and returns the results.
	virtual GData* doit(GData* pIn);

	/// This assumes that train has already been called, and transforms all the rows in pIn.
	virtual GData* transformBatch(GData* pIn);
};




/// This is the base class of algorithms that can transform data
/// one row at a time without supervision, and can (un)transform
/// a row back to its original form if necessary.
class GTwoWayIncrementalTransform : public GIncrementalTransform
{
public:
	GTwoWayIncrementalTransform() : GIncrementalTransform() {}
	GTwoWayIncrementalTransform(GTwtNode* pNode) : GIncrementalTransform(pNode) {}
	virtual ~GTwoWayIncrementalTransform() {}

	/// pIn is a previously transformed row, and pOut is a buffer that will hold the untransformed row.
	/// train must be called before this method is used
	virtual void untransform(const double* pIn, double* pOut) = 0;

	/// This assumes train was previously called, and untransforms all the rows in pIn and returns the results.
	virtual GData* untransformBatch(GData* pIn);
};




/// This class enables functionality similar to filters in Weka.
/// It wraps a modeler and transforms data before it is passed to
/// the modeler, and/or after the modeler passes data back. For example,
/// suppose that you wish to use GNaiveBayes (which only supports
/// discrete attributes) with a dataset that contains continuous attributes.
/// You could use this class to wrap GDiscretize around GNaiveBayes to
/// create a naive bayes that can operate on both continuous and discrete
/// data. As another example, if you have a modeler that only supports real
/// attributes, you could wrap it with GNominalToCat to create a modeler
/// that can operate on nominal data too. Or, if your modeler expects values
/// within a certain range, you could use GNormalize to ensure
/// that the modeler receives values within that range while leaving
/// the caller free to pass in values with whatever range it prefers.
class GFilter : public GIncrementalLearner
{
protected:
	int m_beforeLabelDims, m_afterLabelDims;
	sp_relation m_pBeforeRelation;
	sp_relation m_pAfterRelation;
	GIncrementalTransform* m_pFeatureTransform;
	GTwoWayIncrementalTransform* m_pLabelTransform;
	GSupervisedLearner* m_pModeler;
	bool m_ownFeatureTransform, m_ownLabelTransform, m_ownModeler;
	double* m_pScratchVector;
	double* m_pScratchVector2;
	GPrediction* m_pScratchLabels;

public:
	/// if own is true, this will delete pModeler when it's done with it.
	GFilter(GSupervisedLearner* pModeler, bool own);

	/// Load from a text-based format
	GFilter(GTwtNode* pNode, GRand* pRand, GLearnerLoader* pLoader);

	virtual ~GFilter();

	/// Returns true (because this is the GFilter class)
	virtual bool isFilter() { return true; }

	/// Returns the modeler
	GSupervisedLearner* modeler() { return m_pModeler; }

	/// Save to a text-based format.
	virtual GTwtNode* toTwt(GTwtDoc* pDoc);

	/// Returns the number of feature dims before the transformation.
	/// See the comment for GSupervisedLearner::featureDims
	virtual int featureDims();

	/// Returns the number of label dims before the transformation.
	/// See the comment for GSupervisedLearner::labelDims
	virtual int labelDims();

	/// See the comment for GSupervisedLearner::clear
	virtual void clear();

	/// Set the feature transform. (You don't need to call this if you don't want to transform the features)
	/// If takeOwnership is true, then this object will take care of deleting pFeatureTransform.
	void setFeatureTransform(GIncrementalTransform* pFeatureTransform, bool takeOwnership);

	/// Set the label transform. (You don't need to call this if you don't want to transform the labels)
	/// If takeOwnership is true, then this object will take care of deleting pLabelTransform.
	void setLabelTransform(GTwoWayIncrementalTransform* pLabelTransform, bool takeOwnership);

	/// Transforms both the features and the labels, and then trains the wrapped modeler
	/// with the transformed data.
	virtual void train(GData* pData, int labelDims);

	/// Transforms the features, evaluates with the wrapped modeler on the transformed
	/// features, and then un-transforms the labels to return a prediction in the original
	/// form.
	virtual void predictDistribution(const double* pIn, GPrediction* pOut);

	/// see GSupervisedLearner::predict
	virtual void predict(const double* pIn, double* pOut);

	/// Calls the same method on the wrapped model. You are responsible to
	/// ensure that the transforms are correctly prepped to handle incremental learning.
	virtual void enableIncrementalLearning(sp_relation& pRelation, int labelDims, double* pMins, double* pRanges);

	/// see GIncrementalLearner::trainIncremental
	virtual void trainIncremental(const double* pIn, const double* pOut);

	/// see GIncrementalLearner::trainSparse
	virtual void trainSparse(GSparseMatrix* pData, int labelDims);

protected:
	void makeScratchBuffers();
};




/// Principal Component Analysis. (Computes the principal components about
/// the mean of the data when you call train. The transformed (reduced-dimensional)
/// data will have a mean about the origin.)
class GPCA : public GIncrementalTransform
{
protected:
	int m_targetDims;
	GData* m_pBasisVectors;
	double* m_pEigVals;
	GRand* m_pRand;

public:
	GPCA(int targetDims, GRand* pRand);

	/// Load from a text-based format
	GPCA(GTwtNode* pNode, GRand* pRand);

	virtual ~GPCA();

	/// Save to a text-based format
	virtual GTwtNode* toTwt(GTwtDoc* pDoc);

	/// Specify to compute the eigenvalues during training. This
	/// method must be called before train is called.
	void computeEigVals();

	/// Returns the eigenvalues. Returns NULL if computeEigVals was not called.
	double* eigVals() { return m_pEigVals; }

	/// Returns the number of principal components
	int targetDims() { return m_targetDims; }

	/// Computes a (lossy) high-dimensional point that corresponds with the
	/// specified low-dimensional coordinates.
	void reverse(const double* pIn, double* pOut);

	/// Returns the mean of the data used to train this transform
	double* mean() { return m_pBasisVectors->row(0); }

	/// Returns the i'th principal component vector
	double* basis(int i) { return m_pBasisVectors->row(i + 1); }

	/// See the comment for GIncrementalTransform::train
	virtual void train(GData* pData);
	
	/// See the comment for GIncrementalTransform::enableIncrementalTraining
	virtual void enableIncrementalTraining(sp_relation& pRelation, double* pMins, double* pRanges);
	
	/// See the comment for GIncrementalTransform::transform
	virtual void transform(const double* pIn, double* pOut);
};


/// Principle Component Analysis without the projection. It only rotates
/// axes to align with the first few principal components.
class GPCARotateOnly
{
public:
	/// This rotates the data to align the first nComponents axes with the same
	/// number of principle components.
	static GData* transform(int nDims, int nOutputs, GData* pData, int nComponents, GRand* pRand);

#ifndef NO_TEST_CODE
	/// Performs unit tests for this class. Throws an exception if there is a failure.
	static void test();
#endif // !NO_TEST_CODE
};



/// Just generates Gaussian noise
class GNoiseGenerator : public GIncrementalTransform
{
protected:
	GRand* m_pRand;
	double m_mean, m_deviation;

public:
	GNoiseGenerator(GRand* pRand);

	/// Load from a text-based format
	GNoiseGenerator(GTwtNode* pNode, GRand* pRand);

	virtual ~GNoiseGenerator();

	/// Save to a text-based format
	virtual GTwtNode* toTwt(GTwtDoc* pDoc);

	/// See the comment for GIncrementalTransform::train
	virtual void train(GData* pData);

	/// See the comment for GIncrementalTransform::enableIncrementalTraining
	virtual void enableIncrementalTraining(sp_relation& pRelation, double* pMins, double* pRanges);
	virtual sp_relation& relationAfter() { return m_pRelationBefore; }
	
	/// See the comment for GIncrementalTransform::transform
	virtual void transform(const double* pIn, double* pOut);

	void setMeanAndDeviation(double m, double d) { m_mean = m; m_deviation = d; }
};



/// Generates data by computing the product of each pair of
/// attributes. This is useful for augmenting data.
class GPairProduct : public GIncrementalTransform
{
protected:
	int m_maxDims;

public:
	GPairProduct(int nMaxDims);

	/// Load from a text-based format
	GPairProduct(GTwtNode* pNode);

	virtual ~GPairProduct();

	/// Save to a text-based format
	virtual GTwtNode* toTwt(GTwtDoc* pDoc);

	/// See the comment for GIncrementalTransform::train
	virtual void train(GData* pData);
	
	/// See the comment for GIncrementalTransform::enableIncrementalTraining
	virtual void enableIncrementalTraining(sp_relation& pRelation, double* pMins, double* pRanges);
	
	/// See the comment for GIncrementalTransform::transform
	virtual void transform(const double* pIn, double* pOut);
};



/// Generates subsets of data that contain only the most relevant features for predicting the labels.
/// The train method of this class produces a ranked ordering of the feature attributes by training
/// a single-layer neural network, and deselecting the weakest attribute until all attributes have been
/// deselected. The transform method uses only the highest-ranked attributes.
class GAttributeSelector : public GIncrementalTransform
{
protected:
	int m_labelDims;
	int m_targetFeatures;
	std::vector<int> m_ranks;
	GRand* m_pRand;

public:
	GAttributeSelector(int labelDims, int targetFeatures, GRand* pRand) : GIncrementalTransform(), m_labelDims(labelDims), m_targetFeatures(targetFeatures), m_pRand(pRand)
	{
	}

	GAttributeSelector(GTwtNode* pNode, GRand* pRand);

	virtual ~GAttributeSelector()
	{
	}

#ifndef NO_TEST_CODE
	static void test();
#endif

	virtual GTwtNode* toTwt(GTwtDoc* pDoc);
	
	/// See the comment for GIncrementalTransform::train
	virtual void train(GData* pData);
	
	/// See the comment for GIncrementalTransform::enableIncrementalTraining
	virtual void enableIncrementalTraining(sp_relation& pRelation, double* pMins, double* pRanges);
	
	/// See the comment for GIncrementalTransform::transform
	virtual void transform(const double* pIn, double* pOut);

	/// Specifies the number of features to select
	void setTargetFeatures(int n);

	/// Returns a list of attributes in ranked-order. Most important attributes are first. Weakest attributes are last.
	/// (The results are undefined until after train is called.)
	std::vector<int>& ranks() { return m_ranks; }
};



/// This is sort-of the opposite of discretize. It converts each nominal attribute to a categorical
/// distribution by representing each value using the corresponding row of the identity matrix. For
/// example, if a certain nominal attribute has 4 possible values, then a value of 3 would be encoded
/// as the vector 0 0 1 0. When predictions are converted back to nominal values, the mode of the
/// categorical distribution is used as the predicted value. (This is similar to Weka's
/// NominalToBinaryFilter.)
class GNominalToCat : public GTwoWayIncrementalTransform
{
protected:
	int m_valueCap;
	GRand* m_pRand;
	std::vector<int> m_ranks;

public:
	GNominalToCat(int valueCap = 12);

	/// Load from a text-based format
	GNominalToCat(GTwtNode* pNode);

	virtual ~GNominalToCat();

	/// Save to a text-based format
	virtual GTwtNode* toTwt(GTwtDoc* pDoc);

	/// See the comment for GIncrementalTransform::train
	virtual void train(GData* pData);
	
	/// See the comment for GIncrementalTransform::enableIncrementalTraining
	virtual void enableIncrementalTraining(sp_relation& pRelation, double* pMins, double* pRanges);
	
	/// See the comment for GIncrementalTransform::transform
	virtual void transform(const double* pIn, double* pOut);
	
	/// See the comment for GTwoWayIncrementalTransform::untransform
	virtual void untransform(const double* pIn, double* pOut);

	/// Makes a mapping from the post-transform attribute indexes to the pre-transform attribute indexes
	void reverseAttrMap(std::vector<int>& rmap);

protected:
	void init(sp_relation& pRelationBefore);
};



/// This transform scales and shifts continuous values
/// to make them fall within a specified range.
class GNormalize : public GTwoWayIncrementalTransform
{
protected:
	double m_min, m_max;
	double* m_pMins;
	double* m_pRanges;

public:
	/// min and max specify the target range. (The input domain is determined
	/// automatically when train is called.)
	GNormalize(double min = 0.0, double max = 1.0);

	/// Load from a text-based format
	GNormalize(GTwtNode* pNode);

	virtual ~GNormalize();

	/// Save to a text-based format
	virtual GTwtNode* toTwt(GTwtDoc* pDoc);

	/// See the comment for GIncrementalTransform::train
	virtual void train(GData* pData);
	
	/// See the comment for GIncrementalTransform::enableIncrementalTraining
	virtual void enableIncrementalTraining(sp_relation& pRelation, double* pMins, double* pRanges);
	
	/// See the comment for GIncrementalTransform::transform
	virtual void transform(const double* pIn, double* pOut);
	
	/// See the comment for GTwoWayIncrementalTransform::untransform
	virtual void untransform(const double* pIn, double* pOut);

	void setMinsAndRanges(sp_relation& pRel, const double* pMins, const double* pRanges);
};



/// This transform uses buckets to convert continuous data into discrete data.
/// It is common to use GFilter to combine this with your favorite modeler
/// (which only supports discrete values) to create a modeler that can also support
/// continuous values as well.
class GDiscretize : public GTwoWayIncrementalTransform
{
protected:
	int m_bucketsIn, m_bucketsOut;
	double* m_pMins;
	double* m_pRanges;

public:
	/// if buckets is less than 0, then it will use the floor of the square root of the number of rows in the data
	GDiscretize(int buckets = -1);

	/// Load from a text-based format
	GDiscretize(GTwtNode* pNode);

	virtual ~GDiscretize();

	/// Save to a text-based format
	virtual GTwtNode* toTwt(GTwtDoc* pDoc);

	/// See the comment for GIncrementalTransform::train
	virtual void train(GData* pData);
	
	/// See the comment for GIncrementalTransform::enableIncrementalTraining
	virtual void enableIncrementalTraining(sp_relation& pRelation, double* pMins, double* pRanges);
	
	/// See the comment for GIncrementalTransform::transform
	virtual void transform(const double* pIn, double* pOut);
	
	/// See the comment for GTwoWayIncrementalTransform::untransform
	virtual void untransform(const double* pIn, double* pOut);
};


} // namespace GClasses

#endif // __GTRANSFORM_H__

