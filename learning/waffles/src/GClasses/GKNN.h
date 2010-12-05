/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#ifndef __GKNN_H__
#define __GKNN_H__

#include "GLearner.h"

namespace GClasses {

class GNeighborFinderGeneralizing;
class GRand;
class GAutoScaleInstance;
class GAgglomerativeClusterer;
class GFuzzyKarnaughHyperRect;
class GKnnScaleFactorCritic;
class GOptimizer;
class GRowDistanceScaled;


/// The k-Nearest Neighbor learning algorithm
class GKNN : public GIncrementalLearner
{
public:
	enum InterpolationMethod
	{
		Linear,
		Mean,
		Learner,
	};

protected:
	/// Settings
	GRand* m_pRand;
	GData* m_pInstances;
	sp_relation m_pRelation;
	int m_nNeighbors, m_featureDims, m_labelDims;
	InterpolationMethod m_eInterpolationMethod;
	GSupervisedLearner* m_pLearner;
	bool m_bOwnLearner;
	double m_dElbowRoom;

	/// Scale Factor Optimization
	bool m_optimizeScaleFactors;
	GRowDistanceScaled* m_pDistanceMetric;
	GKnnScaleFactorCritic* m_pCritic;
	GOptimizer* m_pScaleFactorOptimizer;

	/// Working Buffers
	size_t* m_pEvalNeighbors;
	double* m_pEvalDistances;
	double* m_pValueCounts;

	/// Neighbor Finding
	GNeighborFinderGeneralizing* m_pNeighborFinder; // used for evaluation
	GNeighborFinderGeneralizing* m_pNeighborFinder2; // used for incremental training

public:
	/// nOutputCout specifies the number of attribute dimensions.
	/// nNeighbors specifies the number of neighbors to use for evaluation.
	GKNN(int nNeighbors, GRand* pRand);

	/// Load from a text-based format
	GKNN(GTwtNode* pNode, GRand* pRand);

	virtual ~GKNN();

#ifndef NO_TEST_CODE
	/// Performs unit tests for this class. Throws an exception if there is a failure.
	static void test();
#endif

	/// See the comment for GSupervisedLearner::labelDims
	virtual int featureDims();

	/// See the comment for GSupervisedLearner::labelDims
	virtual int labelDims();

	/// Save to a text-based format
	virtual GTwtNode* toTwt(GTwtDoc* pDoc);

	/// Train with all the points in pData
	virtual void train(GData* pData, int labelDims);

	/// See the comment for GIncrementalLearner::trainSparse
	virtual void trainSparse(GSparseMatrix* pData, int labelDims);

	/// See the comment for GSupervisedLearner::predictDistribution
	virtual void predictDistribution(const double* pIn, GPrediction* pOut);

	/// See the comment for GSupervisedLearner::predict
	virtual void predict(const double* pIn, double* pOut);

	/// Discard any training (but not any settings) so it can be trained again
	virtual void clear();

	/// Sets the technique for interpolation. (If you want to use the "Learner" method,
	/// you should call SetInterpolationLearner instead of this method.)
	void setInterpolationMethod(InterpolationMethod eMethod);

	/// Sets the interpolation method to "Learner" and sets the learner to use. If
	/// bTakeOwnership is true, it will delete the learner when this object is deleted.
	void setInterpolationLearner(GSupervisedLearner* pLearner, bool bTakeOwnership);

	/// Adds a copy of pVector to the internal set.
	size_t addVector(const double* pIn, const double* pOut);

	/// Sets the value for elbow room. (This value is only used with incremental training.)
	void setElbowRoom(double d) { m_dElbowRoom = d * d; }

	/// See the comment for GIncrementalLearner::enableIncrementalLearning
	virtual void enableIncrementalLearning(sp_relation& pRelation, int labelDims, double* pMins, double* pRanges);

	/// Adds a vector to the internal set. Also, if the (k+1)th nearest
	/// neighbor of that vector is less than "elbow room" from it, then
	/// the closest neighbor is deleted from the internal set. (You might
	/// be wondering why the decision to delete the closest neighbor is
	/// determined by the distance of the (k+1)th neigbor. This enables a
	/// clump of k points to form in the most frequently sampled locations.
	/// Also, If you make this decision based on a closer neighbor, then big
	/// holes may form in the model if points are sampled in a poor order.)
	/// Call SetElbowRoom to specify the elbow room distance.
	virtual void trainIncremental(const double* pIn, const double* pOut);

	/// Returns the number of neighbors
	int neighborCount() { return m_nNeighbors; }

	/// Returns the internal data set of known instances
	GData* instances() { return m_pInstances; }

	/// Returns the dissimilarity metric
	GRowDistanceScaled* metric() { return m_pDistanceMetric; }

	/// If you set this to true, it will use a hill-climber to optimize the
	/// attribute scaling factors. If you set it to false (the default), it won't.
	void setOptimizeScaleFactors(bool b);

	/// Get the random number generator that was passed to the constructor
	GRand* getRand() { return m_pRand; }

protected:
	/// Finds the nearest neighbors of pVector
	void findNeighbors(const double* pVector);

	/// Interpolate with each neighbor having equal vote
	void interpolateMean(const double* pIn, GPrediction* pOut, double* pOut2);

	/// Interpolate with each neighbor having a linear vote. (Actually it's linear with
	/// respect to the squared distance instead of the distance, because this is faster
	/// to compute.)
	void interpolateLinear(const double* pIn, GPrediction* pOut, double* pOut2);

	/// Interpolates with the provided supervised learning algorithm
	void interpolateLearner(const double* pIn, GPrediction* pOut, double* pOut2);
};


/*
/// This learner considers all training instances with a weight based on the
/// Gaussian of the distance to the instance. It divides the training set
/// into two parts and globally scales each dimension to yield optimal results
/// on the other half. Then it scales the dimensions of each instance locally
/// to optimize further. This isn't a very fast algorithm.
class GKernelInstanceLearner : public GSupervisedLearner
{
protected:
	GRand* m_pRand;
	size_t m_nInstanceCount;
	GAutoScaleInstance* m_pInstances;
	int m_labelDims;
	GData* m_pSet1;
	GData* m_pSet2;
	bool m_bOptimizeGlobally;
	bool m_bOptimizeLocally;
	sp_relation m_pRelation;

public:
	/// nOutputCount specifies the number of output dimensions
	GKernelInstanceLearner(GRand* pRand, int labelDims = 1);

	/// Load from a text-based format
	GKernelInstanceLearner(GTwtNode* pNode, GRand* pRand);

	virtual ~GKernelInstanceLearner();

	/// Save to a text-based format
	virtual GTwtNode* toTwt(GTwtDoc* pDoc);

	/// Train with all the points in pData
	virtual void train(GData* pData, int labelDims);

	/// See the comment for GSupervisedLearner::predictDistribution
	virtual void predictDistribution(const double* pIn, GPrediction* pOut);

	/// Discard any training (but not any settings) so it can be trained again
	virtual void clear();

	/// Specify whether global optimization should be performed
	void SetOptimizeGlobally(bool b) { m_bOptimizeGlobally = b; }
	
	/// Specify whether local optimization should be performed
	void SetOptimizeLocally(bool b) { m_bOptimizeLocally = b; }

	void SetGlobalVars(const double* pVector, GData* pValidationSet);
	void SetClusterVars(const double* pVector, GData* pValidationSet, GAgglomerativeClusterer* pClust, int nCluster);
	void HalfEval(const double* pIn, GPrediction* pOut, GData* pSet);
};
*/

/// An instance-based transduction algorithm
class GNeighborTransducer : public GTransducer
{
protected:
	int m_friendCount;
	GRand* m_pRand;
	int m_intrinsicDims;
	double m_alpha, m_beta;
	bool m_prune;

public:
	GNeighborTransducer(int neighborCount, GRand* pRand);
	
	/// See the comment for GTransducer::transduce
	virtual void transduce(GData* pDataLabeled, GData* pDataUnlabeled, int labelDims);

	/// Specify to use a custom neighbor-finding algorithm
	void useFriends(int intrinsicDims, double alpha, double beta, bool prune);
};


/// This represents a grid of values. It might be useful as a Q-table with Q-learning.
class GInstanceTable : public GIncrementalLearner
{
protected:
	int m_dims;
	size_t* m_pDims;
	size_t* m_pScales;
	double* m_pTable;
	size_t m_product;
	GRand* m_pRand;
	sp_relation m_pRelation;

public:
	/// dims specifies the number of feature dimensions.
	/// pDims specifies the number of discrete zero-based values for each feature dim.
	GInstanceTable(int dims, size_t* pDims, GRand* pRand);
	virtual ~GInstanceTable();

	/// Serialize this table
	virtual GTwtNode* toTwt(GTwtDoc* pDoc);

	/// Returns the number of feature dims
	virtual int featureDims() { return m_dims; }
	
	/// Returns 1
	virtual int labelDims() { return 1; }

	/// Trains this model
	virtual void train(GData* pData, int labelDims);

	/// See the comment for GIncrementalLearner::trainSparse
	virtual void trainSparse(GSparseMatrix* pData, int labelDims);

	/// Predicts a distribution
	virtual void predictDistribution(const double* pIn, GPrediction* pOut);

	/// Makes a prediction
	virtual void predict(const double* pIn, double* pOut);

	/// Clears the internal model
	virtual void clear();

	/// See the comment for GIncrementalLearner::enableIncrementalLearning
	virtual void enableIncrementalLearning(sp_relation& pRelation, int labelDims, double* pMins, double* pRanges);

	/// See the comment for GIncrementalLearner::trainIncremental
	virtual void trainIncremental(const double* pIn, const double* pOut);
};

} // namespace GClasses

#endif // __GKNN_H__
