/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#ifndef __GMANIFOLD_H__
#define __GMANIFOLD_H__

#include "GTransform.h"
#include "GMacros.h"
#include "GRand.h"
#include <vector>
#include <deque>
#include "GVec.h"

namespace GClasses {

struct GManifoldSculptingNeighbor;
class GNeighborFinder;
class GNeuralNet;
class GNeuralNetLayer;


/// This class stores static methods that are useful for manifold learning
class GManifold
{
public:
	/// Computes a set of weights for each neighbor to linearly approximate this point
	static void computeNeighborWeights(GData* pData, size_t point, int k, const size_t* pNeighbors, double* pOutWeights);

	/// Aligns and averages two local neighborhoods together. The results will be
	/// centered around the neighborhood mean. The first point will be
	/// the index point, and the rest will be neighborhood points with
	/// an index that is not INVALID_INDEX.
	static GData* blendNeighborhoods(size_t index, GData* pA, double ratio, GData* pB, int neighborCount, size_t* pHood);

	/// Combines two embeddings to form an "average" embedding. pRatios is an array that specifies how much to weight the
	/// neighborhoods around each point. If the ratio for a point is close to zero, pA will be emphasized more. If
	/// the ratio for the point is close to 1, pB will be emphasized more. "seed" specifies a starting point. It will
	/// blend outward in a breadth-first manner.
	static GData* blendEmbeddings(GData* pA, double* pRatios, GData* pB, int neighborCount, size_t* pNeighborTable, size_t seed);

	/// Performs classic MDS. pDistances must be a square matrix, but only the upper-triangle is used.
	/// Each row in the results is one of the result points.
	/// If useSquaredDistances is true, then the values in pDistances are assumed to be squared
	/// distances, rather than normal Euclidean distances.
	static GData* multiDimensionalScaling(GData* pDistances, int targetDims, GRand* pRand, bool useSquaredDistances);

#ifndef NO_TEST_CODE
	static void test();
#endif
};



/// This is the base class of manifold learning (aka non-linear
/// dimensionality reducing) algorithms.
class GManifoldLearner : public GTransform
{
public:
	GManifoldLearner() : GTransform() {}
	GManifoldLearner(GTwtNode* pNode) : GTransform(pNode) {}
	virtual ~GManifoldLearner() {}
};




/// Manifold Sculpting. A non-linear dimensionality reduction algorithm.
/// (See Gashler, Michael S. and Ventura, Dan and Martinez, Tony. Iterative
/// non-linear dimensionality reduction with manifold sculpting. In Advances
/// in Neural Information Processing Systems 20, pages 513–520, MIT Press,
/// Cambridge, MA, 2008.)
class GManifoldSculpting : public GManifoldLearner
{
protected:
	int m_nDimensions;
	int m_nNeighbors;
	int m_nStuffIndex;
	int m_nRecordSize;
	int m_nCurrentDimension;
	int m_nTargetDims;
	int m_nPass;
	double m_scale;
	double m_dAveNeighborDist;
	double m_dSquishingRate;
	double m_dLearningRate;
	double m_minNeighborDist, m_maxNeighborDist;
	std::deque<size_t> m_q;
	sp_relation m_pRelationAfter;
	GRand* m_pRand;
	GData* m_pData;
	unsigned char* m_pMetaData;
	GNeighborFinder* m_pNF;

public:
	GManifoldSculpting(int nNeighbors, int targetDims, GRand* pRand);
	virtual ~GManifoldSculpting();

	/// Perform NLDR.
	virtual GData* doit(GData* pIn);

	virtual sp_relation& relationAfter() { return m_pRelationAfter; }

	GData& data() { return *m_pData; }

	/// Call this before calling SquishPass. pRealSpaceData should be a dataset
	/// of all real values.
	void beginTransform(GData* pRealSpaceData);

	/// Perform one iteration of squishing. Returns a heuristical error value
	double squishPass(size_t nSeedDataPoint);

	/// Set the rate of squishing. (.99 is a good value)
	void setSquishingRate(double d) { m_dSquishingRate = d; }

	/// Returns the current learning rate
	double learningRate() { return m_dLearningRate; }

	/// Returns the average distance between neighbors
	double aveNeighborDist() { return m_dAveNeighborDist; }

	/// Counts the number of times that a point has a neighbor with an
	/// index that is >= nThreshold away from this points index. (If
	/// the manifold is sampled in order such that points are expected
	/// to find neighbors with indexes close to their own, this can serve
	/// to identify when parts of the manifold are too close to other
	/// parts for so many neighbors to be used.)
	int countShortcuts(int nThreshold);

	/// This will takes ownership of pData. (It will modify pData too.)
	/// Specifies reduced dimensionality data created by another
	/// algorithm, which should be refined using Manifold Sculpting
	void setPreprocessedData(GData* pData);

	/// Neighbors that are closer than min or farther than max will be
	/// ignored. The default is [0, 1e20]. It is common to make nNeighbors
	/// big and max small so that the hyper-sphere will define the neighborhood.
	void setMinAndMaxNeighborDist(double min, double max)
	{
		m_minNeighborDist = min;
		m_maxNeighborDist = max;
	}

	/// Partially supervise the specified point.
	void clampPoint(int n);

	/// Specifies to use the neighborhoods determined by the specified neighbor-finder instead of the nearest
	/// Euclidean-distance neighbors. If this method is called, pNF should have the same number of
	/// neighbors and the same dataset as is passed into this class.
	void setNeighborFinder(GNeighborFinder* pNF) { m_pNF = pNF; }

protected:
	inline struct GManifoldSculptingStuff* stuff(int n)
	{
		return (struct GManifoldSculptingStuff*)&m_pMetaData[n * m_nRecordSize + m_nStuffIndex];
	}

	inline struct GManifoldSculptingNeighbor* record(int n)
	{
		return (struct GManifoldSculptingNeighbor*)&m_pMetaData[n * m_nRecordSize];
	}

	/// You can overload this to add some intelligent supervision to the heuristic
	virtual double supervisedError(int nPoint) { return 0; }

	void calculateMetadata(GData* pData);
	double vectorCorrelation(double* pdA, double* pdV, double* pdB);
	double vectorCorrelation2(double squaredScale, int a, int vertex, struct GManifoldSculptingNeighbor* pNeighborRec);
	double computeError(int nPoint);
	int adjustDataPoint(int nPoint, double* pError);
	double averageNeighborDistance(int nDims);
	void plotData(float radius);
	void moveMeanToOrigin();
};

/*
/// An experimental algorithm for using Manifold Sculpting to model dynamic systems.
class GManifoldSculptingForControl : public GManifoldSculpting
{
protected:
	GData* m_pControlData;
	double m_squaredLambda;
	bool m_alignConsequences;

public:
	GManifoldSculptingForControl(int nNeighbors, int nTargetDims, GRand* pRand, GData* pControlData, double lambda)
	: GManifoldSculpting(nNeighbors, nTargetDims, pRand), m_pControlData(pControlData), m_squaredLambda(lambda * lambda)
	{
	}

	virtual ~GManifoldSculptingForControl()
	{
	}

protected:
	virtual double supervisedError(int nPoint);
};
*/

/// Isomap. (A well-known manifold learning algorithm.)
class GIsomap : public GManifoldLearner
{
protected:
	int m_neighborCount;
	int m_targetDims;
	GNeighborFinder* m_pNF;
	GRand* m_pRand;

public:
	GIsomap(int neighborCount, int targetDims, GRand* pRand);
	GIsomap(GTwtNode* pNode);
	virtual ~GIsomap();

	/// Serializes this object
	GTwtNode* toTwt(GTwtDoc* pDoc);

	/// Specifies to use the neighborhoods determined by the specified neighbor-finder instead of the nearest
	/// Euclidean-distance neighbors to establish local linearity. If this method is called, it will also
	/// use the number of neighbors and the data associated with pNF, and ignore the number of neighbors
	/// specified to the constructor, and ignore the data passed to the "transform" method.
	void setNeighborFinder(GNeighborFinder* pNF);

	/// Performs NLDR
	virtual GData* doit(GData* pIn);
};


/// Locally Linear Embedding. (A well-known manifold learning algorithm.)
class GLLE : public GManifoldLearner
{
protected:
	int m_neighborCount;
	int m_targetDims;
	GNeighborFinder* m_pNF;
	GRand* m_pRand;

public:
	GLLE(int neighborCount, int targetDims, GRand* pRand);
	GLLE(GTwtNode* pNode);
	virtual ~GLLE();
	
	/// Serialize this object
	GTwtNode* toTwt(GTwtDoc* pDoc);

	/// Specifies to use the neighborhoods determined by the specified neighbor-finder instead of the nearest
	/// Euclidean-distance neighbors to establish local linearity. If this method is called, it will also
	/// use the number of neighbors and the data associated with pNF, and ignore the number of neighbors
	/// specified to the constructor, and ignore the data passed to the "transform" method.
	void setNeighborFinder(GNeighborFinder* pNF);

	/// Performs NLDR
	virtual GData* doit(GData* pIn);
};

/*
/// An experimental manifold learning algorithm. (Doesn't really work yet.)
class GManifoldUnfolder : public GManifoldLearner
{
protected:
	int m_neighborCount;
	int m_targetDims;
	GNeighborFinder* m_pNF;
	GRand* m_pRand;

public:
	GManifoldUnfolder(int neighborCount, int targetDims, GRand* pRand);
	GManifoldUnfolder(GTwtNode* pNode);
	virtual ~GManifoldUnfolder();
	GTwtNode* toTwt(GTwtDoc* pDoc);
	void setNeighborFinder(GNeighborFinder* pNF);

	/// Performs NLDR
	virtual GData* doit(GData* pIn);

protected:
	GData* unfold(GNeighborFinder* pNF, int targetDims, GRand* pRand);
};
*/


/// A manifold learning algorithm that reduces dimensionality in local
/// neighborhoods, and then stitches the reduced local neighborhoods together
/// using the Kabsch algorithm.
class GBreadthFirstUnfolding : public GManifoldLearner
{
protected:
	int m_reps;
	int m_neighborCount;
	int m_targetDims;
	GNeighborFinder* m_pNF;
	bool m_useMds;
	GRand* m_pRand;

public:
	/// reps specifies the number of times to compute the embedding, and blend the
	/// results together. If you just want fast results, use reps=1.
	GBreadthFirstUnfolding(int reps, int neighborCount, int targetDims, GRand* pRand);
	GBreadthFirstUnfolding(GTwtNode* pNode, GRand* pRand);
	virtual ~GBreadthFirstUnfolding();

	/// Serialize this object
	GTwtNode* toTwt(GTwtDoc* pDoc);

	/// Specify the neighbor finder to use to pick neighbors for this algorithm
	void setNeighborFinder(GNeighborFinder* pNF);

	/// Perform NLDR
	virtual GData* doit(GData* pIn);

	/// Specify to use multi-dimensional scaling instead of PCA to reduce in local patches.
	void useMds(bool b) { m_useMds = b; }

protected:
	void refineNeighborhood(GData* pLocal, size_t rootIndex, size_t* pNeighborTable, double* pDistanceTable);
	GData* reduceNeighborhood(GData* pIn, size_t index, size_t* pNeighborhoods, double* pSquaredDistances);
	GData* unfold(GData* pIn, size_t* pNeighborTable, double* pSquaredDistances, size_t seed, double* pOutWeights);
};


/// A manifold learning algorithm that uses back-propagation to train a neural net model
/// to map from low-dimensional space to high-dimensional space.
class GUnsupervisedBackProp : public GManifoldLearner
{
protected:
	size_t m_intrinsicDims;
	GNeuralNet* m_pNN;
	GRand* m_pRand;
	size_t m_paramDims;
	size_t* m_pParamRanges;
	GCoordVectorIterator m_cvi;
	double m_rate;

public:
	GUnsupervisedBackProp(size_t intrinsicDims, GRand* pRand);
	GUnsupervisedBackProp(GTwtNode* pNode, GRand* pRand);
	virtual ~GUnsupervisedBackProp();

	/// Returns a pointer to the neural network used to model the manifold. Typically, this
	/// is used to add layers to the neural network, or set the learning rate (etc.) before
	/// calling doit.
	GNeuralNet* neuralNet() { return m_pNN; }

	/// An experimental feature that allows one to parameterize the output values.
	void setParams(std::vector<size_t>& paramRanges);

	/// Sets the rate at which the threshold is grown. A typical value is 0.1. Smaller
	/// values will take longer, but may yield better results. Larger values will go faster,
	/// but suffer moew from local optima problems.
	void setRate(double d) { m_rate = d; }
	
	/// Perform NLDR. (This also trains the internal neural network to map from
	/// low-dimensional space to high-dimensional space.)
	virtual GData* doit(GData* pIn);
};


} // namespace GClasses

#endif // __GMANIFOLD_H__
