/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#ifndef __GNEIGHBORFINDER_H__
#define __GNEIGHBORFINDER_H__

#include "GData.h"
#include <vector>
#include <set>
#include <map>
#include <utility>

namespace GClasses {

class GData;
class GRelation;
class GRand;
class GKdNode;
class GBitTable;
class GDissimilarityMetric;
class GSupervisedLearner;

/// This class enables you to define dissimilarity (distance) metrics between two vectors.
/// pScaleFactors is an optional parameter (it can be NULL) that lets the calling class
/// scale the significance of each dimension. Distance metrics that do not mix with
/// this concept may simply ignore any scale factors.
/// Typically, classes that use this should be able to assume that the triangle
/// inequality will hold, but do not necessarily enforce the parallelogram law.
class GDissimilarityMetric
{
protected:
	sp_relation m_pRelation;

public:
	GDissimilarityMetric() {}
	GDissimilarityMetric(GTwtNode* pNode);
	virtual ~GDissimilarityMetric() {}

	/// Serialize this metric to a text-based format
	virtual GTwtNode* toTwt(GTwtDoc* pDoc) = 0;

	/// This must be called before dissimilarity can be called
	virtual void init(sp_relation& pRelation) = 0;

	/// Computes the dissimilarity between the two specified vectors
	virtual double dissimilarity(const double* pA, const double* pB) = 0;

	/// Returns the relation that specifies the meaning of the vector elements
	sp_relation& relation() { return m_pRelation; }

	/// Deserializes a dissimilarity metric
	static GDissimilarityMetric* fromTwt(GTwtNode* pNode);

	/// Returns a pointer to the vector of scale factors
	virtual double* scaleFactors() { return NULL; }

protected:
	GTwtNode* baseTwtNode(GTwtDoc* pDoc, const char* szClassName);
};



/// This uses Euclidean distance for continuous attributes, and
/// Hamming distance for nominal attributes.
class GRowDistance : public GDissimilarityMetric
{
public:
	GRowDistance()
	: GDissimilarityMetric()
	{
	}

	GRowDistance(GTwtNode* pNode);

	virtual ~GRowDistance() {}

	/// See the comment for GDissimilarityMetric::toTwt
	virtual GTwtNode* toTwt(GTwtDoc* pDoc);

	/// See the comment for GDissimilarityMetric::init
	virtual void init(sp_relation& pRelation);

	/// Returns the distance between pA and pB
	virtual double dissimilarity(const double* pA, const double* pB);
};




/// This uses Euclidean distance for continuous attributes, and
/// Hamming distance for nominal attributes.
class GRowDistanceScaled : public GDissimilarityMetric
{
protected:
	double* m_pScaleFactors;

public:
	GRowDistanceScaled() : m_pScaleFactors(NULL) {}
	GRowDistanceScaled(GTwtNode* pNode);

	virtual ~GRowDistanceScaled()
	{
		delete[] m_pScaleFactors;
	}

	/// See the comment for GDissimilarityMetric::toTwt
	virtual GTwtNode* toTwt(GTwtDoc* pDoc);

	/// See the comment for GDissimilarityMetric::init
	virtual void init(sp_relation& pRelation);

	/// Returns the scaled distance between pA and pB
	virtual double dissimilarity(const double* pA, const double* pB);

	/// Returns the vector of scalar values associated with each dimension
	virtual double* scaleFactors() { return m_pScaleFactors; }
};




/// Interpolates between manhattan distance (norm=1), Euclidean distance (norm=2),
/// and Chebyshev distance (norm=infinity). Throws an exception if any of the
/// attributes are nominal.
class GMinkowskiDistance : public GDissimilarityMetric
{
protected:
	double m_norm;

public:
	GMinkowskiDistance(double norm)
	: GDissimilarityMetric(), m_norm(norm)
	{
	}

	GMinkowskiDistance(GTwtNode* pNode);

	/// See the comment for GDissimilarityMetric::toTwt
	virtual GTwtNode* toTwt(GTwtDoc* pDoc);

	/// See the comment for GDissimilarityMetric::init
	virtual void init(sp_relation& pRelation);

	/// Returns the distance (using the norm passed to the constructor) between pA and pB
	virtual double dissimilarity(const double* pA, const double* pB);
};



/// Finds the k-nearest neighbors of each vector in a dataset
class GNeighborFinder
{
protected:
	GData* m_pData;
	int m_neighborCount;

public:
	/// pRelation, pData, and pScaleFactors are all expected to remain
	/// valid for the duration of this object (because it doesn't make
	/// a copy of any of it). If pScaleFactors is NULL, it will use a
	/// value of one for all dimensions.
	/// pData is the dataset that for which all indexes refer. If you
	/// call addVector or releaseVector, it will add or release the vector
	/// from this dataset.
	GNeighborFinder(GData* pData, int neighborCount)
	: m_pData(pData), m_neighborCount(neighborCount)
	{
	}

	virtual ~GNeighborFinder()
	{
	}

	/// Returns the data passed to the constructor of this object
	GData* data() { return m_pData; }

	/// Returns the number of neighbors to find
	int neighborCount() { return m_neighborCount; }

	/// Returns true if this neighbor finder can operate on points that
	/// are not in the dataset passed to the constructor
	virtual bool canGeneralize() { return false; }

	/// Returns true iff the neighbors and distances are pre-computed
	virtual bool isCached() { return false; }

	/// pOutNeighbors should be an array of size neighborCount.
	/// index refers to the point/vector whose neighbors you want to obtain.
	/// The value INVALID_INDEX may be used to fill slots with no point
	/// if necessary.
	virtual void neighbors(size_t* pOutNeighbors, size_t index) = 0;

	/// pOutNeighbors and pOutDistances should both be arrays of size neighborCount.
	/// index refers to the point/vector whose neighbors you want to obtain.
	/// The neighbors are not necessarily sorted, but you can call GNeighborFinder::sortNeighbors
	/// if you want them to be sorted.
	/// If there are not enough points in the data set to fill the
	/// neighbor array, the empty ones will have an index of INVALID_INDEX.
	virtual void neighbors(size_t* pOutNeighbors, double* pOutDistances, size_t index) = 0;

	/// Uses Quick Sort to sort the neighbors from least to most
	/// dissimilar, followed by any slots for with INVALID_INDEX for the index.
	static void sortNeighbors(int neighborCount, size_t* pNeighbors, double* pDistances);

	/// Uses Quick Sort to sort the neighbors from least to most
	/// dissimilar, followed by any slots for with INVALID_INDEX for the index.
	void sortNeighbors(size_t* pNeighbors, double* pDistances);
};




/// This wraps a neighbor finding algorithm. It caches the queries for neighbors
/// for the purpose of improving runtime performance.
class GNeighborFinderCacheWrapper : public GNeighborFinder
{
protected:
	GNeighborFinder* m_pNF;
	bool m_own;
	size_t* m_pCache;
	double* m_pDissims;

public:
	/// If own is true, then this will take ownership of pNF
	GNeighborFinderCacheWrapper(GNeighborFinder* pNF, bool own);
	virtual ~GNeighborFinderCacheWrapper();
	virtual void neighbors(size_t* pOutNeighbors, size_t index);
	virtual void neighbors(size_t* pOutNeighbors, double* pOutDistances, size_t index);

	/// See the comment for GNeighborFinder::isCached
	virtual bool isCached() { return true; }

	/// Returns a pointer to the neighbor finder that this wraps
	GNeighborFinder* wrappedNeighborFinder() { return m_pNF; }

	/// Returns the cache of neighbors. (You should probably call fillCache before calling this.)
	size_t* cache() { return m_pCache; }

	/// Returns the table of squared dissimilarities
	double* squaredDistanceTable() { return m_pDissims; }

	/// Ensures that the cache is populated with data for every index in the dataset
	void fillCache();

	/// Uses CycleCut to remove shortcut connections. (Assumes fillCache has already been called.)
	size_t cutShortcuts(int cycleLen);

	/// Patches any missing neighbors by randomly selecting another of its neighbors to fill both spots.
	void patchMissingSpots(GRand* pRand);

	/// (Re)computes all neighbor distances using the specified metric.
	void fillDistances(GDissimilarityMetric* pMetric);

	/// Normalizes all the neighborhoods so that all neighbor distances are approximately 1.
	void normalizeDistances();

	/// Returns true iff the neighbors form a connected graph. Assumes that fillCache has
	/// already been called;
	bool isConnected();
};



/// Finds the k-nearest neighbors (in a dataset) of an arbitrary vector (which may or may not
/// be in the dataset).
class GNeighborFinderGeneralizing : public GNeighborFinder
{
protected:
	GDissimilarityMetric* m_pMetric;
	bool m_ownMetric;

public:
	GNeighborFinderGeneralizing(GData* pData, int labelDims, int neighborCount, GDissimilarityMetric* pMetric, bool ownMetric);

	virtual ~GNeighborFinderGeneralizing();

	/// Returns true. See the comment for GNeighborFinder::canGeneralize.
	virtual bool canGeneralize() { return true; }

	/// Add a reference (not a copy) of pVector to the data set.
	virtual size_t addVector(double* pVector) = 0;

	/// Releases a vector from the dataset. You are responsible to delete[] the vector
	/// that this returns.
	virtual double* releaseVector(size_t nIndex) = 0;

	/// If you make major changes, you can call this to tell it to rebuild
	/// any optimization structures.
	virtual void reoptimize() = 0;

	/// pOutNeighbors and pOutDistances should both be arrays of size neighborCount.
	/// pInputVector is the vector whose neighbors will be found.
	/// The neighbors are not necessarily sorted, but you can call GNeighborFinder::sortNeighbors
	/// if you want them to be sorted.
	/// If there are not enough points in the data set to fill the
	/// neighbor array, the empty ones will have an index of INVALID_INDEX.
	virtual void neighbors(size_t* pOutNeighbors, double* pOutDistances, const double* pInputVector) = 0;

	using GNeighborFinder::neighbors;
};



/// Finds neighbors by measuring the distance to all points. This one should work properly even if
/// the distance metric does not support the triangle inequality.
class GBruteForceNeighborFinder : public GNeighborFinderGeneralizing
{
public:
	GBruteForceNeighborFinder(GData* pData, int labelDims, int neighborCount, GDissimilarityMetric* pMetric, bool ownMetric);
	virtual ~GBruteForceNeighborFinder();

	/// Add a point-vector
	virtual size_t addVector(double* pVector);

	/// Returns a point-vector (and removes it from the internal set). You are responsible to delete it.
	virtual double* releaseVector(size_t nIndex);

	/// This is a no-op method in this class.
	virtual void reoptimize();

	/// See the comment for GNeighborFinder::neighbors
	virtual void neighbors(size_t* pOutNeighbors, size_t index);

	/// See the comment for GNeighborFinder::neighbors
	virtual void neighbors(size_t* pOutNeighbors, double* pOutDistances, size_t index);

	/// See the comment for GNeighborFinderGeneralizing::neighbors
	virtual void neighbors(size_t* pOutNeighbors, double* pOutDistances, const double* pInputVector);
};




/// An efficient algorithm for finding neighbors.
class GKdTree : public GNeighborFinderGeneralizing
{
protected:
	int m_maxLeafSize;
	int m_size;
	GKdNode* m_pRoot;

public:
	GKdTree(GData* pData, int labelDims, int neighborCount, GDissimilarityMetric* pMetric, bool ownMetric);
	virtual ~GKdTree();

#ifndef NO_TEST_CODE
	/// Performs unit tests for this class. Throws an exception if there is a failure.
	static void test();
#endif

	/// Add a new point-vector to the internal set
	virtual size_t addVector(double* pVector);

	/// Release a point-vector from the internal set
	virtual double* releaseVector(size_t nIndex);

	/// Rebuilds the tree to improve subsequent performance. This should be called after
	/// a significant number of point-vectors are added to or released from the internal set.
	virtual void reoptimize();

	/// See the comment for GNeighborFinder::neighbors
	virtual void neighbors(size_t* pOutNeighbors, size_t index);

	/// See the comment for GNeighborFinder::neighbors
	virtual void neighbors(size_t* pOutNeighbors, double* pOutDistances, size_t index);

	/// See the comment for GNeighborFinderGeneralizing::neighbors
	virtual void neighbors(size_t* pOutNeighbors, double* pOutDistances, const double* pInputVector);

	/// Specify the max number of point-vectors to store in each leaf node.
	void setMaxLeafSize(int n) { m_maxLeafSize = n; }
	
	/// Returns the root node of the kd-tree.
	GKdNode* root() { return m_pRoot; }
	
	/// Build the tree
	GKdNode* buildTree(size_t count, size_t* pIndexes);

	/// Returns true iff the specified point-vector is on the >= side of the specified pivot
	bool isGreaterOrEqual(const double* pPat, int attr, double pivot);

protected:
	/// This is the helper method that finds the neighbors
	void findNeighbors(size_t* pOutNeighbors, double* pOutDistances, const double* pInputVector, size_t nExclude);

	/// Computes a good pivot for the specified attribute, and the goodness of splitting on
	/// that attribute. For continuous attributes, the pivot is the (not scaled) mean and the goodness is
	/// the scaled variance. For nominal attributes, the pivot is the most common value and the
	/// goodness is scaled entropy.
	void computePivotAndGoodness(size_t count, size_t* pIndexes, int attr, double* pOutPivot, double* pOutGoodness);

	/// Moves all the indexes that refer to rows that have a value less than pivot in
	/// the specified attribute to the beginning of the list, and the rest to the end. Returns
	/// the number of rows with a value less than the pivot. For nominal values, not-equal
	/// values are moved to the beginning, and equal values are moved to the end.
	size_t splitIndexes(size_t count, size_t* pIndexes, int attr, double pivot);
};





/// This finds the shortcuts in a table of neighbors and replaces them with INVALID_INDEX.
class GShortcutPruner
{
protected:
	size_t* m_pNeighborhoods;
	size_t m_n;
	int m_k;
	int m_cycleThresh;
	int m_subGraphRange;
	size_t m_cuts;

public:
	/// pNeighborMap is expected to be an array of size n*k, where n is the
	/// number of points, and k is the number of neighbors.
	GShortcutPruner(size_t* pNeighborhoods, size_t n, int k);
	~GShortcutPruner();

#ifndef NO_TEST_CODE
	/// Performs unit tests for this class. Throws an exception if there is a failure.
	static void test();
#endif // NO_TEST_CODE

	/// Sets the cycle-length threshold. (The default is 14.)
	void setCycleThreshold(int cycleThresh) { m_cycleThresh = cycleThresh; }

	/// Sets the sub graph range. (The default is 6.)
	void setSubGraphRange(int range) { m_subGraphRange = range; }

	/// Do the pruning. Returns the number of shortcuts that were removed.
	/// Any atomic cycles in the graph (where neighbors are treated as bi-directional)
	/// with a cycle-length of cycleThresh or bigger indicates the existence of a shortcut
	/// that must be cut. To determine which edge in the cycle is the shortcut, it will
	/// make a subgraph containing all nodes withing "subGraphRange" hops of any vertex
	/// in the cycle, and compute the betweenness centrality of every edge in the sub-graph.
	/// The edge on the cycle with the largest betweenness is determed to be the shortcut,
	/// and is replaced with INVALID_INDEX.
	size_t prune();

	/// Internal method
	void onDetectBigAtomicCycle(std::vector<size_t>& cycle);

protected:
	bool isEveryNodeReachable();
};





/// This finds the shortcuts in a table of neighbors and replaces them with INVALID_INDEX.
class GCycleCut
{
protected:
	size_t* m_pNeighborhoods;
	GData* m_pPoints;
	std::map<std::pair<size_t, size_t>, double> m_capacities;
	std::vector<size_t> m_cuts;
	int m_k;
	int m_cycleThresh;
//	double m_aveDist;
	size_t m_cutCount;

public:
	/// pNeighborMap is expected to be an array of size n*k, where n is the
	/// number pPoints->rows(), and k is the number of neighbors.
	GCycleCut(size_t* pNeighborhoods, GData* pPoints, int k);
	~GCycleCut();

#ifndef NO_TEST_CODE
	/// Performs unit tests for this class. Throws an exception if there is a failure.
	static void test();
#endif // NO_TEST_CODE

	/// Sets the cycle-length threshold. (The default is 14.)
	void setCycleThreshold(int cycleThresh) { m_cycleThresh = cycleThresh; }

	/// Do the cutting. Returns the number of edges that were removed.
	/// Any atomic cycles in the graph (where neighbors are treated as bi-directional)
	/// with a cycle-length of cycleThresh or bigger will be cut. It will make the
	/// smallest cut that removes all big atomic cycles
	size_t cut();

	/// Internal method
	void onDetectBigAtomicCycle(std::vector<size_t>& cycle);

protected:
	bool doAnyBigAtomicCyclesExist();
};





/// This class intelligently selects neighbors for each point in a dataset, such that the neighbors
/// define a good neighborhood for manifold learning. A relaxation technique is used to ensure
/// that neighbors lie on a consistent tangent-space while remaining close to the point. This makes
/// manifold learning possible with difficult (somtimes even self-intersecting) manifolds.
class GManifoldNeighborFinder : public GNeighborFinder
{
protected:
	int m_littleK;
	size_t* m_pNeighborhoods;
	double m_learningRate, m_meanSquaredDist, m_a, m_b;
	int m_windowSize;
	double m_minImprovement;
	bool m_prune;

public:
#ifndef NO_TEST_CODE
	/// Performs unit tests for this class. Throws an exception if there is a failure.
	static void test();
#endif
	/// pData is the set of points for which you want to find a set of neighbors.
	/// littleK is the number of neighbors that you wish to obtain for each point.
	/// bigK is the neighborhood size of candidate neighbors. (4*littleK might be a good value for bigK.)
	/// intrinsicDims is the number of dimensions in the tangent space of the manifold on which the points in pData lie.
	/// alpha specifies how much to weight the penalty of the sine of the angle between tangent-spaces. (2.0 might be a good starting value.)
	/// beta specifies how much to weight the penalty of the projection distance of a point onto the neighbors tangent-space. (1.0 might be a good starting value.)
	/// (The distance between points is inherently given a penalizing weight of 1.0, so you can emphasize distance by
	/// decreasing alpha and beta, and you can de-emphasize distance by increasing alpha and beta.)
	/// prune specifies whether you want it to perform post-processing step that prunes shortcut connections. If prune is true,
	/// then shortcut connections will be replaced with INVALID_INDEX, so some points will have fewer than littleK neighbors.
	GManifoldNeighborFinder(GData* pData, int littleK, int bigK, int intrinsicDims, double alpha, double beta, bool prune, GRand* pRand);
	virtual ~GManifoldNeighborFinder();

	/// See the comment for GNeighborFinder::neighbors
	virtual void neighbors(size_t* pOutNeighbors, size_t index);

	/// See the comment for GNeighborFinder::neighbors
	virtual void neighbors(size_t* pOutNeighbors, double* pOutDistances, size_t index);
};



/// A neighbor finder that specializes in dynamical systems. It determines
/// neighbors by searching for the shortest path of actions between observations,
/// and computes the distance as the length of the path.
class GDynamicSystemNeighborFinder : public GNeighborFinder
{
protected:
	bool m_ownActionsData;
	GData* m_pActions;
	std::vector<GSupervisedLearner*> m_consequenceMaps;
	GRand* m_pRand;

public:
	GDynamicSystemNeighborFinder(GData* pObservations, GData* pActions, bool ownActionsData, int neighborCount, GRand* pRand);
	virtual ~GDynamicSystemNeighborFinder();

	/// Computes the neighbors of the specified vector
	virtual void neighbors(size_t* pOutNeighbors, size_t index);

	/// Computes the neighbors and distances of the specified vector
	virtual void neighbors(size_t* pOutNeighbors, double* pOutDistances, size_t index);

protected:
	/// Returns false if distCap is exceeded, or if the results
	/// are too imprecise to be reliable. Otherwise, returns true, and path is set
	/// to contain the number of times that each action must be performed to travel
	/// from point "from" to point "to".
	bool findPath(size_t from, size_t to, double* path, double distCap);
};

/*
/// An experimental neighbor finder designed to normalize away non-uniform
/// distances in observation space.
class GTemporalNeighborFinder : public GNeighborFinder
{
protected:
	GSupervisedLearner* m_pMap;
	GRand* m_pRand;
	double* m_pBuf;

public:
	GTemporalNeighborFinder(GData* pObservations, int neighborCount, GRand* pRand);
	virtual ~GTemporalNeighborFinder();

	/// Computes the neighbors of the specified vector
	virtual void neighbors(size_t* pOutNeighbors, size_t index);

	/// Computes the neighbors and distances of the specified vector
	virtual void neighbors(size_t* pOutNeighbors, double* pOutDistances, size_t index);
};
*/



/// A simple neighbor-finder that reports the nearest neighbors in the sequence.
/// (That is, the previous and next rows are the closest neighbors.) The distance
/// is sequential distance to the neighbor (not squared).
class GSequenceNeighborFinder : public GNeighborFinder
{
public:
	GSequenceNeighborFinder(GData* pData, int neighborCount);
	virtual ~GSequenceNeighborFinder();
	/// Computes the neighbors of the specified vector
	virtual void neighbors(size_t* pOutNeighbors, size_t index);

	/// Computes the neighbors and distances of the specified vector
	virtual void neighbors(size_t* pOutNeighbors, double* pOutDistances, size_t index);
};


} // namespace GClasses

#endif // __GNEIGHBORFINDER_H__
