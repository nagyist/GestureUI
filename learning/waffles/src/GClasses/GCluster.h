/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#ifndef __GCLUSTER_H__
#define __GCLUSTER_H__

#include "GTransform.h"
#include "GLearner.h"
#include <deque>

namespace GClasses {

class GDissimilarityMetric;
class GNeuralNet;

/// The base class for clustering algorithms. Classes that inherit from this
/// class must implement a method named "cluster" which performs clustering, and
/// a method named "whichCluster" which reports which cluster the specified row
/// is determined to be a member of.
class GClusterer : public GTransform
{
protected:
	int m_clusterCount;

public:
	GClusterer(int nClusterCount)
	: GTransform(), m_clusterCount(nClusterCount)
	{
	}

	virtual ~GClusterer()
	{
	}

	/// Clusters pIn and outputs a dataset with one column that specifies
	/// the cluster number for each row.
	virtual GData* doit(GData* pIn)
	{
		cluster(pIn);
		sp_relation pRel = new GUniformRelation(1, m_clusterCount);
		GData* pOut = new GData(pRel);
		size_t nCount = pIn->rows();
		pOut->newRows(nCount);
		for(size_t i = 0; i < nCount; i++)
			pOut->row(i)[0] = (double)whichCluster(i);
		return pOut;
	}

	/// Performs clustering.
	virtual void cluster(GData* pData) = 0;

	/// Reports which cluster the specified row is a member of.
	virtual int whichCluster(size_t nVector) = 0;
};





/// This merges each cluster with its closest neighbor. (The distance between
/// clusters is computed as the distance between the closest members of the
/// clusters times (n^b), where n is the total number of points from both
/// clusters, and b is a balancing factor.
class GAgglomerativeClusterer : public GClusterer
{
protected:
	int m_nDims;
	GData* m_pData;
	int m_nVectorCount;
	int m_nCurrentClusterCount;
	int m_nTargetClusterCount;
	int m_nFirstCluster; // The index of the first cluster-head-point.
	int* m_pNextCluster; // Maps from every cluster-head point to the next cluster-head point. -2 for non-head points. Terminated with a -1.
	int* m_pPrevCluster; // Maps from every cluster-head point to the prev cluster-head point. -2 for non-head points. Terminated with a -1.
	int* m_pNextNeighbor; // Maps from every point to the next point in the same cluster. The list has no end. It cycles back to the head point.
	int* m_pClusters; // Maps from every point to the cluster in which that point belongs. Generated lazily when GetCluster is first called.
	std::deque<int> m_q;
	double m_dBalancingFactor;

public:
	GAgglomerativeClusterer(int nClusterCount);
	virtual ~GAgglomerativeClusterer();

#ifndef NO_TEST_CODE
	/// Performs unit tests for this class. Throws an exception if there is a failure.
	static void test();
#endif // !NO_TEST_CODE

	/// Performs clustering
	virtual void cluster(GData* pData);

	/// Identifies the cluster of the specified row
	virtual int whichCluster(size_t nVector);

	/// (See the class comment about how the distance between clusters is
	/// computed.) So a large value for b will discourage
	/// the merging of large clusters, which will cause the final clusters
	/// to be closer in their final size. If b is close to 0, then it is very
	/// likely that there will be clusters of vastly different sizes in the end.
	void setBalancingFactor(double d) { m_dBalancingFactor = d; }

protected:
	double computeClusterDistance(int a, int b, bool* pAIsBigger);
	void findMerges();
	void merge();
};


/// This is a semi-supervised agglomerative clusterer. It can only handle
/// one output, and it must be nominal. All inputs must be continuous. Also,
/// it assumes that all output values are represented in the training set.
class GAgglomerativeTransducer : public GTransducer
{
protected:
	GData* m_pData;
	int m_featureDims;
	int m_nVectorCount;
	int m_nCurrentClusterCount;
	int m_nFirstCluster; // The index of the first cluster-head-point.
	int* m_pNextCluster; // Maps from every cluster-head point to the next cluster-head point. -2 for non-head points. Terminated with a -1.
	int* m_pPrevCluster; // Maps from every cluster-head point to the prev cluster-head point. -2 for non-head points. Terminated with a -1.
	int* m_pNextNeighbor; // Maps from every point to the next point in the same cluster. The list has no end. It cycles back to the head point.
	std::deque<int> m_q;
	double m_dBalancingFactor;

public:
	GAgglomerativeTransducer();
	virtual ~GAgglomerativeTransducer();

	/// See the comment for GTransducer::transduce.
	/// labelDims must be 1
	virtual void transduce(GData* pDataLabeled, GData* pDataUnlabeled, int labelDims);

	/// (See the class comment about how the distance between clusters is
	/// computed.) So a large value for b will discourage
	/// the merging of large clusters, which will cause the final clusters
	/// to be closer in their final size. If b is close to 0, then it is very
	/// likely that there will be clusters of vastly different sizes in the end.
	void setBalancingFactor(double d) { m_dBalancingFactor = d; }

protected:
#ifndef NO_TEST_CODE
	void makeStateImage(const char* szFilename);
#endif // !NO_TEST_CODE

	double computeClusterDistance(int a, int b, bool* pAIsBigger);
	void findMerges();
	bool merge();
	bool mergeRows(int a1, int b1);
	void mergeLikeClasses(GData* pDataLabeled);
};


/// An implementation of the K-means clustering algorithm.
class GKMeans : public GClusterer
{
protected:
	int m_nDims;
	GData* m_pData;
	size_t m_nClusters;
	int* m_pClusters;
	GRand* m_pRand;

public:
	GKMeans(size_t nClusters, GRand* pRand);
	~GKMeans();

	/// Performs clustering
	virtual void cluster(GData* pData);

	/// Identifies the cluster of the specified row
	virtual int whichCluster(size_t nVector);

protected:
	bool clusterAttempt(int nMaxIterations);
	bool selectSeeds(GData* pSeeds);
};


/// An implementation of the K-medoids clustering algorithm
class GKMedoids : public GClusterer
{
protected:
	size_t* m_pMedoids;
	GDissimilarityMetric* m_pMetric;
	GData* m_pData;
	double m_d;

public:
	GKMedoids(int clusters);
	virtual ~GKMedoids();

	/// Takes ownership of pMetric
	void setMetric(GDissimilarityMetric* pMetric);

	/// Performs clustering
	virtual void cluster(GData* pData);

	/// Identifies the cluster of the specified row
	virtual int whichCluster(size_t nVector);

protected:
	double curErr();
};



class GGraphCutTransducer : public GTransducer
{
protected:
	int m_neighborCount;
	GRand* m_pRand;
	size_t* m_pNeighbors;
	double* m_pDistances;

public:
	GGraphCutTransducer(int neighborCount, GRand* pRand);
	virtual ~GGraphCutTransducer();

	/// See the comment for GTransducer::transduce.
	/// labelDims must be 1
	virtual void transduce(GData* pDataLabeled, GData* pDataUnlabeled, int labelDims);
};


/*
class GNeuralTransducer : public GTransducer
{
protected:
	GRand* m_pRand;
	GNeuralNet* m_pNN;
	std::vector<size_t> m_paramRanges;

public:
	GNeuralTransducer(GRand* pRand);
	virtual ~GNeuralTransducer();
	GNeuralNet* neuralNet() { return m_pNN; }

	void setParams(std::vector<size_t>& ranges);

	/// See the comment for GTransducer::transduce.
	/// labelDims must be 1
	virtual void transduce(GData* pDataLabeled, GData* pDataUnlabeled, int labelDims);
};
*/
} // namespace GClasses

#endif // __GCLUSTER_H__
