/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#ifndef __GGRAPHCUT_H__
#define __GGRAPHCUT_H__

#include "GOptimizer.h"
#include <stdlib.h>
#include <deque>

namespace GClasses {

struct GGraphCutNode;
struct GGraphCutEdge;
class GHeap;
class GRegionAjacencyGraph;
class GGraphEdgeIterator;
class GRand;


/// This implements an optimized max-flow/min-cut algorithm described in
/// "An experimental comparison of min-cut/max-flow algorithms for energy minimization in vision"
/// by Boykov, Y. and Kolmogorov, V.
/// This implementation assumes that edges are undirected
class GGraphCut
{
friend class GGraphEdgeIterator;
protected:
	std::deque<int> m_q;
	GHeap* m_pHeap;
	int m_nNodes;
	struct GGraphCutNode* m_pNodes;
	struct GGraphCutNode* m_pSource;
	struct GGraphCutNode* m_pSink;

public:
	GGraphCut(int nNodes);
	~GGraphCut();

#ifndef NO_TEST_CODE
	/// Performs unit tests for this class. Throws an exception if there is a failure.
	static void test();
#endif // !NO_TEST_CODE

	/// Returns the number of nodes in the graph
	int nodeCount() { return m_nNodes; }

	/// Adds an edge to the graph. You must add all the edges
	/// before calling "Cut". The edge will be stored internally
	/// as a directed edge (from nNode1 to nNode2), but the
	/// Cut method will treat them as undirected edges.
	void addEdge(int nNode1, int nNode2, float fCapacity);

	/// Creates an edge from the node that represents each region
	/// to the node for each of its neighbor regions.
	void getEdgesFromRegionList(GRegionAjacencyGraph* pRegionList);

	/// This computes the cut. nSourceNode is the node that
	/// represents the source, and nSinkNode is the node that
	/// represents the sink.
	void cut(int nSourceNode, int nSinkNode);

	/// Determine whether the specified node is on the source-side
	/// or the sink-side of the cut. (You must call "Cut" before
	/// calling this method.)
	bool isSource(int nNode);

	/// Returns true if the specified node borders the cut
	bool doesBorderTheCut(int nNode);

protected:
	void growNode(int nNode);
	void augmentPath(struct GGraphCutEdge* pEdge);
	void recycleTree(int nChild, int nParent);
	void findAHome(int nNode);
};


/// Iterates over the edges that connect to the
/// specified node
class GGraphEdgeIterator
{
protected:
	GGraphCut* m_pGraph;
	int m_nNode;
	struct GGraphCutEdge* m_pCurrentEdge;

public:
	GGraphEdgeIterator(GGraphCut* pGraph, int nNode);
	~GGraphEdgeIterator();

	/// Starts over with a new node
	void reset(int nNode);

	/// Gets the next edge. Returns false if there isn't
	/// one left to get. If it returns true, pNode, pEdgeWeight, and pOutgoing
	/// will contain the values of the edge. (Note that pOutgoing tells you
	/// which direction the edge is going, even though my implementation of
	/// graph-cut treats them as undirected edges.)
	bool next(int* pNode, float* pEdgeWeight, bool* pOutgoing);
};


/// Computes the shortest-cost path between all pairs of
/// vertices in a graph. Takes O(n^3) time.
class GFloydWarshall
{
protected:
	size_t m_nodes;
	GData* m_pCosts;
	size_t* m_pPaths;

public:
	/// nodes specifies the number of nodes in the graph.
	GFloydWarshall(size_t nodes);
	~GFloydWarshall();

#ifndef NO_TEST_CODE
	/// Performs unit tests for this class. Throws an exception if there is a failure.
	static void test();
#endif

	/// Returns the number of nodes in the graph
	size_t nodeCount() { return m_nodes; }

	/// Adds a directed edge to the graph. (You must call this to add all
	/// the edges before calling compute.)
	void addDirectedEdge(size_t from, size_t to, double cost);

	/// Computes the shortest-cost path between every pair of points
	/// (You must add all the edges before you call compute.)
	void compute();

	/// Returns the smallest cost to get from node "from" to node "to".
	/// (You must call compute before calling this)
	double cost(size_t from, size_t to);

	/// Returns the next node on the shortest path from node "from"
	/// to get to node "goal"
	/// (You must call compute before calling this)
	size_t next(size_t from, size_t goal);

	/// Returns a pointer to the cost matrix
	GData* costMatrix() { return m_pCosts; }

	/// Returns the cost matrix. You are responsible to delete it.
	/// (This class should not be used after this method is called.)
	GData* releaseCostMatrix();

	/// Scans the entire cost matrix, and returns false if any node
	/// is unreachable from any other node. (Assumes compute has already
	/// been called.)
	bool isConnected();
};


/// Finds the shortest path from an origin vertex to all other vertices.
/// Implemented with a binary-heap priority-queue. If the graph is sparse
/// on edges, it will run in about O(n log(n)) time. If the graph is dense,
/// it runs in about O(n^2 log(n))
class GDijkstra
{
protected:
	size_t m_nodes;
	double* m_pCosts;
	size_t* m_pPrevious;
	std::vector<size_t>* m_pNeighbors;
	std::vector<double>* m_pEdgeCosts;

public:
	GDijkstra(size_t nodes);
	~GDijkstra();

#ifndef NO_TEST_CODE
	/// Performs unit tests for this class. Throws an exception if there is a failure.
	static void test();
#endif

	/// Returns the number of nodes in the graph
	size_t nodeCount() { return m_nodes; }

	/// Adds a directed edge to the graph. (You must call this to add all
	/// the edges before calling compute.)
	void addDirectedEdge(size_t from, size_t to, double cost);

	/// Finds the shortest-cost path from the specified origin to every
	/// other point in the graph
	void compute(size_t origin);

	/// Returns the total cost to travel from the origin to the specified target node
	double cost(size_t target);

	/// Returns the previous node on the shortest path from the origin
	/// to the specified vertex
	size_t previous(size_t vertex);
};


/// Computes the number of times that the shortest-path between
/// every pair of points passes over each edge and vertex
class GBrandesBetweennessCentrality
{
protected:
	size_t m_nodeCount;
	std::vector<size_t>* m_pNeighbors;
	double* m_pVertexBetweenness;
	std::vector<double>* m_pEdgeBetweenness;

public:
	GBrandesBetweennessCentrality(size_t nodes);
	~GBrandesBetweennessCentrality();

#ifndef NO_TEST_CODE
	/// Performs unit tests for this class. Throws an exception if there is a failure.
	static void test();
#endif

	/// Returns the number of nodes in the graph
	size_t nodeCount();

	/// Adds a directed edge to the graph. (You must call this to add all
	/// the edges before calling compute.)
	void addDirectedEdge(size_t from, size_t to);

	/// Adds a directed edge if the specified edge does not already
	/// exist. (This method is inefficient if there are a lot of edges.
	/// If you want efficiency, keep track of the edges yourself.)
	void addDirectedEdgeIfNotDupe(size_t from, size_t to);

	/// Computes the betweenness for all nodes. (You must add all the edges
	/// before calling this.)
	void compute();

	/// Returns the betweenness of the specified vertex. (You must call compute
	/// before calling this.) Note that for undirected graphs, you should divide
	/// this value by 2.
	double vertexBetweenness(size_t vertex);

	/// Returns the betweenness of the specified edge. (You must call compute
	/// before calling this.) Note that for undirected graphs, you should divide
	/// this value by 2.
	double edgeBetweenness(size_t vertex, size_t neighborIndex);

	/// Returns the index of the specified neighbor "to" (by iterating over all the
	/// neighbors of "from" until it finds "to"). Returns INVALID_INDEX if not found.
	size_t neighborIndex(size_t from, size_t to);
};


/// This finds all of the atomic cycles (cycles that cannot be divided
/// into two smaller cycles) in a graph
class GAtomicCycleFinder
{
protected:
	size_t m_nodeCount;
	std::vector<size_t>* m_pNeighbors;
	std::vector<size_t>* m_pMirrors;

public:
	GAtomicCycleFinder(size_t nodes);
	virtual ~GAtomicCycleFinder();

#ifndef NO_TEST_CODE
	/// Performs unit tests for this class. Throws an exception if there is a failure.
	static void test();
#endif

	/// Returns the number of nodes in the graph
	size_t nodeCount();

	/// Adds an undirected edge to the graph. (You must call this to add all
	/// the edges before calling compute.)
	void addEdge(size_t a, size_t b);

	/// Adds an undirected edge to the graph if a duplicate edge does not
	/// already exist. (This method is inefficient if there are a lot of edges.
	/// If you want efficiency, keep track of the edges yourself.)
	void addEdgeIfNotDupe(size_t a, size_t b);

	/// Finds all the atomic cycles in the graph, and calls onDetectAtomicCycle
	/// for each one
	void compute();

	/// You must overload this method to receive the cycles as they are detected
	/// (The edges are every torroidally adjacent pair of vertices in cycle.)
	/// If true is returned, it will continue to find more atomic cycles.
	/// If false is returned, it will stop immediately.
	virtual bool onDetectAtomicCycle(std::vector<size_t>& cycle) = 0;
};

} // namespace GClasses

#endif // __GGRAPHCUT_H__
