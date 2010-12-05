// ----------------------------------------------------------------
// The contents of this file are distributed under the CC0 license.
// See http://creativecommons.org/publicdomain/zero/1.0/
// ----------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <exception>
#include <stdio.h>
#include <math.h>
#include <wchar.h>
#include "../GClasses/GMacros.h"
#include "../GClasses/GApp.h"
#include "../GSup/GBezier.h"
#include "../GClasses/GBitTable.h"
#include "../GClasses/GCluster.h"
#include "../GClasses/GData.h"
#include "../GSup/GDate.h"
#include "../GClasses/GDecisionTree.h"
#include "../GSup/GDiff.h"
#include "../GClasses/GEnsemble.h"
#include "../GClasses/GFile.h"
#include "../GClasses/GFourier.h"
#include "../GClasses/GGraph.h"
#include "../GClasses/GHashTable.h"
#include "../GClasses/GHiddenMarkovModel.h"
#include "../GClasses/GHillClimber.h"
#include "../GClasses/GKNN.h"
#include "../GClasses/GLinear.h"
#include "../GClasses/GMacros.h"
#include "../GClasses/GManifold.h"
#include "../GClasses/GMath.h"
#include "../GClasses/GMixtureOfGaussians.h"
#include "../GClasses/GNaiveBayes.h"
#include "../GClasses/GNaiveInstance.h"
#include "../GClasses/GNeighborFinder.h"
#include "../GClasses/GNeuralNet.h"
#include "../GClasses/GPolynomial.h"
#include "../GClasses/GPriorityQueue.h"
#include "../GClasses/GRand.h"
#include "../GSup/GRayTrace.h"
#include "../GClasses/GRegion.h"
#include "../GSup/GSocket.h"
#include "../GClasses/GSparseMatrix.h"
#include "../GClasses/GSpinLock.h"
#include "../GClasses/GStabSearch.h"
#include "../GClasses/GThread.h"
#include "../GClasses/GTime.h"
#include "../GClasses/GTransform.h"
#include "../GClasses/GTwt.h"
#include "../GClasses/GVec.h"

using namespace GClasses;
using std::cerr;

typedef void (*ClassTestFunc)();

struct ClassTest
{
	const char* szName;
	ClassTestFunc pTest;
};

static struct ClassTest testTable[] =
{
	{"GAgglomerativeClusterer", GAgglomerativeClusterer::test},
	{"GAtomicCycleFinder", GAtomicCycleFinder::test},
	{"GAttributeSelector", GAttributeSelector::test},
	{"GBag", GBag::test},
	{"GBaselineLearner", GBaselineLearner::test},
	{"GBezier", GBezier::test},
	{"GBitTable", GBitTable::test},
	{"GBrandesBetweenness", GBrandesBetweennessCentrality::test},
	{"GBucket", GBucket::test},
	{"GCompressor", GCompressor::test},
	{"GCoordVectorIterator", GCoordVectorIterator::test},
	{"GCycleCut", GCycleCut::test},
	{"GData", GData::test},
	{"GDate", TestGDate},
	{"GDecisionTree", GDecisionTree::test},
	{"GDiff", GDiff::test},
	{"GDijkstra", GDijkstra::test},
	{"GFloydWarshall", GFloydWarshall::test},
	{"GFourier", GFourier::test},
	{"GGraphCut", GGraphCut::test},
	{"GHashTable", GHashTable::test},
	{"GHiddenMarkovModel", GHiddenMarkovModel::test},
	{"GKdTree", GKdTree::test},
	{"GKNN", GKNN::test},
	{"GLinearProgramming", GLinearProgramming::test},
	{"GLinearRegressor", GLinearRegressor::test},
	{"GMath", GMath::test},
	{"GManifold", GManifold::test},
	{"GManifoldNeighborFinder", GManifoldNeighborFinder::test},
	{"GMeanMarginsTree", GMeanMarginsTree::test},
	{"GMixtureOfGaussians", GMixtureOfGaussians::test},
	{"GNaiveBayes", GNaiveBayes::test},
	{"GNaiveInstance", GNaiveInstance::test},
	{"GNeuralNet", GNeuralNet::test},
	{"GNeuralNetPseudoInverse", GNeuralNetPseudoInverse::test},
	{"GPCARotateOnly", GPCARotateOnly::test},
	{"GPolynomial", GPolynomial::test},
	{"GPriorityQueue", GPriorityQueue::test},
	{"GProbeSearch", GProbeSearch::test},
	{"GRand", GRand::test},
	{"GShortcutPruner", GShortcutPruner::test},
	{"GSocket", GSocketClient::test},
//	{"GSparseMatrix", GSparseMatrix::test},
	{"GSpinLock", GSpinLock::test},
	{"GSubImageFinder", GSubImageFinder::test},
	{"GSubImageFinder2", GSubImageFinder2::test},
	{"GTwt", GTwtDoc::test},
	{"GVec", GVec::test},
};

bool RunTest(int nTest)
{
	const char* szName = testTable[nTest].szName;
	printf("%s", szName);
	int nSpaces = 70 - strlen(szName);
	for( ; nSpaces > 0; nSpaces--)
		printf(" ");
	fflush(stdout);
	ClassTestFunc pTest = testTable[nTest].pTest;
	bool bPass = false;
	try
	{
		pTest();
		bPass = true;
	}
	catch(...)
	{
	}
	if(bPass)
		printf("Passed\n");
	else
		printf("FAILED!!!\n");
	return bPass;
}

void RunTests()
{
	int testCount = (sizeof(testTable) / sizeof(struct ClassTest));
	for(int i = 0; i < testCount; i++)
		RunTest(i);
	printf("Done.\n");
}

int main(int argc, char *argv[])
{
	GApp::enableFloatingPointExceptions();
	int nRet = 0;
	try
	{
		RunTests();
	}
	catch(const std::exception& e)
	{
		cerr << e.what() << "\n";
		nRet = 1;
	}

	return nRet;
}

