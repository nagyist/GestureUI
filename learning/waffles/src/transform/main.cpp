// ----------------------------------------------------------------
// The contents of this file are distributed under the CC0 license.
// See http://creativecommons.org/publicdomain/zero/1.0/
// ----------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include "../GClasses/GApp.h"
#include "../GClasses/GBits.h"
#include "../GClasses/GCluster.h"
#include "../GClasses/GMacros.h"
#include "../GClasses/GData.h"
#include "../GClasses/GImage.h"
#include "../GClasses/GRand.h"
#include "../GClasses/GFile.h"
#include "../GClasses/GTransform.h"
#include "../GClasses/GVec.h"
#include "../GClasses/GHashTable.h"
#include "../GClasses/GHillClimber.h"
#include "../GClasses/GManifold.h"
#include "../GClasses/GNeighborFinder.h"
#include "../GClasses/GNeuralNet.h"
#include "../GClasses/GHeap.h"
#include "../GClasses/GRect.h"
#include "../GClasses/GSparseMatrix.h"
#include "../GClasses/GMath.h"
#include <time.h>
#include <iostream>
#include <string>
#include <set>
#ifdef WIN32
#	include <direct.h>
#	include <process.h>
#endif
#include <exception>
#include "../wizard/usage.h"

using namespace GClasses;
using std::cout;
using std::cerr;
using std::vector;
using std::string;
using std::set;

void parseAttributeList(vector<size_t>& list, GArgReader& args, size_t attrCount)
{
	const char* szList = args.pop_string();
	set<size_t> attrSet;
	while(true)
	{
		// Skip whitespace
		while(*szList <= ' ' && *szList != '\0')
			szList++;

		// Find the next ',' or the end of string, and
		int i;
		int j = -1;
		for(i = 0; szList[i] != '\0' && szList[i] != ','; i++)
		{
			if(j < 0 && szList[i] == '-')
				j = i;
		}
		if(j >= 0)
		{
			while(szList[j + 1] <= ' ' && szList[j + 1] != '\0')
				j++;
		}

		// Add the attributes to the list
		if(i > 0)
		{
			if(*szList < '0' || *szList > '9')
				ThrowError("Expected a number");
			if(j < 0)
			{
#ifdef WIN32
				size_t val = (size_t)_strtoui64(szList, (char**)NULL, 10);
#else
				size_t val = strtoull(szList, (char**)NULL, 10);
#endif
				if(val >= attrCount)
					ThrowError("Invalid attribute index: ", gformat(val), ". (Attributes are zero-indexed.) Valid values are from 0 to ", gformat(attrCount - 1));
				if(attrSet.find(val) != attrSet.end())
					ThrowError("Attribute ", gformat(val), " is listed multiple times");
				attrSet.insert(val);
				list.push_back(val);
			}
			else
			{
				if(szList[j + 1] < '0' || szList[j + 1] > '9')
					ThrowError("Expected a number");
#ifdef WIN32
				size_t beg = (size_t)_strtoui64(szList, (char**)NULL, 10);
				size_t end = (size_t)_strtoui64(szList + j + 1, (char**)NULL, 10);
#else
				size_t beg = strtoull(szList, (char**)NULL, 10);
				size_t end = strtoull(szList + j + 1, (char**)NULL, 10);
#endif
				int step = 1;
				if(end < beg)
					step = -1;
				for(size_t val = beg; true; val += step)
				{
					if(attrSet.find(val) != attrSet.end())
						ThrowError("Attribute ", gformat(val), " is listed multiple times");
					attrSet.insert(val);
						list.push_back(val);
					if(val == end)
						break;
				}
			}
		}

		// Advance
		szList += i;
		if(*szList == '\0')
			break;
		szList++;
	}
}

GNeighborFinder* instantiateNeighborFinder(GData* pData, GRand* pRand, GArgReader& args)
{
	// Get the algorithm name
	GNeighborFinder* pNF = NULL;
	const char* alg = args.pop_string();

	// Parse the options
	int cutCycleLen = 0;
	bool normalize = false;
	while(args.next_is_flag())
	{
		if(args.if_pop("-cyclecut"))
			cutCycleLen = args.pop_uint();
		else if(args.if_pop("-normalize"))
			normalize = true;
		else
			ThrowError("Invalid neighbor finder option: ", args.peek());
	}

	// Parse required algorithms
	if(_stricmp(alg, "bruteforce") == 0)
	{
		int neighbors = args.pop_uint();
		pNF = new GBruteForceNeighborFinder(pData, 0, neighbors, NULL, true);
	}
	else if(_stricmp(alg, "kdtree") == 0)
	{
		int neighbors = args.pop_uint();
		pNF = new GKdTree(pData, 0, neighbors, NULL, true);
	}
	else if(_stricmp(alg, "manifold") == 0)
	{
		int neighbors = args.pop_uint();
		int tangentSpaceDims = args.pop_uint();
		double alpha = args.pop_double();
		double beta = args.pop_double();
		pNF = new GManifoldNeighborFinder(pData, neighbors, neighbors * 4, tangentSpaceDims, alpha, beta, false, pRand);
	}
	else if(_stricmp(alg, "system") == 0)
	{
		GData* pControlData = GData::loadArff(args.pop_string());
		Holder<GData> hControlData(pControlData);
		if(pControlData->rows() != pData->rows())
			ThrowError("mismatching number of rows");
		int neighbors = args.pop_uint();
		pNF = new GDynamicSystemNeighborFinder(pData, hControlData.release(), true, neighbors, pRand);
	}
	else
		ThrowError("Unrecognized neighbor finding algorithm: ", alg);

	// Normalize
	if(normalize)
	{
		GNeighborFinderCacheWrapper* pNF2 = new GNeighborFinderCacheWrapper(pNF, true);
		pNF2->fillCache();
		pNF2->normalizeDistances();
		pNF = pNF2;
	}

	// Apply CycleCut
	if(cutCycleLen > 0)
	{
		GNeighborFinderCacheWrapper* pNF2 = new GNeighborFinderCacheWrapper(pNF, true);
		pNF2->fillCache();
		pNF2->cutShortcuts(cutCycleLen);
		pNF = pNF2;
	}

	return pNF;
}

void agglomerativeclusterer(GArgReader& args)
{
	// Load the file and params
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);
	int clusters = args.pop_uint();

	// Do the clustering
	GAgglomerativeClusterer clusterer(clusters);
	GData* pOut = clusterer.doit(pData);
	Holder<GData> hOut(pOut);
	pOut->print(cout);
}

void AddIndexAttribute(GArgReader& args)
{
	// Parse args
	const char* filename = args.pop_string();
	int nStartValue = 0;
	int nIncrement = 1;
	while(args.size() > 0)
	{
		if(args.if_pop("-start"))
			nStartValue = args.pop_uint();
		else if(args.if_pop("-increment"))
			nIncrement = args.pop_uint();
		else
			ThrowError("Invalid option: ", args.peek());
	}

	GData* pData = GData::loadArff(filename);
	Holder<GData> hData(pData);
	GData indexes(1);
	indexes.newRows(pData->rows());
	for(size_t i = 0; i < pData->rows(); i++)
		indexes.row(i)[0] = nStartValue + i * nIncrement;
	GData* pUnified = GData::mergeHoriz(&indexes, pData);
	Holder<GData> hUnified(pUnified);
	pUnified->print(cout);
}

void addMatrices(GArgReader& args)
{
	GData* pA = GData::loadArff(args.pop_string());
	Holder<GData> hA(pA);
	GData* pB = GData::loadArff(args.pop_string());
	Holder<GData> hB(pB);
	pA->add(pB, false);
	pA->print(cout);
}

void addNoise(GArgReader& args)
{
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);
	double dev = args.pop_double();

	// Parse the options
	unsigned int seed = getpid() * (unsigned int)time(NULL);
	int excludeLast = 0;
	while(args.next_is_flag())
	{
		if(args.if_pop("-seed"))
			seed = args.pop_uint();
		else if(args.if_pop("-excludelast"))
			excludeLast = args.pop_uint();
		else
			ThrowError("Invalid neighbor finder option: ", args.peek());
	}

	GRand prng(seed);
	size_t cols = pData->cols() - excludeLast;
	for(size_t r = 0; r < pData->rows(); r++)
	{
		double* pRow = pData->row(r);
		for(size_t c = 0; c < cols; c++)
			*(pRow++) += dev * prng.normal();
	}
	pData->print(cout);
}

void align(GArgReader& args)
{
	GData* pA = GData::loadArff(args.pop_string());
	Holder<GData> hA(pA);
	GData* pB = GData::loadArff(args.pop_string());
	Holder<GData> hB(pB);
	GData* pC = GData::align(pA, pB);
	Holder<GData> hC(pC);
	pC->print(cout);
}

void attributeSelector(GArgReader& args)
{
	// Load the data
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);

	// Parse the options
	unsigned int seed = getpid() * (unsigned int)time(NULL);
	int labelDims = 1;
	int targetFeatures = 1;
	string outFilename = "";
	while(args.next_is_flag())
	{
		if(args.if_pop("-seed"))
			seed = args.pop_uint();
		else if(args.if_pop("-out"))
		{
			targetFeatures = args.pop_uint();
			outFilename = args.pop_string();
		}
		else if(args.if_pop("-labeldims"))
			labelDims = args.pop_uint();
		else
			ThrowError("Invalid neighbor finder option: ", args.peek());
	}

	// Do the attribute selection
	GRand prng(seed);
	GAttributeSelector as(labelDims, targetFeatures, &prng);
	if(outFilename.length() > 0)
	{
		as.train(pData);
		GData* pDataOut = as.transformBatch(pData);
		Holder<GData> hDataOut(pDataOut);
		cout << "Reduced data saved to " << outFilename.c_str() << ".\n";
		pDataOut->saveArff(outFilename.c_str());
	}
	else
		as.train(pData);
	cout << "\nAttribute rankings from most salient to least salient. (Attributes are zero-indexed.)\n";
	GArffRelation* pRel = (GArffRelation*)pData->relation().get();
	for(size_t i = 0; i < as.ranks().size(); i++)
		cout << as.ranks()[i] << " " << pRel->attrName(as.ranks()[i]) << "\n";
}

void autoCorrelation(GArgReader& args)
{
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);
	size_t lag = MIN((size_t)256, pData->rows() / 2);
	int dims = pData->cols();
	GTEMPBUF(double, mean, dims);
	pData->centroid(mean);
	GData ac(dims + 1);
	for(size_t i = 1; i <= lag; i++)
	{
		double* pRow = ac.newRow();
		*(pRow++) = i;
		for(int j = 0; j < dims; j++)
		{
			*pRow = 0;
			size_t k;
			for(k = 0; k + i < pData->rows(); k++)
			{
				double* pA = pData->row(k);
				double* pB = pData->row(k + i);
				*pRow += (pA[j] - mean[j]) * (pB[j] - mean[j]);
			}
			*pRow /= k;
			pRow++;
		}
	}
	ac.print(cout);
}

void blendEmbeddings(GArgReader& args)
{
	// Load the files and params
	GData* pDataOrig = GData::loadArff(args.pop_string());
	Holder<GData> hDataOrig(pDataOrig);
	unsigned int seed = getpid() * (unsigned int)time(NULL);
	GRand prng(seed);
	GNeighborFinder* pNF = instantiateNeighborFinder(pDataOrig, &prng, args);
	Holder<GNeighborFinder> hNF(pNF);
	GData* pDataA = GData::loadArff(args.pop_string());
	Holder<GData> hDataA(pDataA);
	GData* pDataB = GData::loadArff(args.pop_string());
	Holder<GData> hDataB(pDataB);
	if(pDataA->rows() != pDataOrig->rows() || pDataB->rows() != pDataOrig->rows())
		ThrowError("mismatching number of rows");
	if(pDataA->cols() != pDataB->cols())
		ThrowError("mismatching number of cols");

	// Parse Options
	while(args.size() > 0)
	{
		if(args.if_pop("-seed"))
			prng.setSeed(args.pop_uint());
		else
			ThrowError("Invalid option: ", args.peek());
	}

	// Get a neighbor table
	if(!pNF->isCached())
	{
		GNeighborFinderCacheWrapper* pNF2 = new GNeighborFinderCacheWrapper(hNF.release(), true);
		hNF.reset(pNF2);
		pNF = pNF2;
	}
	((GNeighborFinderCacheWrapper*)pNF)->fillCache();
	size_t* pNeighborTable = ((GNeighborFinderCacheWrapper*)pNF)->cache();

	// Do the blending
	size_t startPoint = (size_t)prng.next(pDataA->rows());
	double* pRatios = new double[pDataA->rows()];
	ArrayHolder<double> hRatios(pRatios);
	GVec::setAll(pRatios, 0.5, pDataA->rows());
	GData* pDataC = GManifold::blendEmbeddings(pDataA, pRatios, pDataB, pNF->neighborCount(), pNeighborTable, startPoint);
	Holder<GData> hDataC(pDataC);
	pDataC->print(cout);
}

void breadthFirstUnfolding(GArgReader& args)
{
	// Load the file and params
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);
	unsigned int nSeed = getpid() * (unsigned int)time(NULL);
	GRand prng(nSeed);
	GNeighborFinder* pNF = instantiateNeighborFinder(pData, &prng, args);
	Holder<GNeighborFinder> hNF(pNF);
	int targetDims = args.pop_uint();

	// Parse Options
	int reps = 1;
	Holder<GData> hControlData(NULL);
	while(args.size() > 0)
	{
		if(args.if_pop("-seed"))
			prng.setSeed(args.pop_uint());
		else if(args.if_pop("-reps"))
			reps = args.pop_uint();
		else
			ThrowError("Invalid option: ", args.peek());
	}

	// Transform the data
	GBreadthFirstUnfolding transform(reps, pNF->neighborCount(), targetDims, &prng);
	transform.setNeighborFinder(pNF);
	GData* pDataAfter = transform.doit(pData);
	Holder<GData> hDataAfter(pDataAfter);
	pDataAfter->print(cout);
}

void center(GArgReader& args)
{
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);
	unsigned int r = args.pop_uint();
	int cols = pData->cols();
	double* pRow = pData->row(r);
	for(size_t i = 0; i < r; r++)
		GVec::subtract(pData->row(i), pRow, cols);
	for(size_t i = r + 1; i < pData->rows(); i++)
		GVec::subtract(pData->row(i), pRow, cols);
	GVec::setAll(pRow, 0.0, cols);
	pData->print(cout);
}

void cholesky(GArgReader& args)
{
	GData* pA = GData::loadArff(args.pop_string());
	Holder<GData> hA(pA);
	pA->cholesky();
	pA->print(cout);
}

void correlation(GArgReader& args)
{
	GData* pA = GData::loadArff(args.pop_string());
	Holder<GData> hA(pA);
	int attr1 = args.pop_uint();
	int attr2 = args.pop_uint();

	// Parse Options
	bool aboutorigin = false;
	while(args.size() > 0)
	{
		if(args.if_pop("-aboutorigin"))
			aboutorigin = true;
		else
			ThrowError("Invalid option: ", args.peek());
	}

	double m1, m2;
	if(aboutorigin)
	{
		m1 = 0;
		m2 = 0;
	}
	else
	{
		m1 = pA->mean(attr1);
		m2 = pA->mean(attr2);
	}
	double corr = pA->linearCorrelationCoefficient(attr1, m1, attr2, m2);
	cout.precision(14);
	cout << corr << "\n";
}

void determinant(GArgReader& args)
{
	GData* pA = GData::loadArff(args.pop_string());
	Holder<GData> hA(pA);
	double d = pA->determinant();
	cout.precision(14);
	cout << d << "\n";
}

void Discretize(GArgReader& args)
{
	// Load the file
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);

	// Parse Options
	int nFirst = 0;
	int nLast = pData->relation()->size() - 1;
	int nBuckets = MAX(2, (int)floor(sqrt((double)pData->rows() + 0.5)));
	while(args.size() > 0)
	{
		if(args.if_pop("-buckets"))
			nBuckets = args.pop_uint();
		else if(args.if_pop("-colrange"))
		{
			nFirst = args.pop_uint();
			nLast = args.pop_uint();
		}
		else
			ThrowError("Invalid option: ", args.peek());
	}
	if(nFirst < 0 || nLast >= pData->relation()->size() || nLast < nFirst)
		ThrowError("column index out of range");

	// Discretize the continuous attributes in the specified range
	for(int i = nFirst; i <= nLast; i++)
	{
		if(pData->relation()->valueCount(i) != 0)
			continue;
		double min, range;
		pData->minAndRange(i, &min, &range);
		for(size_t j = 0; j < pData->rows(); j++)
		{
			double* pPat = pData->row(j);
			pPat[i] = MAX(0, MIN(nBuckets - 1, (int)floor(((pPat[i] - min) * nBuckets) / range)));
		}
		((GArffRelation*)pData->relation().get())->setAttrValueCount(i, nBuckets);
	}

	// Print results
	pData->print(cout);
}

void dropColumns(GArgReader& args)
{
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);
	vector<size_t> colList;
	int attrCount = pData->cols();
	parseAttributeList(colList, args, attrCount);
	std::sort(colList.begin(), colList.end());
	std::reverse(colList.begin(), colList.end());
	for(size_t i = 0; i < colList.size(); i++)
		pData->deleteColumn(colList[i]);
	pData->print(cout);
}

void DropMissingValues(GArgReader& args)
{
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);
	GRelation* pRelation = pData->relation().get();
	int dims = pRelation->size();
	for(size_t i = pData->rows() - 1; i < pData->rows(); i--)
	{
		double* pPat = pData->row(i);
		bool drop = false;
		for(int j = 0; j < dims; j++)
		{
			if(pRelation->valueCount(j) == 0)
			{
				if(pPat[j] == UNKNOWN_REAL_VALUE)
				{
					drop = true;
					break;
				}
			}
			else
			{
				if(pPat[j] == UNKNOWN_DISCRETE_VALUE)
				{
					drop = true;
					break;
				}
			}
		}
		if(drop)
			pData->deleteRow(i);
	}
	pData->print(cout);
}

void dropRows(GArgReader& args)
{
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);
	size_t newSize = args.pop_uint();
	while(pData->rows() > newSize)
		pData->deleteRow(pData->rows() - 1);
	pData->print(cout);
}

void Export(GArgReader& args)
{
	// Load
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);

	// Parse options
	const char* separator = ",";
	while(args.size() > 0)
	{
		if(args.if_pop("-tab"))
			separator = "	";
		else if(args.if_pop("-space"))
			separator = " ";
		else
			ThrowError("Invalid option: ", args.peek());
	}

	// Print
	for(size_t i = 0; i < pData->rows(); i++)
		pData->relation()->printRow(cout, pData->row(i), separator);
}

void Import(GArgReader& args)
{
	// Load the file
	size_t len;
	const char* filename = args.pop_string();
	char* pFile = GFile::loadFile(filename, &len);
	ArrayHolder<char> hFile(pFile);

	// Parse Options
	char separator = ',';
	bool tolerant = false;
	bool columnNamesInFirstRow = false;
	while(args.size() > 0)
	{
		if(args.if_pop("-tab"))
			separator = '\t';
		else if(args.if_pop("-space"))
			separator = ' ';
		else if(args.if_pop("-whitespace"))
			separator = '\0';
		else if(args.if_pop("-semicolon"))
			separator = ';';
		else if(args.if_pop("-separator"))
			separator = args.pop_string()[0];
		else if(args.if_pop("-tolerant"))
			tolerant = true;
		else if(args.if_pop("-columnnames"))
			columnNamesInFirstRow = true;
		else
			ThrowError("Invalid option: ", args.peek());
	}

	// Parse the file
	GData* pData = GData::parseCsv(pFile, len, separator, columnNamesInFirstRow, tolerant);
	Holder<GData> hData(pData);
	((GArffRelation*)pData->relation().get())->setName(filename);

	// Print the data
	pData->print(cout);
}

void isomap(GArgReader& args)
{
	// Load the file and params
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);
	unsigned int nSeed = getpid() * (unsigned int)time(NULL);
	GRand prng(nSeed);
	GNeighborFinder* pNF = instantiateNeighborFinder(pData, &prng, args);
	Holder<GNeighborFinder> hNF(pNF);
	int targetDims = args.pop_uint();

	// Parse Options
	while(args.size() > 0)
	{
		if(args.if_pop("-seed"))
			prng.setSeed(args.pop_uint());
		else
			ThrowError("Invalid option: ", args.peek());
	}

	// Transform the data
	GIsomap transform(pNF->neighborCount(), targetDims, &prng);
	transform.setNeighborFinder(pNF);
	GData* pDataAfter = transform.doit(pData);
	Holder<GData> hDataAfter(pDataAfter);
	pDataAfter->print(cout);
}

void kmeans(GArgReader& args)
{
	// Load the file and params
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);
	int clusters = args.pop_uint();

	// Parse Options
	unsigned int nSeed = getpid() * (unsigned int)time(NULL);
	while(args.size() > 0)
	{
		if(args.if_pop("-seed"))
			nSeed = args.pop_uint();
		else
			ThrowError("Invalid option: ", args.peek());
	}

	// Do the clustering
	GRand prng(nSeed);
	GKMeans clusterer(clusters, &prng);
	GData* pOut = clusterer.doit(pData);
	Holder<GData> hOut(pOut);
	pOut->print(cout);
}

void kmedoids(GArgReader& args)
{
	// Load the file and params
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);
	int clusters = args.pop_uint();

	// Do the clustering
	GKMedoids clusterer(clusters);
	GData* pOut = clusterer.doit(pData);
	Holder<GData> hOut(pOut);
	pOut->print(cout);
}

void lle(GArgReader& args)
{
	// Load the file and params
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);
	unsigned int nSeed = getpid() * (unsigned int)time(NULL);
	GRand prng(nSeed);
	GNeighborFinder* pNF = instantiateNeighborFinder(pData, &prng, args);
	Holder<GNeighborFinder> hNF(pNF);
	int targetDims = args.pop_uint();

	// Parse Options
	while(args.size() > 0)
	{
		if(args.if_pop("-seed"))
			prng.setSeed(args.pop_uint());
		else
			ThrowError("Invalid option: ", args.peek());
	}

	// Transform the data
	GLLE transform(pNF->neighborCount(), targetDims, &prng);
	transform.setNeighborFinder(pNF);
	GData* pDataAfter = transform.doit(pData);
	Holder<GData> hDataAfter(pDataAfter);
	pDataAfter->print(cout);
}

void ManifoldSculpting(GArgReader& args)
{
	// Load the file and params
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);
	unsigned int nSeed = getpid() * (unsigned int)time(NULL);
	GRand prng(nSeed);
	GNeighborFinder* pNF = instantiateNeighborFinder(pData, &prng, args);
	Holder<GNeighborFinder> hNF(pNF);
	int targetDims = args.pop_uint();

	// Parse Options
	const char* szPreprocessedData = NULL;
	double scaleRate = 0.999;
	while(args.size() > 0)
	{
		if(args.if_pop("-seed"))
			prng.setSeed(args.pop_uint());
		else if(args.if_pop("-continue"))
			szPreprocessedData = args.pop_string();
		else if(args.if_pop("-scalerate"))
			scaleRate = args.pop_double();
		else
			ThrowError("Invalid option: ", args.peek());
	}

	// Load the hint data
	GData* pDataHint = NULL;
	Holder<GData> hDataHint(NULL);
	if(szPreprocessedData)
	{
		pDataHint = GData::loadArff(szPreprocessedData);
		hDataHint.reset(pDataHint);
		if(pDataHint->relation()->size() != targetDims)
			ThrowError("Wrong number of dims in the hint data");
		if(pDataHint->rows() != pData->rows())
			ThrowError("Wrong number of patterns in the hint data");
	}

	// Transform the data
	GManifoldSculpting transform(pNF->neighborCount(), targetDims, &prng);
	transform.setSquishingRate(scaleRate);
	if(pDataHint)
		transform.setPreprocessedData(hDataHint.release());
	transform.setNeighborFinder(pNF);
	GData* pDataAfter = transform.doit(pData);
	Holder<GData> hDataAfter(pDataAfter);
	pDataAfter->print(cout);
}
/*
void manifoldSculptingForControl(GArgReader& args)
{
	// Load the file and params
	GData* pDataObs = GData::loadArff(args.pop_string());
	Holder<GData> hDataObs(pDataObs);
	GData* pDataControl = GData::loadArff(args.pop_string());
	Holder<GData> hDataControl(pDataControl);
	int neighbors = args.pop_uint();
	int targetDims = args.pop_uint();

	// Parse Options
	unsigned int nSeed = getpid() * (unsigned int)time(NULL);
	const char* szPreprocessedData = NULL;
	double scaleRate = 0.999;
	double lambda = 0;
	while(args.size() > 0)
	{
		if(args.if_pop("-seed"))
			nSeed = args.pop_uint();
		else if(args.if_pop("-continue"))
			szPreprocessedData = args.pop_string();
		else if(args.if_pop("-scalerate"))
			scaleRate = args.pop_double();
		else if(args.if_pop("-alignconsequences"))
			lambda = args.pop_double();
		else
			ThrowError("Invalid option: ", args.peek());
	}

	// Load the hint data
	GData* pDataHint = NULL;
	Holder<GData> hDataHint(NULL);
	if(szPreprocessedData)
	{
		pDataHint = GData::loadArff(szPreprocessedData);
		hDataHint.reset(pDataHint);
		if(pDataHint->relation()->size() != targetDims)
			ThrowError("Wrong number of dims in the hint data");
		if(pDataHint->rows() != pDataObs->rows())
			ThrowError("Wrong number of patterns in the hint data");
	}

	// Transform the data
	GRand prng(nSeed);
	GManifoldSculptingForControl transform(neighbors, targetDims, &prng, pDataControl, lambda);
	transform.setSquishingRate(scaleRate);
	if(pDataHint)
		transform.setPreprocessedData(hDataHint.release());

	GNeighborFinder* pNF = new GDynamicSystemNeighborFinder(pDataObs, pDataControl, false, neighbors, &prng);
	Holder<GNeighborFinder> hNF(pNF);
	transform.setNeighborFinder(pNF);
	GData* pDataAfter = transform.doit(pDataObs);
	Holder<GData> hDataAfter(pDataAfter);
	pDataAfter->print(cout);
}

void manifoldUnfolder(GArgReader& args)
{
	// Load the file and params
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);
	unsigned int nSeed = getpid() * (unsigned int)time(NULL);
	GRand prng(nSeed);
	GNeighborFinder* pNF = instantiateNeighborFinder(pData, &prng, args);
	Holder<GNeighborFinder> hNF(pNF);
	int targetDims = args.pop_uint();

	// Parse Options
	while(args.size() > 0)
	{
		if(args.if_pop("-seed"))
			prng.setSeed(args.pop_uint());
		else
			ThrowError("Invalid option: ", args.peek());
	}

	// Transform the data
	GManifoldUnfolder transform(pNF->neighborCount(), targetDims, &prng);
	transform.setNeighborFinder(pNF);
	GData* pDataAfter = transform.doit(pData);
	Holder<GData> hDataAfter(pDataAfter);
	pDataAfter->print(cout);
}
*/
void ComputeMeanSquaredError(GData* pData1, GData* pData2, int dims, double* pResults)
{
	GVec::setAll(pResults, 0.0, dims);
	for(size_t i = 0; i < pData1->rows(); i++)
	{
		double* pPat1 = pData1->row(i);
		double* pPat2 = pData2->row(i);
		for(int j = 0; j < dims; j++)
		{
			double d = (*pPat1 - *pPat2);
			pResults[j] += (d * d);
			pPat1++;
			pPat2++;
		}
	}
	GVec::multiply(pResults, 1.0 / pData1->rows(), dims);
}

class FitDataCritic : public GTargetFunction
{
protected:
	GData* m_pData1;
	GData* m_pData2;
	int m_attrs;
	GData m_transformed;
	GData m_transform;
	double* m_pResults;

public:
	FitDataCritic(GData* pData1, GData* pData2, int attrs)
	: GTargetFunction(attrs + attrs * attrs), m_pData1(pData1), m_pData2(pData2), m_attrs(attrs), m_transformed(attrs), m_transform(attrs)
	{
		m_transform.newRows(attrs);
		m_transformed.newRows(pData1->rows());
		m_transform.makeIdentity();
		m_pResults = new double[attrs];
	}

	virtual ~FitDataCritic()
	{
		delete[] m_pResults;
	}

	virtual bool isStable() { return true; }
	virtual bool isConstrained() { return false; }

	virtual void initVector(double* pVector)
	{
		GVec::setAll(pVector, 0.0, m_attrs);
		m_transform.toVector(pVector + m_attrs);
	}

	void TransformData(const double* pVector)
	{
		m_transform.fromVector(pVector + m_attrs, m_attrs);
		for(size_t i = 0; i < m_pData2->rows(); i++)
		{
			double* pPatIn = m_pData2->row(i);
			double* pPatOut = m_transformed.row(i);
			m_transform.multiply(pPatIn, pPatOut);
			GVec::add(pPatOut, pVector, m_attrs);
		}
	}

	virtual double computeError(const double* pVector)
	{
		TransformData(pVector);
		ComputeMeanSquaredError(m_pData1, &m_transformed, m_attrs, m_pResults);
		double sum = GVec::sumElements(m_pResults, m_attrs);
		return sum;
	}

	void ShowResults(const double* pVector)
	{
		TransformData(pVector);
		ComputeMeanSquaredError(m_pData1, &m_transformed, m_attrs, m_pResults);
		//GVec::print(stdout, m_pResults, m_attrs);
		cout.precision(14);
		cout << GVec::sumElements(m_pResults, m_attrs) << "\n";
	}

	const double* GetResults() { return m_pResults; }
};

void MeasureMeanSquaredError(GArgReader& args)
{
	// Load the first file
	GData* pData1 = GData::loadArff(args.pop_string());
	Holder<GData> hData1(pData1);

	// Load the second file
	GData* pData2 = GData::loadArff(args.pop_string());
	Holder<GData> hData2(pData2);

	// check sizes
	if(pData1->relation()->size() != pData2->relation()->size())
		ThrowError("The datasets must have the same number of dims");
	if(pData1->rows() != pData2->rows())
		ThrowError("The datasets must have the same size");

	// Parse Options
	bool fit = false;
	while(args.size() > 0)
	{
		if(args.if_pop("-fit"))
			fit = true;
		else
			ThrowError("Invalid option: ", args.peek());
	}

	int dims = pData1->relation()->size();
	if(fit)
	{
		FitDataCritic critic(pData1, pData2, dims);
		GHillClimber search(&critic);

		double dPrevError;
		double dError = search.iterate();
		cerr.precision(14);
		cerr << dError << "\n";
		cerr.flush();
		while(true)
		{
			dPrevError = dError;
			for(int i = 1; i < 30; i++)
				search.iterate();
			dError = search.iterate();
			cerr << dError << "\n";
			cerr.flush();
			if((dPrevError - dError) / dPrevError < 1e-10)
				break;
		}
		critic.ShowResults(search.currentVector());
	}
	else
	{
		// Compute mean squared error
		GTEMPBUF(double, results, dims);
		ComputeMeanSquaredError(pData1, pData2, dims, results);
		GVec::print(cout, 14, results, dims);
	}
	cout << "\n";
}

void mergeHoriz(GArgReader& args)
{
	GData* pData1 = GData::loadArff(args.pop_string());
	Holder<GData> hData1(pData1);
	GData* pMerged = pData1;
	Holder<GData> hMerged(NULL);
	while(args.size() > 0)
	{
		GData* pData2 = GData::loadArff(args.pop_string());
		Holder<GData> hData2(pData2);
		if(pMerged->rows() != pData2->rows())
			ThrowError("The datasets must have the same number of rows");
		pMerged = GData::mergeHoriz(pMerged, pData2);
		hMerged.reset(pMerged);
	}
	pMerged->print(cout);
}

void mergeVert(GArgReader& args)
{
	GData* pData1 = GData::loadArff(args.pop_string());
	Holder<GData> hData1(pData1);
	GData* pData2 = GData::loadArff(args.pop_string());
	Holder<GData> hData2(pData2);
	if(pData1->cols() != pData2->cols())
		ThrowError("The datasets must have the same number of columns");
	pData1->mergeVert(pData2);
	pData1->print(cout);
}

void multiDimensionalScaling(GArgReader& args)
{
	GRand prng(0);
	GData* pDistances = GData::loadArff(args.pop_string());
	int targetDims = args.pop_uint();

	// Parse Options
	bool useSquaredDistances = false;
	while(args.size() > 0)
	{
		if(args.if_pop("-squareddistances"))
			useSquaredDistances = true;
		else
			ThrowError("Invalid option: ", args.peek());
	}

	GData* pResults = GManifold::multiDimensionalScaling(pDistances, targetDims, &prng, useSquaredDistances);
	Holder<GData> hResults(pResults);
	pResults->print(cout);
}

void multiplyMatrices(GArgReader& args)
{
	GData* pA = GData::loadArff(args.pop_string());
	Holder<GData> hA(pA);
	GData* pB = GData::loadArff(args.pop_string());
	Holder<GData> hB(pB);

	// Parse Options
	bool transposeA = false;
	bool transposeB = false;
	while(args.size() > 0)
	{
		if(args.if_pop("-transposea"))
			transposeA = true;
		else if(args.if_pop("-transposeb"))
			transposeB = true;
		else
			ThrowError("Invalid option: ", args.peek());
	}

	GData* pC = GData::multiply(*pA, *pB, transposeA, transposeB);
	Holder<GData> hC(pC);
	pC->print(cout);
}

void multiplyScalar(GArgReader& args)
{
	GData* pA = GData::loadArff(args.pop_string());
	Holder<GData> hA(pA);
	double scale = args.pop_double();
	pA->multiply(scale);
	pA->print(cout);
}

void normalize(GArgReader& args)
{
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);

	double min = 0.0;
	double max = 1.0;
	while(args.size() > 0)
	{
		if(args.if_pop("-range"))
		{
			min = args.pop_double();
			max = args.pop_double();
		}
		else
			ThrowError("Invalid option: ", args.peek());
	}

	GNormalize transform(min, max);
	transform.train(pData);
	GData* pOut = transform.transformBatch(pData);
	Holder<GData> hOut(pOut);
	pOut->print(cout);
}

void neighbors(GArgReader& args)
{
	// Load the data
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);
	int neighborCount = args.pop_uint();

	// Find the neighbors
	GKdTree neighborFinder(pData, 0, neighborCount, NULL, true);
	GTEMPBUF(size_t, neighbors, neighborCount);
	GTEMPBUF(double, distances, neighborCount);
	double sumClosest = 0;
	double sumAll = 0;
	for(size_t i = 0; i < pData->rows(); i++)
	{
		neighborFinder.neighbors(neighbors, distances, i);
		neighborFinder.sortNeighbors(neighbors, distances);
		sumClosest += sqrt(distances[0]);
		for(int j = 0; j < neighborCount; j++)
			sumAll += sqrt(distances[j]);
	}
	cout.precision(14);
	cout << "average closest neighbor distance = " << (sumClosest / pData->rows()) << "\n";
	cout << "average neighbor distance = " << (sumAll / (pData->rows() * neighborCount)) << "\n";
}

void nominalToCat(GArgReader& args)
{
	// Load the file
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);

	// Parse Options
	int maxValues = 12;
	while(args.size() > 0)
	{
		if(args.if_pop("-maxvalues"))
			maxValues = args.pop_uint();
		else
			ThrowError("Invalid option: ", args.peek());
	}

	// Transform the data
	GNominalToCat transform(maxValues);
	transform.train(pData);
	GData* pDataNew = transform.transformBatch(pData);
	Holder<GData> hDataNew(pDataNew);

	// Print results
	pDataNew->print(cout);
}

void principalComponents(GArgReader& args)
{
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);
	GData results(pData->cols());
	int k = args.pop_uint();
	results.newRows(k);

	// Parse Options
	bool aboutOrigin = false;
	while(args.size() > 0)
	{
		if(args.if_pop("-aboutorigin"))
			aboutOrigin = true;
		else
			ThrowError("Invalid option: ", args.peek());
	}

	// Compute the mean
	GTEMPBUF(double, mean, pData->cols());
	if(aboutOrigin)
		GVec::setAll(mean, 0.0, pData->cols());
	else
		pData->centroid(mean);

	// Compute principal components
	GRand prng(0);
	for(int i = 0; i < k; i++)
	{
		pData->principalComponent(results.row(i), pData->cols(), mean, &prng);
		if(i < k - 1)
			pData->removeComponent(mean, results.row(i), pData->cols());
	}

	results.print(cout);
}

void PrincipleComponentAnalysis(GArgReader& args)
{
	// Load the file
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);

	// Transform the data
	int nTargetDims = args.pop_uint();
	GRand prng(0);
	GPCA transform(nTargetDims, &prng);
	transform.train(pData);
	GData* pDataAfter = transform.transformBatch(pData);
	Holder<GData> hDataAfter(pDataAfter);
	pDataAfter->print(cout);
}

void pseudoInverse(GArgReader& args)
{
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);
	GData* pInverted = pData->pseudoInverse();
	Holder<GData> hInverted(pInverted);
	pInverted->print(cout);
}

void reducedRowEchelonForm(GArgReader& args)
{
	GData* pA = GData::loadArff(args.pop_string());
	Holder<GData> hA(pA);
	pA->toReducedRowEchelonForm();
	pA->print(cout);
}

void ShowBriefUsage(const char* appName)
{
	UsageNode* pUsageTree = makeTransformUsageTree();
	Holder<UsageNode> hUsageTree(pUsageTree);
	pUsageTree->print(0, 3, 76, false);
	cout << "__________________________________\nTo see the full usage info, enter:\n    " << appName << " usage\n";
	cout.flush();
}

void ShowInvalidSyntaxError(const char* appName)
{
	cout << "\nInvalid syntax.\n\n";
	ShowBriefUsage(appName);
}

void ShowUsage(const char* appName)
{
	UsageNode* pUsageTree = makeTransformUsageTree();
	Holder<UsageNode> hUsageTree(pUsageTree);
	pUsageTree->print(0, 3, 76, true);
	cout.flush();
}

void significance(GArgReader& args)
{
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);
	int attr1 = args.pop_uint();
	int attr2 = args.pop_uint();

	// Parse options
	double tolerance = 0.001;
	while(args.size() > 0)
	{
		if(args.if_pop("-tol"))
			tolerance = args.pop_double();
		else
			ThrowError("Invalid option: ", args.peek());
	}

	// Perform the significance tests
	cout.precision(8);
	cout << "\n";
	{
		cout << "Means=" << pData->mean(attr1) << ", " << pData->mean(attr2) << "\n";
		int v;
		double t;
		pData->pairedTTest(&v, &t, attr1, attr2, false);
		double p = GMath::tTestAlphaValue(v, t);
		cout << "Paired T Test: v=" << v << ", t=" << t << ", p=" << p << "\n";
	}
	{
		int v;
		double t;
		pData->pairedTTest(&v, &t, attr1, attr2, true);
		double p = GMath::tTestAlphaValue(v, t);
		cout << "Paired T Test with normalized values: v=" << v << ", t=" << t << ", p=" << p << "\n";
	}
	cout << "\n";
	{
		int less = 0;
		int eq = 0;
		int more = 0;
		for(size_t i = 0; i < pData->rows(); i++)
		{
			double* pRow = pData->row(i);
			if(ABS(pRow[attr1] - pRow[attr2]) < tolerance)
				eq++;
			else if(pRow[attr1] < pRow[attr2])
				less++;
			else
				more++;
		}
		cout << ((double)less / pData->rows()) << "%% less, " << ((double)eq / pData->rows()) << "%% same, " << ((double)more / pData->rows()) << "%% greater\n";
		double t = pData->wilcoxonSignedRanksTest(attr1, attr2, tolerance);
		int n = (int)pData->rows();
		double a = GMath::wilcoxonAlphaValue(n, t);
		cout << "Wilcoxon Signed Ranks Test: n=" << n << ", t=" << t << ", alpha=" << a << "\n";
	}
	cout << "\n";
}

void sparseShuffle(GArgReader& args)
{
	// Load
	GSparseMatrix* pData = GSparseMatrix::load(args.pop_string());
	Holder<GSparseMatrix> hData(pData);

	// Parse options
	unsigned int nSeed = getpid() * (unsigned int)time(NULL);
	while(args.size() > 0)
	{
		if(args.if_pop("-seed"))
			nSeed = args.pop_uint();
		else
			ThrowError("Invalid option: ", args.peek());
	}

	// Shuffle and print
	GRand prng(nSeed);
	pData->shuffle(&prng);
	pData->print(cout);
}

void sparseSplit(GArgReader& args)
{
	// Load
	GSparseMatrix* pData = GSparseMatrix::load(args.pop_string());
	Holder<GSparseMatrix> hData(pData);
	int pats1 = args.pop_uint();
	int pats2 = (int)pData->rows() - pats1;
	if(pats2 < 0)
		ThrowError("out of range. The data only has ", gformat(pData->rows()), " rows.");
	const char* szFilename1 = args.pop_string();
	const char* szFilename2 = args.pop_string();

	// Split
	GSparseMatrix* pPart1 = pData->subMatrix(0, 0, pData->cols(), pats1);
	Holder<GSparseMatrix> hPart1(pPart1);
	GSparseMatrix* pPart2 = pData->subMatrix(0, pats1, pData->cols(), pats2);
	Holder<GSparseMatrix> hPart2(pPart2);
	pPart1->save(szFilename1);
	pPart2->save(szFilename2);
}

void split(GArgReader& args)
{
	// Load
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);
	int pats = (int)pData->rows() - args.pop_uint();
	if(pats < 0)
		ThrowError("out of range. The data only has ", gformat(pData->rows()), " rows.");
	const char* szFilename1 = args.pop_string();
	const char* szFilename2 = args.pop_string();

	// Split
	GData other(pData->relation());
	pData->splitBySize(&other, pats);
	pData->saveArff(szFilename1);
	other.saveArff(szFilename2);
}

void squaredDistance(GArgReader& args)
{
	GData* pA = GData::loadArff(args.pop_string());
	Holder<GData> hA(pA);
	GData* pB = GData::loadArff(args.pop_string());
	Holder<GData> hB(pB);
	double d = pA->sumSquaredDifference(*pB, false);
	cout << "Sum squared distance: " << d << "\n";
	cout << "Mean squared distance: " << (d / pA->rows()) << "\n";
}

void ReplaceMissingValues(GArgReader& args)
{
	// Load
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);

	// Parse options
	unsigned int nSeed = getpid() * (unsigned int)time(NULL);
	while(args.size() > 0)
	{
		if(args.if_pop("-seed"))
			nSeed = args.pop_uint();
		else
			ThrowError("Invalid option: ", args.peek());
	}

	// Replace missing values and print
	GRand prng(nSeed);
	for(int i = 0; i < pData->relation()->size(); i++)
		pData->randomlyReplaceMissingValues(i, &prng);
	pData->shuffle(&prng);
	pData->print(cout);
}

void Shuffle(GArgReader& args)
{
	// Load
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);

	// Parse options
	unsigned int nSeed = getpid() * (unsigned int)time(NULL);
	while(args.size() > 0)
	{
		if(args.if_pop("-seed"))
			nSeed = args.pop_uint();
		else
			ThrowError("Invalid option: ", args.peek());
	}

	// Shuffle and print
	GRand prng(nSeed);
	pData->shuffle(&prng);
	pData->print(cout);
}

void singularValueDecomposition(GArgReader& args)
{
	// Load
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);

	// Parse options
	string ufilename = "u.arff";
	string sigmafilename;
	string vfilename = "v.arff";
	int maxIters = 100;
	while(args.size() > 0)
	{
		if(args.if_pop("-ufilename"))
			ufilename = args.pop_string();
		else if(args.if_pop("-sigmafilename"))
			sigmafilename = args.pop_string();
		else if(args.if_pop("-vfilename"))
			vfilename = args.pop_string();
		else if(args.if_pop("-maxiters"))
			maxIters = args.pop_uint();
		else
			ThrowError("Invalid option: ", args.peek());
	}

	GData* pU;
	double* pDiag;
	GData* pV;
	pData->singularValueDecomposition(&pU, &pDiag, &pV, false, maxIters);
	Holder<GData> hU(pU);
	ArrayHolder<double> hDiag(pDiag);
	Holder<GData> hV(pV);
	pU->saveArff(ufilename.c_str());
	pV->saveArff(vfilename.c_str());
	if(sigmafilename.length() > 0)
	{
		GData sigma(pV->rows());
		sigma.newRows(pU->rows());
		sigma.setAll(0.0);
		size_t m = MIN(sigma.rows(), (size_t)sigma.cols());
		for(size_t i = 0; i < m; i++)
			sigma.row(i)[i] = pDiag[i];
		sigma.saveArff(sigmafilename.c_str());
	}
	else
	{
		GVec::print(cout, 14, pDiag, MIN(pU->rows(), pV->rows()));
		cout << "\n";
	}
}

void SortByAttribute(GArgReader& args)
{
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);
	int nAttr = args.pop_uint();
	int attrCount = pData->relation()->size();
	if(nAttr < 0 || nAttr >= attrCount)
		ThrowError("Index out of range");

	// Parse options
	bool descending = false;
	while(args.size() > 0)
	{
		if(args.if_pop("-descending"))
			descending = true;
		else
			ThrowError("Invalid option: ", args.peek());
	}

	pData->sort(nAttr);
	if(descending)
		pData->reverseRows();
	pData->print(cout);
}

void SwapAttributes(GArgReader& args)
{
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);
	int nAttr1 = args.pop_uint();
	int nAttr2 = args.pop_uint();
	int attrCount = pData->relation()->size();
	if(nAttr1 < 0 || nAttr1 >= attrCount)
		ThrowError("Index out of range");
	if(nAttr2 < 0 || nAttr2 >= attrCount)
		ThrowError("Index out of range");
	pData->swapColumns(nAttr1, nAttr2);
	pData->print(cout);
}

void Transpose(GArgReader& args)
{
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);
	GData* pTransposed = pData->transpose();
	Holder<GData> hTransposed(pTransposed);
	pTransposed->print(cout);
}

void transition(GArgReader& args)
{
	// Load the input data
	GData* pActions = GData::loadArff(args.pop_string());
	Holder<GData> hActions(pActions);
	GData* pState = GData::loadArff(args.pop_string());
	Holder<GData> hState(pState);
	if(pState->rows() != pActions->rows())
		ThrowError("Expected the same number of rows in both datasets");

	// Parse options
	bool delta = false;
	while(args.size() > 0)
	{
		if(args.if_pop("-delta"))
			delta = true;
		else
			ThrowError("Invalid option: ", args.peek());
	}

	// Make the output data
	int actionDims = pActions->cols();
	int stateDims = pState->cols();
	GMixedRelation* pRelation = new GMixedRelation();
	sp_relation pRel = pRelation;
	pRelation->addAttrs(pActions->relation().get());
	pRelation->addAttrs(stateDims + stateDims, 0);
	GData* pTransition = new GData(pRel);
	pTransition->newRows(pActions->rows() - 1);
	for(size_t i = 0; i < pActions->rows() - 1; i++)
	{
		double* pOut = pTransition->row(i);
		GVec::copy(pOut, pActions->row(i), actionDims);
		GVec::copy(pOut + actionDims, pState->row(i), stateDims);
		GVec::copy(pOut + actionDims + stateDims, pState->row(i + 1), stateDims);
		if(delta)
			GVec::subtract(pOut + actionDims + stateDims, pState->row(i), stateDims);
	}
	pTransition->print(cout);
}

void unsupervisedBackProp(GArgReader& args)
{
	// Load the file and params
	GData* pData = GData::loadArff(args.pop_string());
	Holder<GData> hData(pData);
	int targetDims = args.pop_uint();

	// Parse Options
	unsigned int nSeed = getpid() * (unsigned int)time(NULL);
	GRand prng(nSeed);
	GUnsupervisedBackProp ubp(targetDims, &prng);
	vector<size_t> paramRanges;
	while(args.size() > 0)
	{
		if(args.if_pop("-seed"))
			prng.setSeed(args.pop_uint());
		else if(args.if_pop("-addlayer"))
			ubp.neuralNet()->addLayer(args.pop_uint());
		else if(args.if_pop("-learningrate"))
			ubp.neuralNet()->setLearningRate(args.pop_double());
		else if(args.if_pop("-rate"))
			ubp.setRate(args.pop_double());
		else if(args.if_pop("-params"))
		{
			size_t paramDims = args.pop_uint();
			for(size_t i = 0; i < paramDims; i++)
				paramRanges.push_back(args.pop_uint());
		}
		else
			ThrowError("Invalid option: ", args.peek());
	}
	ubp.setParams(paramRanges);

	// Transform the data
	GData* pDataAfter = ubp.doit(pData);
	Holder<GData> hDataAfter(pDataAfter);
	pDataAfter->print(cout);
}

int DoIt(int argc, char *argv[])
{
	PathData pd;
	GFile::parsePath(argv[0], &pd);
	const char* appName = argv[0] + pd.fileStart;

	GArgReader args(argc, argv);
	args.pop_string();
	int ret = 0;
	if(args.size() >= 1)
	{
		if(args.if_pop("usage")) ShowUsage(appName);
		else if(args.if_pop("add")) addMatrices(args);
		else if(args.if_pop("addindexcolumn")) AddIndexAttribute(args);
		else if(args.if_pop("addnoise")) addNoise(args);
		else if(args.if_pop("agglomerative")) agglomerativeclusterer(args);
		else if(args.if_pop("align")) align(args);
		else if(args.if_pop("attributeselector")) attributeSelector(args);
		else if(args.if_pop("autocorrelation")) autoCorrelation(args);
		else if(args.if_pop("blendembeddings")) blendEmbeddings(args);
		else if(args.if_pop("breadthfirstunfolding")) breadthFirstUnfolding(args);
		else if(args.if_pop("nominalToCat")) nominalToCat(args);
		else if(args.if_pop("center")) center(args);
		else if(args.if_pop("cholesky")) cholesky(args);
		else if(args.if_pop("correlation")) correlation(args);
		else if(args.if_pop("determinant")) determinant(args);
		else if(args.if_pop("discretize")) Discretize(args);
		else if(args.if_pop("dropcolumns")) dropColumns(args);
		else if(args.if_pop("dropmissingvalues")) DropMissingValues(args);
		else if(args.if_pop("droprows")) dropRows(args);
		else if(args.if_pop("export")) Export(args);
		else if(args.if_pop("import")) Import(args);
		else if(args.if_pop("isomap")) isomap(args);
		else if(args.if_pop("kmeans")) kmeans(args);
		else if(args.if_pop("kmedoids")) kmedoids(args);
		else if(args.if_pop("lle")) lle(args);
		else if(args.if_pop("manifoldsculpting")) ManifoldSculpting(args);
//		else if(args.if_pop("msfc")) manifoldSculptingForControl(args);
//		else if(args.if_pop("manifoldunfolder")) manifoldUnfolder(args);
		else if(args.if_pop("measuremeansquarederror")) MeasureMeanSquaredError(args);
		else if(args.if_pop("mergehoriz")) mergeHoriz(args);
		else if(args.if_pop("mergevert")) mergeVert(args);
		else if(args.if_pop("multidimensionalscaling")) multiDimensionalScaling(args);
		else if(args.if_pop("multiply")) multiplyMatrices(args);
		else if(args.if_pop("multiplyscalar")) multiplyScalar(args);
		else if(args.if_pop("normalize")) normalize(args);
		else if(args.if_pop("neighbors")) neighbors(args);
		else if(args.if_pop("pca")) PrincipleComponentAnalysis(args);
		else if(args.if_pop("principalcomponents")) principalComponents(args);
		else if(args.if_pop("pseudoinverse")) pseudoInverse(args);
		else if(args.if_pop("replacemissingvalues")) ReplaceMissingValues(args);
		else if(args.if_pop("reducedrowechelonform")) reducedRowEchelonForm(args);
		else if(args.if_pop("unsupervisedbackprop")) unsupervisedBackProp(args);
		else if(args.if_pop("shuffle")) Shuffle(args);
		else if(args.if_pop("significance")) significance(args);
		else if(args.if_pop("sortcolumn")) SortByAttribute(args);
		else if(args.if_pop("sparseshuffle")) sparseShuffle(args);
		else if(args.if_pop("sparsesplit")) sparseSplit(args);
		else if(args.if_pop("split")) split(args);
		else if(args.if_pop("squaredDistance")) squaredDistance(args);
		else if(args.if_pop("svd")) singularValueDecomposition(args);
		else if(args.if_pop("swapcolumns")) SwapAttributes(args);
		else if(args.if_pop("transition")) transition(args);
		else if(args.if_pop("transpose")) Transpose(args);
		else
		{
			ret = 1;
			ShowInvalidSyntaxError(appName);
		}
	}
	else
	{
		ret = 1;
		ShowInvalidSyntaxError(appName);
	}
	return ret;
}

int main(int argc, char *argv[])
{
#ifdef _DEBUG
	GApp::enableFloatingPointExceptions();
#endif
	int nRet = 0;
	try
	{
		nRet = DoIt(argc, argv);
	}
	catch(const std::exception& e)
	{
		cerr << e.what() << "\n";
		nRet = 100;
	}

	return nRet;
}

