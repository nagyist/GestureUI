// ----------------------------------------------------------------
// The contents of this file are distributed under the CC0 license.
// See http://creativecommons.org/publicdomain/zero/1.0/
// ----------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include "../GClasses/GActivation.h"
#include "../GClasses/GApp.h"
#include "../GClasses/GData.h"
#include "../GClasses/GCluster.h"
#include "../GClasses/GDecisionTree.h"
#include "../GClasses/GDistribution.h"
#include "../GClasses/GEnsemble.h"
#include "../GClasses/GFile.h"
#include "../GClasses/GImage.h"
#include "../GClasses/GKernelTrick.h"
#include "../GClasses/GKNN.h"
#include "../GClasses/GMacros.h"
#include "../GClasses/GManifold.h"
#include "../GClasses/GNaiveBayes.h"
#include "../GClasses/GNaiveInstance.h"
#include "../GClasses/GNeuralNet.h"
#include "../GClasses/GRand.h"
#include "../GClasses/GSparseMatrix.h"
#include "../GClasses/GSystemLearner.h"
#include "../GClasses/GTime.h"
#include "../GClasses/GTransform.h"
#include "../GClasses/GTwt.h"
#include "../GClasses/GVec.h"
#include "../wizard/usage.h"
#include <time.h>
#include <iostream>
#ifdef WIN32
#	include <direct.h>
#	include <process.h>
#endif
#include <exception>
#include <string>
#include <vector>
#include <set>

using namespace GClasses;
using std::cout;
using std::cerr;
using std::string;
using std::vector;
using std::set;


GTransducer* InstantiateAlgorithm(GRand* pRand, GArgReader& args);

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

GData* loadData(GArgReader& args, size_t* pLabelDims)
{
  // Load the dataset by extension
  const char* szFilename = args.pop_string();
  PathData pd;
  GFile::parsePath(szFilename, &pd);
  GData* pData = NULL;;
  if(_stricmp(szFilename + pd.extStart, ".arff") == 0)
    pData = GData::loadArff(szFilename);
  else if(_stricmp(szFilename + pd.extStart, ".csv") == 0)
    pData = GData::loadCsv(szFilename, ',', false, false);
  else if(_stricmp(szFilename + pd.extStart, ".dat") == 0)
    pData = GData::loadCsv(szFilename, '\0', false, false);
  else
    ThrowError("Unsupported file format: ", szFilename + pd.extStart);
  Holder<GData> hData(pData);

  // Parse params
  vector<size_t> ignore;
  vector<size_t> labels;
  while(args.next_is_flag())
    {
      if(args.if_pop("-labels"))
	parseAttributeList(labels, args, pData->cols());
      else if(args.if_pop("-ignore"))
	parseAttributeList(ignore, args, pData->cols());
      else
	ThrowError("Invalid agglomerativetransducer option: ", args.peek());
    }

  // Throw out the ignored attributes
  std::sort(ignore.begin(), ignore.end());
  for(size_t i = ignore.size() - 1; i < ignore.size(); i--)
    {
      pData->deleteColumn(ignore[i]);
      for(size_t j = 0; j < labels.size(); j++)
	{
	  if(labels[j] >= ignore[i])
	    {
	      if(labels[j] == ignore[i])
		ThrowError("Attribute ", gformat(labels[j]), " is both ignored and used as a label");
	      labels[j]--;
	    }
	}
    }

  // Swap label columns to the end
  *pLabelDims = MAX((size_t)1, labels.size());
  for(size_t i = 0; i < labels.size(); i++)
    {
      size_t src = labels[i];
      size_t dst = pData->cols() - *pLabelDims + i;
      if(src != dst)
	{
	  pData->swapColumns(src, dst);
	  for(size_t j = i + 1; j < labels.size(); j++)
	    {
	      if(labels[j] == dst)
		{
		  labels[j] = src;
		  break;
		}
	    }
	}
    }

  return hData.release();
}

GAgglomerativeTransducer* InstantiateAgglomerativeTransducer(GRand* pRand, GArgReader& args)
{
  GAgglomerativeTransducer* pTransducer = new GAgglomerativeTransducer();
  while(args.next_is_flag())
    {
      if(args.if_pop("-balancingfactor"))
	pTransducer->setBalancingFactor(args.pop_double());
      else
	ThrowError("Invalid agglomerativetransducer option: ", args.peek());
    }
  return pTransducer;
}

GBag* InstantiateBag(GRand* pRand, GArgReader& args)
{
  GBag* pEnsemble = new GBag(pRand);
  while(args.size() > 0)
    {
      if(args.if_pop("end"))
	break;
      int instance_count = args.pop_uint();
      int arg_pos = args.get_pos();
      for(int i = 0; i < instance_count; i++)
	{
	  args.set_pos(arg_pos);
	  GTransducer* pLearner = InstantiateAlgorithm(pRand, args);
	  if(!pLearner->canGeneralize())
	    {
	      delete(pLearner);
	      ThrowError("bag does not support algorithms that cannot generalize.");
	    }
	  pEnsemble->addLearner((GSupervisedLearner*)pLearner);
	}
    }
  return pEnsemble;
}

GBucket* InstantiateBucket(GRand* pRand, GArgReader& args)
{
  GBucket* pEnsemble = new GBucket(pRand);
  while(args.size() > 0)
    {
      if(args.if_pop("end"))
	break;
      GTransducer* pLearner = InstantiateAlgorithm(pRand, args);
      if(!pLearner->canGeneralize())
	{
	  delete(pLearner);
	  ThrowError("crossvalidationselector does not support algorithms that cannot generalize.");
	}
      pEnsemble->addLearner((GSupervisedLearner*)pLearner);
    }
  return pEnsemble;
}

GBaselineLearner* InstantiateBaseline(GRand* pRand, GArgReader& args)
{
  GBaselineLearner* pModel = new GBaselineLearner(pRand);
  return pModel;
}

GDecisionTree* InstantiateDecisionTree(GRand* pRand, GArgReader& args)
{
  GDecisionTree* pModel = new GDecisionTree(pRand);
  while(args.next_is_flag())
    {
      if(args.if_pop("-random"))
	pModel->useRandomDivisions(args.pop_uint());
      else if(args.if_pop("-leafthresh"))
	pModel->setLeafThresh(args.pop_uint());
      else
	ThrowError("Invalid decisiontree option: ", args.peek());
    }
  return pModel;
}

GFilter* InstantiateDiscretize(GRand* pRand, GArgReader& args)
{
  bool features = true;
  bool labels = true;
  int buckets = -1;
  while(args.next_is_flag())
    {
      if(args.if_pop("-notfeatures"))
	features = false;
      else if(args.if_pop("-notlabels"))
	labels = false;
      else if(args.if_pop("-buckets"))
	buckets = args.pop_uint();
      else
	ThrowError("Invalid discretize option: ", args.peek());
    }
  if(args.size() < 1)
    ThrowError("expected an algorithm after \"discretize\".");
  GTransducer* pLearner = InstantiateAlgorithm(pRand, args);
  if(!pLearner->canGeneralize())
    {
      delete(pLearner);
      ThrowError("discretize does not support algorithms that cannot generalize.");
    }
  GFilter* pTL = new GFilter((GSupervisedLearner*)pLearner, true);
  if(features)
    pTL->setFeatureTransform(new GDiscretize(buckets), true);
  if(labels)
    pTL->setLabelTransform(new GDiscretize(buckets), true);
  return pTL;
}

GGraphCutTransducer* InstantiateGraphCutTransducer(GRand* pRand, GArgReader& args)
{
  if(args.size() < 1)
    ThrowError("The number of neighbors must be specified for graphcuttransducer");
  int neighborCount = args.pop_uint();
  GGraphCutTransducer* pTransducer = new GGraphCutTransducer(neighborCount, pRand);
  return pTransducer;
}

GKernelMachine* InstantiateKernelMachine(GRand* pRand, GArgReader& args)
{
  GKernelMachine* pModel = new GKernelMachine();
  return pModel;
}

GKNN* InstantiateKNN(GRand* pRand, GArgReader& args)
{
  if(args.size() < 1)
    ThrowError("The number of neighbors must be specified for knn");
  int neighborCount = args.pop_uint();
  GKNN* pModel = new GKNN(neighborCount, pRand);
  while(args.next_is_flag())
    {
      if(args.if_pop("-equalweight"))
	pModel->setInterpolationMethod(GKNN::Mean);
      else if(args.if_pop("-scalefeatures"))
	pModel->setOptimizeScaleFactors(true);
      else
	ThrowError("Invalid knn option: ", args.peek());
    }
  return pModel;
}

GMeanMarginsTree* InstantiateMeanMarginsTree(GRand* pRand, GArgReader& args)
{
  GMeanMarginsTree* pModel = new GMeanMarginsTree(pRand);
  return pModel;
}
/*
  GModerateNet* InstantiateModerateNet(GRand* pRand, GArgReader& args)
  {
  GModerateNet* pModel = new GModerateNet(pRand);
  while(args.next_is_flag())
  {
  if(args.if_pop("-addlayer"))
  pModel->addLayer(args.pop_uint());
  else if(args.if_pop("-learningrate"))
  pModel->setLearningRate(args.pop_double());
  else if(args.if_pop("-momentum"))
  pModel->setMomentum(args.pop_double());
  else if(args.if_pop("-windowepochs"))
  pModel->setIterationsPerValidationCheck(args.pop_uint());
  else if(args.if_pop("-minwindowimprovement"))
  pModel->setMinImprovement(args.pop_double());
  else if(args.if_pop("-lambda"))
  pModel->setLambda(args.pop_double());
  else if(args.if_pop("-squash"))
  {
  const char* szSF = args.pop_string();
  GActivationFunction* pSF = NULL;
  if(strcmp(szSF, "logistic") == 0)
  pSF = new GActivationLogistic();
  else if(strcmp(szSF, "arctan") == 0)
  pSF = new GActivationArcTan();
  else if(strcmp(szSF, "tanh") == 0)
  pSF = new GActivationTanH();
  else if(strcmp(szSF, "algebraic") == 0)
  pSF = new GActivationAlgebraic();
  else if(strcmp(szSF, "identity") == 0)
  pSF = new GActivationIdentity();
  else if(strcmp(szSF, "bend") == 0)
  pSF = new GActivationBend();
  else if(strcmp(szSF, "bidir") == 0)
  pSF = new GActivationBiDir();
  else
  ThrowError("Unrecognized activation function: ", szSF);
  pModel->setActivationFunction(pSF, true);
  }
  else if(args.if_pop("-crossentropy"))
  pModel->setBackPropTargetFunction(GNeuralNet::cross_entropy);
  else if(args.if_pop("-rootcubed"))
  pModel->setBackPropTargetFunction(GNeuralNet::root_cubed_error);
  else
  ThrowError("Invalid neuralnet option: ", args.peek());
  }
  return pModel;
  }
*/
GNaiveBayes* InstantiateNaiveBayes(GRand* pRand, GArgReader& args)
{
  GNaiveBayes* pModel = new GNaiveBayes(pRand);
  while(args.next_is_flag())
    {
      if(args.if_pop("-ess"))
	pModel->setEquivalentSampleSize(args.pop_double());
      else
	ThrowError("Invalid naivebayes option: ", args.peek());
    }
  return pModel;
}

GNaiveInstance* InstantiateNaiveInstance(GRand* pRand, GArgReader& args)
{
  if(args.size() < 1)
    ThrowError("The number of neighbors must be specified for naiveinstance");
  int neighborCount = args.pop_uint();
  GNaiveInstance* pModel = new GNaiveInstance(neighborCount);
  return pModel;
}

GNeighborTransducer* InstantiateNeighborTransducer(GRand* pRand, GArgReader& args)
{
  if(args.size() < 1)
    ThrowError("The number of neighbors must be specified for neighbortransducer");
  int friendCount = args.pop_uint();
  GNeighborTransducer* pTransducer = new GNeighborTransducer(friendCount, pRand);
  bool prune = false;
  double alpha, beta;
  int intrinsicDims;
  while(args.next_is_flag())
    {
      if(args.if_pop("-prune"))
	prune = true;
      else if(args.if_pop("-friends"))
	{
	  intrinsicDims = args.pop_uint();
	  alpha = args.pop_double();
	  beta = args.pop_double();
	}
      else
	ThrowError("Invalid neighbortransducer option: ", args.peek());
    }
  return pTransducer;
}

GNeuralNet* InstantiateNeuralNet(GRand* pRand, GArgReader& args)
{
  GNeuralNet* pModel = new GNeuralNet(pRand);
  while(args.next_is_flag())
    {
      if(args.if_pop("-addlayer"))
	pModel->addLayer(args.pop_uint());
      else if(args.if_pop("-learningrate"))
	pModel->setLearningRate(args.pop_double());
      else if(args.if_pop("-momentum"))
	pModel->setMomentum(args.pop_double());
      else if(args.if_pop("-windowepochs"))
	pModel->setIterationsPerValidationCheck(args.pop_uint());
      else if(args.if_pop("-minwindowimprovement"))
	pModel->setMinImprovement(args.pop_double());
      else if(args.if_pop("-squash"))
	{
	  const char* szSF = args.pop_string();
	  GActivationFunction* pSF = NULL;
	  if(strcmp(szSF, "logistic") == 0)
	    pSF = new GActivationLogistic();
	  else if(strcmp(szSF, "arctan") == 0)
	    pSF = new GActivationArcTan();
	  else if(strcmp(szSF, "tanh") == 0)
	    pSF = new GActivationTanH();
	  else if(strcmp(szSF, "algebraic") == 0)
	    pSF = new GActivationAlgebraic();
	  else if(strcmp(szSF, "identity") == 0)
	    pSF = new GActivationIdentity();
	  else if(strcmp(szSF, "bend") == 0)
	    pSF = new GActivationBend();
	  else if(strcmp(szSF, "bidir") == 0)
	    pSF = new GActivationBiDir();
	  else if(strcmp(szSF, "piecewise") == 0)
	    pSF = new GActivationPiecewise();
	  else
	    ThrowError("Unrecognized activation function: ", szSF);
	  pModel->setActivationFunction(pSF, true);
	}
      else if(args.if_pop("-crossentropy"))
	pModel->setBackPropTargetFunction(GNeuralNet::cross_entropy);
      else if(args.if_pop("-rootcubed"))
	pModel->setBackPropTargetFunction(GNeuralNet::root_cubed_error);
      else
	ThrowError("Invalid neuralnet option: ", args.peek());
    }
  return pModel;
}
/*
  GNeuralTransducer* InstantiateNeuralTransducer(GRand* pRand, GArgReader& args)
  {
  GNeuralTransducer* pTransducer = new GNeuralTransducer(pRand);
  vector<size_t> paramDims;
  while(args.next_is_flag())
  {
  if(args.if_pop("-addlayer"))
  pTransducer->neuralNet()->addLayer(args.pop_uint());
  else if(args.if_pop("-params"))
  {
  size_t count = args.pop_uint();
  for(size_t i = 0; i < count; i++)
  paramDims.push_back(args.pop_uint());
  }
  else
  ThrowError("Invalid agglomerativetransducer option: ", args.peek());
  }
  pTransducer->setParams(paramDims);
  return pTransducer;
  }
*/
GFilter* InstantiateNormalize(GRand* pRand, GArgReader& args)
{
  bool features = true;
  bool labels = true;
  double min = 0.0;
  double max = 1.0;
  while(args.next_is_flag())
    {
      if(args.if_pop("-notfeatures"))
	features = false;
      else if(args.if_pop("-notlabels"))
	labels = false;
      else if(args.if_pop("-range"))
	{
	  min = args.pop_double();
	  max = args.pop_double();
	}
      else
	ThrowError("Invalid normalize option: ", args.peek());
    }
  if(args.size() < 1)
    ThrowError("expected an algorithm after \"normalize\".");
  GTransducer* pLearner = InstantiateAlgorithm(pRand, args);
  if(!pLearner->canGeneralize())
    {
      delete(pLearner);
      ThrowError("normalize does not support algorithms that cannot generalize.");
    }
  GFilter* pTL = new GFilter((GSupervisedLearner*)pLearner, true);
  if(features)
    pTL->setFeatureTransform(new GNormalize(min, max), true);
  if(labels)
    pTL->setLabelTransform(new GNormalize(min, max), true);
  return pTL;
}

GFilter* InstantiateNominalToCat(GRand* pRand, GArgReader& args)
{
  bool features = true;
  bool labels = true;
  int maxValues = 12;
  while(args.next_is_flag())
    {
      if(args.if_pop("-notfeatures"))
	features = false;
      else if(args.if_pop("-notlabels"))
	labels = false;
      else if(args.if_pop("-maxvalues"))
	maxValues = args.pop_uint();
      else
	ThrowError("Invalid nominaltocat option: ", args.peek());
    }
  if(args.size() < 1)
    ThrowError("expected an algorithm after \"nominaltocat\".");
  GTransducer* pLearner = InstantiateAlgorithm(pRand, args);
  if(!pLearner->canGeneralize())
    {
      delete(pLearner);
      ThrowError("nominaltocat does not support algorithms that cannot generalize.");
    }
  GFilter* pTL = new GFilter((GSupervisedLearner*)pLearner, true);
  if(features)
    pTL->setFeatureTransform(new GNominalToCat(maxValues), true);
  if(labels)
    pTL->setLabelTransform(new GNominalToCat(maxValues), true);
  return pTL;
}

GFilter* InstantiatePCA(GRand* pRand, GArgReader& args)
{
  if(args.size() < 1)
    ThrowError("The number of target dimensions must be specified for pca");
  int targetDims = args.pop_uint();
  if(args.size() < 1)
    ThrowError("expected an algorithm after \"pca ", gformat(targetDims), "\".");
  GTransducer* pLearner = InstantiateAlgorithm(pRand, args);
  if(!pLearner->canGeneralize())
    {
      delete(pLearner);
      ThrowError("pca does not support algorithms that cannot generalize.");
    }
  GFilter* pTL = new GFilter((GSupervisedLearner*)pLearner, true);
  pTL->setFeatureTransform(new GPCA(targetDims, pRand), true);
  return pTL;
}

GTransducer* InstantiateAlgorithm(GRand* pRand, GArgReader& args)
{
  if(args.size() < 1)
    ThrowError("No algorithm specified.");
  else if(args.if_pop("agglomerativetransducer"))
    return InstantiateAgglomerativeTransducer(pRand, args);
  else if(args.if_pop("bag"))
    return InstantiateBag(pRand, args);
  else if(args.if_pop("baseline"))
    return InstantiateBaseline(pRand, args);
  else if(args.if_pop("bucket"))
    return InstantiateBucket(pRand, args);
  else if(args.if_pop("nominaltocat"))
    return InstantiateNominalToCat(pRand, args);
  else if(args.if_pop("decisiontree"))
    return InstantiateDecisionTree(pRand, args);
  else if(args.if_pop("discretize"))
    return InstantiateDiscretize(pRand, args);
  else if(args.if_pop("graphcuttransducer"))
    return InstantiateGraphCutTransducer(pRand, args);
  else if(args.if_pop("kernelmachine"))
    return InstantiateKernelMachine(pRand, args);
  else if(args.if_pop("knn"))
    return InstantiateKNN(pRand, args);
  else if(args.if_pop("meanmarginstree"))
    return InstantiateMeanMarginsTree(pRand, args);
  //	else if(args.if_pop("moderatenet"))
  //		return InstantiateModerateNet(pRand, args);
  else if(args.if_pop("naivebayes"))
    return InstantiateNaiveBayes(pRand, args);
  else if(args.if_pop("naiveinstance"))
    return InstantiateNaiveInstance(pRand, args);
  else if(args.if_pop("neighbortransducer"))
    return InstantiateNeighborTransducer(pRand, args);
  else if(args.if_pop("neuralnet"))
    return InstantiateNeuralNet(pRand, args);
  //	else if(args.if_pop("neuraltransducer"))
  //		return InstantiateNeuralTransducer(pRand, args);
  else if(args.if_pop("normalize"))
    return InstantiateNormalize(pRand, args);
  else if(args.if_pop("pca"))
    return InstantiatePCA(pRand, args);
  ThrowError("Unrecognized algorithm name: ", args.peek());
  return NULL;
}

void ShowUsage(const char* appName)
{
  UsageNode* pUsageTree = makeLearnUsageTree();
  Holder<UsageNode> hUsageTree(pUsageTree);
  pUsageTree->print(0, 3, 76, true);
  UsageNode* pUsageTree2 = makeAlgorithmUsageTree();
  Holder<UsageNode> hUsageTree2(pUsageTree2);
  pUsageTree2->print(0, 3, 76, true);
  cout.flush();
}

void ShowBriefUsage(const char* appName)
{
  UsageNode* pUsageTree = makeLearnUsageTree();
  Holder<UsageNode> hUsageTree(pUsageTree);
  pUsageTree->print(0, 3, 76, false);
  UsageNode* pUsageTree2 = makeAlgorithmUsageTree();
  Holder<UsageNode> hUsageTree2(pUsageTree2);
  pUsageTree2->print(0, 3, 76, false);
  cout << "__________________________________\nTo see the full usage info, enter:\n    " << appName << " usage\n";
  cout.flush();
}

void Train(GArgReader& args)
{
  // Parse options
  unsigned int seed = getpid() * (unsigned int)time(NULL);
  while(args.next_is_flag())
    {
      if(args.if_pop("-seed"))
	seed = args.pop_uint();
      else
	ThrowError("Invalid train option: ", args.peek());
    }

  // Load the data
  GRand prng(seed);
  if(args.size() < 1)
    ThrowError("No dataset specified.");
  size_t labelDims;
  GData* pData = loadData(args, &labelDims);
  Holder<GData> hData(pData);

  // Instantiate the modeler
  GTransducer* pSupLearner = InstantiateAlgorithm(&prng, args);
  Holder<GTransducer> hModel(pSupLearner);
  if(args.size() > 0)
    ThrowError("Superfluous argument: ", args.peek());
  if(!pSupLearner->canGeneralize())
    ThrowError("This algorithm cannot be \"trained\". It can only be used to \"transduce\".");
  GSupervisedLearner* pModel = (GSupervisedLearner*)pSupLearner;

  // Train the modeler
  pModel->train(pData, labelDims);

  // Output the trained model
  GTwtDoc doc;
  GTwtNode* pRoot = pModel->toTwt(&doc);
  doc.setRoot(pRoot);
  doc.write(cout);
}

void predict(GArgReader& args)
{
  // Parse options
  unsigned int seed = getpid() * (unsigned int)time(NULL);
  while(args.next_is_flag())
    {
      if(args.if_pop("-seed"))
	seed = args.pop_uint();
      else
	ThrowError("Invalid predict option: ", args.peek());
    }

  // Load the model
  GRand prng(seed);
  GTwtDoc doc;
  if(args.size() < 1)
    ThrowError("Model not specified.");
  doc.load(args.pop_string());
  GLearnerLoader ll(true);
  GSupervisedLearner* pModeler = ll.loadModeler(doc.root(), &prng);
  Holder<GSupervisedLearner> hModeler(pModeler);

  // Load the data
  if(args.size() < 1)
    ThrowError("No dataset specified.");
  size_t labelDims;
  GData* pData = loadData(args, &labelDims);
  Holder<GData> hData(pData);
  if(labelDims != (size_t)pModeler->labelDims())
    ThrowError("The model was trained with ", gformat(pModeler->labelDims()), " label dims, but the specified dataset has ", gformat(labelDims));

  // Test
  for(size_t i = 0; i < pData->rows(); i++)
    {
      double* pPat = pData->row(i);
      pModeler->predict(pPat, pPat + pData->cols() - pModeler->labelDims());
    }

  // Print results
  pData->print(cout);
}

void predictOnePattern(GArgReader& args)
{
  // Parse options
  unsigned int seed = getpid() * (unsigned int)time(NULL);
  while(args.next_is_flag())
    {
      if(args.if_pop("-seed"))
	seed = args.pop_uint();
      else
	ThrowError("Invalid predictonepattern option: ", args.peek());
    }

  // Load the model
  GRand prng(seed);
  GTwtDoc doc;
  if(args.size() < 1)
    ThrowError("Model not specified.");
  doc.load(args.pop_string());
  GLearnerLoader ll(true);
  GSupervisedLearner* pModeler = ll.loadModeler(doc.root(), &prng);
  Holder<GSupervisedLearner> hModeler(pModeler);

  // Load the dataset
  size_t labelDims;
  GData* pData = loadData(args, &labelDims);
  Holder<GData> hData(pData);
  if(labelDims != (size_t)pModeler->labelDims())
    ThrowError("The model was trained with ", gformat(pModeler->labelDims()), " label dims, but the specified dataset has ", gformat(labelDims));
  if(pData->relation()->type() != GRelation::ARFF)
    ThrowError("Expected a dataset with ARFF metadata");
  GArffRelation* pRel = (GArffRelation*)pData->relation().get();

  // Parse the pattern
  int featureDims = pModeler->featureDims();
  GTEMPBUF(double, pattern, featureDims);
  for(int i = 0; i < featureDims; i++)
    pattern[i] = pRel->parseValue(i, args.pop_string());

  // Predict
  GPrediction* out = new GPrediction[pModeler->labelDims()];
  ArrayHolder<GPrediction> hOut(out);
  pModeler->predictDistribution(pattern, out);

  // Display the prediction
  cout.precision(8);
  for(int i = 0; i < pModeler->labelDims(); i++)
    {
      if(i > 0)
	cout << ", ";
      if(pRel->valueCount(featureDims + i) == 0)
	cout << out[i].mode();
      else
	{
	  string s;
	  pRel->attrValue(&s, featureDims + i, (int)out[i].mode());
	  cout << s.c_str();
	}
    }
  cout << "\n\n";

  // Display the confidence values
  for(int i = 0; i < pModeler->labelDims(); i++)
    {
      if(out[i].isContinuous())
	{
	  GNormalDistribution* pNorm = out[i].asNormal();
	  cout << pRel->attrName(featureDims + i) << ") Normal: predicted mean=" << pNorm->mean() << " predicted variance=" << pNorm->variance() << "\n";
	}
      else
	{
	  GCategoricalDistribution* pCat = out[i].asCategorical();
	  cout << pRel->attrName(featureDims + i) << ") Categorical Confidences: {";
	  double* pValues = pCat->values(pCat->valueCount());
	  for(int j = 0; j < pCat->valueCount(); j++)
	    {
	      if(j > 0)
		cout << ", ";
	      string s;
	      pRel->attrValue(&s, featureDims + i, j);
	      cout << s << "=" << pValues[j];
	    }
	  cout << "}\n";
	}
    }
}

void Test(GArgReader& args)
{
  // Parse options
  unsigned int seed = getpid() * (unsigned int)time(NULL);
  while(args.next_is_flag())
    {
      if(args.if_pop("-seed"))
	seed = args.pop_uint();
      else
	ThrowError("Invalid test option: ", args.peek());
    }

  // Load the model
  GRand prng(seed);
  GTwtDoc doc;
  if(args.size() < 1)
    ThrowError("Model not specified.");
  doc.load(args.pop_string());
  GLearnerLoader ll(true);
  GSupervisedLearner* pModeler = ll.loadModeler(doc.root(), &prng);
  Holder<GSupervisedLearner> hModeler(pModeler);

  // Load the data
  if(args.size() < 1)
    ThrowError("No dataset specified.");
  size_t labelDims;
  GData* pData = loadData(args, &labelDims);
  Holder<GData> hData(pData);
  if(labelDims != (size_t)pModeler->labelDims())
    ThrowError("The model was trained with ", gformat(pModeler->labelDims()), " label dims, but the specified dataset has ", gformat(labelDims));

  // Test
  GTEMPBUF(double, results, pModeler->labelDims());
  pModeler->accuracy(pData, results);
  GVec::print(cout, 14, results, pModeler->labelDims());
  cout << "\n";
}

void Transduce(GArgReader& args)
{
  // Parse options
  unsigned int seed = getpid() * (unsigned int)time(NULL);
  while(args.next_is_flag())
    {
      if(args.if_pop("-seed"))
	seed = args.pop_uint();
      else
	ThrowError("Invalid transduce option: ", args.peek());
    }

  // Load the labeled set
  GRand prng(seed);
  if(args.size() < 1)
    ThrowError("No labeled set specified.");

  size_t labelDims1, labelDims2;
  GData* pDataLabeled = loadData(args, &labelDims1);
  Holder<GData> hDataLabeled(pDataLabeled);
  GData* pDataUnlabeled = loadData(args, &labelDims2);
  Holder<GData> hDataUnlabeled(pDataUnlabeled);
  if(pDataLabeled->cols() != pDataUnlabeled->cols())
    ThrowError("The labeled and unlabeled datasets must have the same number of columns. (The labels in the unlabeled set are just place-holders, and will be overwritten.)");
  if(labelDims1 != labelDims2)
    ThrowError("The labeled and unlabeled datasets must have the same number of label dims. (The labels in the unlabeled set are just place-holders, and will be overwritten.)");

  // Instantiate the modeler
  GTransducer* pSupLearner = InstantiateAlgorithm(&prng, args);
  Holder<GTransducer> hModel(pSupLearner);
  if(args.size() > 0)
    ThrowError("Superfluous argument: ", args.peek());

  // Transduce
  pSupLearner->transduce(pDataLabeled, pDataUnlabeled, labelDims1);

  // Print results
  pDataUnlabeled->print(cout);
}

void TransductiveAccuracy(GArgReader& args)
{
  // Parse options
  unsigned int seed = getpid() * (unsigned int)time(NULL);
  while(args.next_is_flag())
    {
      if(args.if_pop("-seed"))
	seed = args.pop_uint();
      else
	ThrowError("Invalid transacc option: ", args.peek());
    }

  // Load the labeled set
  GRand prng(seed);
  size_t labelDims1, labelDims2;
  GData* pDataTrain = loadData(args, &labelDims1);
  Holder<GData> hDataTrain(pDataTrain);
  GData* pDataTest = loadData(args, &labelDims2);
  Holder<GData> hDataTest(pDataTest);
  if(pDataTrain->cols() != pDataTest->cols())
    ThrowError("The train and test datasets must have the same number of columns. (The labels in the test set are just place-holders, and will be overwritten.)");
  if(labelDims1 != labelDims2)
    ThrowError("The train and test datasets must have the same number of label dims. (The labels in the test set are just place-holders, and will be overwritten.)");

  // Instantiate the modeler
  GTransducer* pSupLearner = InstantiateAlgorithm(&prng, args);
  Holder<GTransducer> hModel(pSupLearner);
  if(args.size() > 0)
    ThrowError("Superfluous argument: ", args.peek());

  // Transduce and measure accuracy
  GTEMPBUF(double, results, labelDims1);
  pSupLearner->trainAndTest(pDataTrain, pDataTest, labelDims1, results);

  // Print results
  GVec::print(cout, 14, results, labelDims1);
  cout << "\n";
}

void SplitTest(GArgReader& args)
{
  // Parse options
  unsigned int seed = getpid() * (unsigned int)time(NULL);
  double trainRatio = 0.5;
  int reps = 1;
  while(args.next_is_flag())
    {
      if(args.if_pop("-seed"))
	seed = args.pop_uint();
      else if(args.if_pop("-trainratio"))
	trainRatio = args.pop_double();
      else if(args.if_pop("-reps"))
	reps = args.pop_uint();
      else
	ThrowError("Invalid splittest option: ", args.peek());
    }
  if(trainRatio < 0 || trainRatio > 1)
    ThrowError("trainratio must be between 0 and 1");

  // Load the data
  GRand prng(seed);
  if(args.size() < 1)
    ThrowError("No dataset specified.");
  size_t labelDims;
  GData* pData = loadData(args, &labelDims);
  Holder<GData> hData(pData);

  // Instantiate the modeler
  GTransducer* pSupLearner = InstantiateAlgorithm(&prng, args);
  Holder<GTransducer> hModel(pSupLearner);
  if(args.size() > 0)
    ThrowError("Superfluous argument: ", args.peek());

  // Do the reps
  int trainingPatterns = MAX((size_t)1, MIN(pData->rows() - 1, (size_t)floor(pData->rows() * trainRatio + 0.5)));
  int testPatterns = pData->rows() - trainingPatterns;
  GTEMPBUF(double, results, 2 * labelDims);
  double* repResults = results + labelDims;
  GVec::setAll(results, 0, labelDims);
  for(int i = 0; i < reps; i++)
    {
      // Shuffle and split the data
      pData->shuffle(&prng);
      GData dataTest(pData->relation(), pData->heap());
      dataTest.reserve(testPatterns);
      pData->splitBySize(&dataTest, testPatterns);

      // Test and print results
      pSupLearner->trainAndTest(pData, &dataTest, labelDims, repResults);
      cout << "rep " << i << ") ";
      GVec::print(cout, 14, repResults, labelDims);
      cout << "\n";
      double weight = 1.0 / (i + 1);
      GVec::multiply(results, 1.0 - weight, labelDims);
      GVec::addScaled(results, weight, repResults, labelDims);
      pData->mergeVert(&dataTest);
    }
  cout << "-----\n";
  GVec::print(cout, 14, results, labelDims);
  cout << "\n";
}

void CrossValidateCallback(void* pSupLearner, int nRep, int nFold, int labelDims, double* pFoldResults)
{
  cout << "Rep: " << nRep << ", Fold: " << nFold <<", Accuracy: ";
  GVec::print(cout, 14, pFoldResults, labelDims);
  cout << "\n";
}

void CrossValidate(GArgReader& args)
{
  // Parse options
  unsigned int seed = getpid() * (unsigned int)time(NULL);
  int reps = 5;
  int folds = 2;
  bool succinct = false;
  while(args.next_is_flag())
    {
      if(args.if_pop("-seed"))
	seed = args.pop_uint();
      else if(args.if_pop("-reps"))
	reps = args.pop_uint();
      else if(args.if_pop("-folds"))
	folds = args.pop_uint();
      else if(args.if_pop("-succinct"))
	succinct = true;
      else
	ThrowError("Invalid crossvalidate option: ", args.peek());
    }
  if(reps < 1)
    ThrowError("There must be at least 1 rep.");
  if(folds < 2)
    ThrowError("There must be at least 2 folds.");

  // Load the data
  if(args.size() < 1)
    ThrowError("No dataset specified.");
  size_t labelDims;
  GData* pData = loadData(args, &labelDims);
  Holder<GData> hData(pData);

  // Instantiate the modeler
  GRand prng(seed);
  GTransducer* pSupLearner = InstantiateAlgorithm(&prng, args);
  Holder<GTransducer> hModel(pSupLearner);
  if(args.size() > 0)
    ThrowError("Superfluous argument: ", args.peek());

  // Test
  cout.precision(8);
  GData* pResults = pSupLearner->repValidate(pData, reps, folds, labelDims, &prng, succinct ? NULL : CrossValidateCallback, pSupLearner);
  Holder<GData> hResults(pResults);
  if(!succinct)
    cout << "-----\n";
  for(size_t i = 0; i < labelDims; i++)
    {
      int attr = pData->cols() - labelDims + i;
      double mean = pResults->mean(i);
      double variance = pResults->variance(i, mean);
      if(!succinct)
	cout << "Attr: " << attr << ", Mean accuracy: ";
      cout << mean;
      if(succinct)
	{
	  if(i + 1 < labelDims)
	    cout << ", ";
	}
      else
	cout << ", Deviation: " << sqrt(variance) << "\n";
    }
  cout << "\n";
}

void trainSparse(GArgReader& args)
{
  // Parse options
  unsigned int seed = getpid() * (unsigned int)time(NULL);
  int labelDims = 1;
  while(args.next_is_flag())
    {
      if(args.if_pop("-seed"))
	seed = args.pop_uint();
      else if(args.if_pop("-labeldims"))
	labelDims = args.pop_uint();
      else
	ThrowError("Invalid trainsparse option: ", args.peek());
    }

  // Load the data
  if(args.size() < 1)
    ThrowError("No dataset specified.");
  GSparseMatrix* pData = GSparseMatrix::load(args.pop_string());
  Holder<GSparseMatrix> hData(pData);
  if(labelDims < 1 || (labelDims > (int)pData->cols()))
    ThrowError("labelDims out of range");

  // Instantiate the modeler
  GRand prng(seed);
  GTransducer* pSupLearner = InstantiateAlgorithm(&prng, args);
  Holder<GTransducer> hModel(pSupLearner);
  if(args.size() > 0)
    ThrowError("Superfluous argument: ", args.peek());
  if(!pSupLearner->canTrainIncrementally())
    ThrowError("This algorithm cannot be trained with a sparse matrix. Only incremental learners (such as naivebayes or neuralnet) support this functionality.");
  GIncrementalLearner* pModel = (GIncrementalLearner*)pSupLearner;

  // Train the modeler
  pModel->trainSparse(pData, labelDims);

  // Output the trained model
  GTwtDoc doc;
  GTwtNode* pRoot = pModel->toTwt(&doc);
  doc.setRoot(pRoot);
  doc.write(cout);
}

void predictSparse(GArgReader& args)
{
  // Parse options
  unsigned int seed = getpid() * (unsigned int)time(NULL);
  while(args.next_is_flag())
    {
      if(args.if_pop("-seed"))
	seed = args.pop_uint();
      else
	ThrowError("Invalid predictsparse option: ", args.peek());
    }

  // Load the model
  GRand prng(seed);
  GTwtDoc doc;
  if(args.size() < 1)
    ThrowError("Model not specified.");
  doc.load(args.pop_string());
  GLearnerLoader ll(true);
  GSupervisedLearner* pModeler = ll.loadModeler(doc.root(), &prng);
  Holder<GSupervisedLearner> hModeler(pModeler);

  // Load the data
  if(args.size() < 1)
    ThrowError("No dataset specified.");
  GSparseMatrix* pData = GSparseMatrix::load(args.pop_string());
  Holder<GSparseMatrix> hData(pData);

  // Test
  double* pFullRow = new double[pData->cols()];
  double* pLabels = pFullRow + pData->cols() - pModeler->labelDims();
  ArrayHolder<double> hFullRow(pFullRow);
  for(unsigned int i = 0; i < pData->rows(); i++)
    {
      pData->fullRow(pFullRow, i);
      pModeler->predict(pFullRow, pLabels);
      cout << i << ") ";
      GVec::print(cout, 7, pLabels, pModeler->labelDims());
      cout << "\n";
    }
}

void vette(string& s)
{
  for(size_t i = 0; i < s.length(); i++)
    {
      if(s[i] <= ' ' || s[i] == '\'' || s[i] == '"')
	s[i] = '_';
    }
}

void PrecisionRecall(GArgReader& args)
{
  // Parse options
  unsigned int seed = getpid() * (unsigned int)time(NULL);
  int reps = 5;
  int samples = 100;
  while(args.next_is_flag())
    {
      if(args.if_pop("-seed"))
	seed = args.pop_uint();
      else if(args.if_pop("-reps"))
	reps = args.pop_uint();
      else if(args.if_pop("-samples"))
	samples = args.pop_uint();
      else
	ThrowError("Invalid precisionrecall option: ", args.peek());
    }
  if(reps < 1)
    ThrowError("There must be at least 1 rep.");
  if(samples < 2)
    ThrowError("There must be at least 2 samples.");

  // Load the data
  if(args.size() < 1)
    ThrowError("No dataset specified.");
  size_t labelDims;
  GData* pData = loadData(args, &labelDims);
  Holder<GData> hData(pData);

  // Instantiate the modeler
  GRand prng(seed);
  GTransducer* pSupLearner = InstantiateAlgorithm(&prng, args);
  Holder<GTransducer> hModel(pSupLearner);
  if(args.size() > 0)
    ThrowError("Superfluous argument: ", args.peek());
  if(!pSupLearner->canGeneralize())
    ThrowError("This algorithm cannot be \"trained\". It can only be used to \"transduce\".");
  GSupervisedLearner* pModel = (GSupervisedLearner*)pSupLearner;

  // Build the relation for the results
  sp_relation pRelation;
  pRelation = new GArffRelation();
  ((GArffRelation*)pRelation.get())->setName("untitled");
  GArffRelation* pRel = (GArffRelation*)pRelation.get();
  pRel->addAttribute("recall", 0, NULL);
  int featureDims = pData->cols() - labelDims;
  for(size_t i = 0; i < labelDims; i++)
    {
      int valCount = MAX(1, pData->relation()->valueCount(featureDims + i));
      for(int val = 0; val < valCount; val++)
	{
	  string s = "precision_";
	  if(pData->relation()->type() == GRelation::ARFF)
	    s += ((GArffRelation*)pData->relation().get())->attrName(i);
	  else
	    {
	      s += "attr";
	      s += i;
	    }
	  if(valCount > 1)
	    {
	      s += "_";
	      pData->relation()->attrValue(&s, featureDims + i, val);
	    }
	  vette(s);
	  pRel->addAttribute(s.c_str(), 0, NULL);
	}
    }

  // Measure precision/recall
  GData results(pRelation);
  results.newRows(samples);
  for(int i = 0; i < samples; i++)
    results.row(i)[0] = (double)i / samples;
  int pos = 1;
  for(size_t i = 0; i < labelDims; i++)
    {
      int valCount = MAX(1, pData->relation()->valueCount(featureDims + i));
      double* precision = new double[valCount * samples];
      ArrayHolder<double> hPrecision(precision);
      pModel->precisionRecall(precision, samples, pData, labelDims, i, reps, &prng);
      for(int j = 0; j < valCount; j++)
	results.setCol(pos++, precision + samples * j);
    }
  GAssert(pos == pRelation->size()); // counting problem
  results.print(cout);
}

class MyRecurrentModel : public GRecurrentModel
{
protected:
  const char* m_stateFilename;
  double m_validateInterval;
  double m_dStart;

public:
  MyRecurrentModel(GSupervisedLearner* pTransition, GSupervisedLearner* pObservation, int actionDims, int contextDims, int obsDims, GRand* pRand, std::vector<size_t>* pParamDims, const char* stateFilename, double validateInterval)
    : GRecurrentModel(pTransition, pObservation, actionDims, contextDims, obsDims, pRand, pParamDims), m_stateFilename(stateFilename), m_validateInterval(validateInterval)
  {
    m_dStart = GTime::seconds();
  }

  virtual ~MyRecurrentModel()
  {
  }

  virtual void onFinishedComputingStateEstimate(GData* pStateEstimate)
  {
    if(m_stateFilename)
      pStateEstimate->saveArff(m_stateFilename);
    cout << "% Computed state estimate in " << GTime::seconds() - m_dStart << " seconds.\n";
    cout.flush();
  }

  virtual void onObtainValidationScore(int timeSlice, double seconds, double squaredError)
  {
    if(m_validateInterval > 0)
      {
	if(squaredError == UNKNOWN_REAL_VALUE)
	  cout << (m_validateInterval * timeSlice) << ", ?\n";
	else
	  cout << (m_validateInterval * timeSlice) << ", " << sqrt(squaredError) << "\n";
	cout.flush();
      }
  }

};

template <class T>
class VectorContentsHolder
{
private:
  std::vector<T*>* m_p;

public:
  VectorContentsHolder(std::vector<T*>* p = NULL)
  {
    m_p = p;
  }

  VectorContentsHolder(const VectorContentsHolder& other)
  {
    ThrowError("tried to copy a holder");
  }

  ~VectorContentsHolder()
  {
    reset();
  }

  const VectorContentsHolder& operator=(const VectorContentsHolder& other)
  {
    ThrowError("tried to copy a holder");
    return *this;
  }

  void reset(std::vector<T*>* p = NULL)
  {
    if(p != m_p)
      {
	for(typename std::vector<T*>::iterator it = m_p->begin(); it != m_p->end(); it++)
	  delete(*it);
	m_p->clear();
	m_p = p;
      }
  }

  std::vector<T*>* get()
  {
    return m_p;
  }

  std::vector<T*>* release()
  {
    std::vector<T*>* pTmp = m_p;
    m_p = NULL;
    return pTmp;
  }

  std::vector<T*>& operator*() const
  {
    return *m_p;
  }

  std::vector<T*>* operator->() const
  {
    return m_p;
  }
};

void trainRecurrent(GArgReader& args)
{
  // Parse options
  unsigned int seed = getpid() * (unsigned int)time(NULL);
  vector<size_t> paramDims;
  const char* stateFilename = NULL;
  double validationInterval = 0;
  vector<string> validationFilenames;
  const char* outFilename = "model.twt";
  double trainTime = 60 * 60; // 1 hour
  bool useIsomap = false;
  while(args.next_is_flag())
    {
      if(args.if_pop("-seed"))
	seed = args.pop_uint();
      else if(args.if_pop("-paramdims"))
	{
	  unsigned int count = args.pop_uint();
	  for(unsigned int i = 0; i < count; i++)
	    paramDims.push_back(args.pop_uint());
	}
      else if(args.if_pop("-state"))
	stateFilename = args.pop_string();
      else if(args.if_pop("-validate"))
	{
	  validationInterval = args.pop_double();
	  int count = args.pop_uint();
	  for(int i = 0; i < count; i++)
	    {
	      validationFilenames.push_back(args.pop_string());
	      validationFilenames.push_back(args.pop_string());
	    }
	}
      else if(args.if_pop("-out"))
	outFilename = args.pop_string();
      else if(args.if_pop("-traintime"))
	trainTime = args.pop_double();
      else if(args.if_pop("-isomap"))
	useIsomap = true;
      else
	ThrowError("Invalid trainRecurrent option: ", args.peek());
    }

  // Parse the algorithm
  const char* alg = args.pop_string();
  int bpttDepth = 0;
  int bpttItersPerGrow = 0;
  double annealDeviation = 0.0;
  double annealDecay = 0.0;
  double annealTimeWindow = 0.0;
  if(strcmp(alg, "moses") == 0)
    {
    }
  else if(strcmp(alg, "aaron") == 0)
    {
    }
  else if(strcmp(alg, "joshua") == 0)
    {
    }
  else if(strcmp(alg, "bptt") == 0)
    {
      bpttDepth = args.pop_uint();
      bpttItersPerGrow = args.pop_uint();
    }
  else if(strcmp(alg, "bpttcal") == 0)
    {
      bpttDepth = args.pop_uint();
      bpttItersPerGrow = args.pop_uint();
    }
  else if(strcmp(alg, "evolutionary") == 0)
    {
    }
  else if(strcmp(alg, "hillclimber") == 0)
    {
    }
  else if(strcmp(alg, "annealing") == 0)
    {
      annealDeviation = args.pop_double();
      annealDecay = args.pop_double();
      annealTimeWindow = args.pop_double();
    }
  else
    ThrowError("Unrecognized recurrent model training algorithm: ", alg);

  // Load the data
  GData* pDataObs = GData::loadArff(args.pop_string());
  Holder<GData> hDataObs(pDataObs);
  GData* pDataAction = GData::loadArff(args.pop_string());
  Holder<GData> hDataAction(pDataAction);

  // Get the number of context dims
  int contextDims = args.pop_uint();

  // Infer remaining values and check that the parts fit together
  int pixels = 1;
  for(vector<size_t>::iterator it = paramDims.begin(); it != paramDims.end(); it++)
    pixels *= *it;
  int channels = pDataObs->cols() / pixels;
  if((channels * pixels) != pDataObs->cols())
    ThrowError("The number of columns in the observation data must be a multiple of the product of the param dims");

  // Instantiate the recurrent model
  GRand prng(seed);
  GTransducer* pTransitionFunc = InstantiateAlgorithm(&prng, args);
  Holder<GTransducer> hTransitionFunc(pTransitionFunc);
  if(!pTransitionFunc->canGeneralize())
    ThrowError("The algorithm specified for the transition function cannot be \"trained\". It can only be used to \"transduce\".");
  GTransducer* pObservationFunc = InstantiateAlgorithm(&prng, args);
  Holder<GTransducer> hObservationFunc(pObservationFunc);
  if(!pObservationFunc->canGeneralize())
    ThrowError("The algorithm specified for the observation function cannot be \"trained\". It can only be used to \"transduce\".");
  if(args.size() > 0)
    ThrowError("Superfluous argument: ", args.peek());
  MyRecurrentModel model((GSupervisedLearner*)hTransitionFunc.release(), (GSupervisedLearner*)hObservationFunc.release(), pDataAction->cols(), contextDims, pDataObs->cols(), &prng, &paramDims, stateFilename, validationInterval);

  // Set it up to do validation during training if specified
  vector<GData*> validationData;
  VectorContentsHolder<GData> hValidationData(&validationData);
  if(validationInterval > 0)
    {
      for(size_t i = 0; i < validationFilenames.size(); i++)
	validationData.push_back(GData::loadArff(validationFilenames[i].c_str()));
      model.validateDuringTraining(validationInterval, &validationData);
      cout << "@RELATION validation_scores\n\n@ATTRIBUTE seconds real\n@ATTRIBUTE " << alg << " real\n\n@DATA\n";
    }

  // Set other flags
  model.setTrainingSeconds(trainTime);
  model.setUseIsomap(useIsomap);

  // Do the training
  if(strcmp(alg, "moses") == 0)
    model.trainMoses(pDataAction, pDataObs);
  else if(strcmp(alg, "aaron") == 0)
    model.trainAaron(pDataAction, pDataObs);
  else if(strcmp(alg, "joshua") == 0)
    model.trainJoshua(pDataAction, pDataObs);
  else if(strcmp(alg, "bptt") == 0)
    model.trainBackPropThroughTime(pDataAction, pDataObs, bpttDepth, bpttItersPerGrow);
  else if(strcmp(alg, "evolutionary") == 0)
    model.trainEvolutionary(pDataAction, pDataObs);
  else if(strcmp(alg, "hillclimber") == 0)
    model.trainHillClimber(pDataAction, pDataObs, 0.0, 0.0, 0.0, true, false);
  else if(strcmp(alg, "annealing") == 0)
    model.trainHillClimber(pDataAction, pDataObs, annealDeviation, annealDecay, annealTimeWindow, false, true);
  GTwtDoc doc;
  doc.setRoot(model.toTwt(&doc));
  doc.save(outFilename);
}

int main(int argc, char *argv[])
{
#ifdef _DEBUG
  GApp::enableFloatingPointExceptions();
#endif
  int nRet = 0;
  try
    {
      PathData pd;
      GFile::parsePath(argv[0], &pd);
      const char* appName = argv[0] + pd.fileStart;

      GArgReader args(argc, argv);
      args.pop_string();
      if(args.size() >= 1)
	{
	  if(args.if_pop("usage"))
	    ShowUsage(appName);
	  else if(args.if_pop("train"))
	    Train(args);
	  else if(args.if_pop("test"))
	    Test(args);
	  else if(args.if_pop("predict"))
	    predict(args);
	  else if(args.if_pop("predictonepattern"))
	    predictOnePattern(args);
	  else if(args.if_pop("transduce"))
	    Transduce(args);
	  else if(args.if_pop("transacc"))
	    TransductiveAccuracy(args);
	  else if(args.if_pop("splittest"))
	    SplitTest(args);
	  else if(args.if_pop("crossvalidate"))
	    CrossValidate(args);
	  else if(args.if_pop("precisionrecall"))
	    PrecisionRecall(args);
	  else if(args.if_pop("trainsparse"))
	    trainSparse(args);
	  else if(args.if_pop("predictsparse"))
	    predictSparse(args);
	  else if(args.if_pop("trainrecurrent"))
	    trainRecurrent(args);
	  else
	    {
	      nRet = 1;
	      ThrowError("Unrecognized command: ", args.pop_string());
	    }
	}
      else
	{
	  nRet = 1;
	  ShowBriefUsage(appName);
	}
    }
  catch(const std::exception& e)
    {
      std::cerr << e.what() << "\n";
      nRet = 1;
    }
  return nRet;
}
