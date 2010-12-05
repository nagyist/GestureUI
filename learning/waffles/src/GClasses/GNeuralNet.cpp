/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#include "GNeuralNet.h"
#include "GMath.h"
#include "GActivation.h"
#include "GDistribution.h"
#include "GMacros.h"
#include "GRand.h"
#include "GVec.h"
#include "GTwt.h"
#include "GHillClimber.h"
#include "GTransform.h"
#include "GSparseMatrix.h"
#include "GImage.h"

namespace GClasses {

using std::vector;

void GNeuron::resetWeights(GRand* pRand, double inputCenter)
{
	for(vector<double>::iterator weight = m_weights.begin(); weight != m_weights.end(); weight++)
		*weight = pRand->normal() * 0.1;

	// Remove all bias (todo: this has very little effect since the weights are small already--why even bother?)
	double& bias = m_weights[0];
	for(vector<double>::iterator weight = m_weights.begin() + 1; weight != m_weights.end(); weight++)
		bias -= inputCenter * (*weight);
}

// ----------------------------------------------------------------------

void GNeuralNetLayer::resetWeights(GRand* pRand, double inputCenter)
{
	for(vector<GNeuron>::iterator neuron = m_neurons.begin(); neuron != m_neurons.end(); neuron++)
		neuron->resetWeights(pRand, inputCenter);
}

// ----------------------------------------------------------------------

GBackProp::GBackProp(GNeuralNet* pNN)
: m_pNN(pNN)
{
	// Initialize structures to mirror the neural network
	m_layers.resize(m_pNN->m_layers.size());
	for(size_t i = 0; i < m_layers.size(); i++)
	{
		GBackPropLayer& layer = m_layers[i];
		layer.m_neurons.resize(pNN->m_layers[i].m_neurons.size());
		for(size_t j = 0; j < layer.m_neurons.size(); j++)
		{
			GBackPropNeuron& neuron = layer.m_neurons[j];
			neuron.m_weights.resize(pNN->m_layers[i].m_neurons[j].m_weights.size());
		}
	}
}

void GBackProp::backPropLayer(GNeuralNetLayer* pNNFromLayer, GNeuralNetLayer* pNNToLayer, GBackPropLayer* pBPFromLayer, GBackPropLayer* pBPToLayer, size_t fromBegin)
{
	vector<GNeuron>::iterator nnFrom, nnCur;
	vector<GBackPropNeuron>::iterator bpFrom, bpCur;
	vector<double>::iterator nn_w;

	// Sum the error times weight for all the children
	nnFrom = pNNFromLayer->m_neurons.begin();
	bpFrom = pBPFromLayer->m_neurons.begin();
	nn_w = nnFrom->m_weights.begin() + 1 + fromBegin;
	bpCur = pBPToLayer->m_neurons.begin();
	while(bpCur != pBPToLayer->m_neurons.end())
	{
		bpCur->m_error = (*nn_w) * bpFrom->m_error; // use "=" for first pass
		nn_w++;
		bpCur++;
	}
	nnFrom++;
	bpFrom++;
	while(bpFrom != pBPFromLayer->m_neurons.end())
	{
		nn_w = nnFrom->m_weights.begin() + 1 + fromBegin;
		bpCur = pBPToLayer->m_neurons.begin();
		while(bpCur != pBPToLayer->m_neurons.end())
		{
			bpCur->m_error += (*nn_w) * bpFrom->m_error; // use "+=" for subsequent passes
			nn_w++;
			bpCur++;
		}
		nnFrom++;
		bpFrom++;
	}

	// Multiply by the derivative of the activation function
	nnCur = pNNToLayer->m_neurons.begin();
	bpCur = pBPToLayer->m_neurons.begin();
	while(bpCur != pBPToLayer->m_neurons.end())
	{
		bpCur->m_error *= pNNToLayer->m_pActivationFunction->derivativeInverse(nnCur->m_value);
		nnCur++;
		bpCur++;
	}
}

void GBackProp::backPropLayer2(GNeuralNetLayer* pNNFromLayer1, GNeuralNetLayer* pNNFromLayer2, GNeuralNetLayer* pNNToLayer, GBackPropLayer* pBPFromLayer1, GBackPropLayer* pBPFromLayer2, GBackPropLayer* pBPToLayer, int pass)
{
	vector<GNeuron>::iterator nnFrom, nnCur;
	vector<GBackPropNeuron>::iterator bpFrom, bpCur;
	vector<double>::iterator nn_w;

	double sum, alpha;
	int w = 1;
	nnCur = pNNToLayer->m_neurons.begin();
	bpCur = pBPToLayer->m_neurons.begin();
	while(bpCur != pBPToLayer->m_neurons.end())
	{
		// Sum error from the first previous layer
		sum = 0;
		nnFrom = pNNFromLayer1->m_neurons.begin();
		bpFrom = pBPFromLayer1->m_neurons.begin();
		while(bpFrom != pBPFromLayer1->m_neurons.end())
		{
			sum += nnFrom->m_weights[w] * bpFrom->m_error;
			nnFrom++;
			bpFrom++;
		}

		// Sum error from the second previous layer
		if(pNNFromLayer2)
		{
			nnFrom = pNNFromLayer2->m_neurons.begin();
			bpFrom = pBPFromLayer2->m_neurons.begin();
			while(bpFrom != pBPFromLayer2->m_neurons.end())
			{
				sum += nnFrom->m_weights[w] * bpFrom->m_error;
				nnFrom++;
				bpFrom++;
			}
		}

		// Multiply by derivative of activation function
		sum *= pNNToLayer->m_pActivationFunction->derivativeInverse(nnCur->m_value);

		// Average with error computed from previous passes
		alpha = 1.0 / pass;
		bpCur->m_error *= (1.0 - alpha);
		bpCur->m_error += (alpha * sum);

		nnCur++;
		bpCur++;
		w++;
	}
}

void GBackProp::adjustWeights(GNeuralNetLayer* pNNFromLayer, GNeuralNetLayer* pNNToLayer, GBackPropLayer* pBPFromLayer, double learningRate, double momentum)
{
	vector<GNeuron>::iterator nnFrom;
	vector<GBackPropNeuron>::iterator bpFrom;
	vector<double>::iterator nn_w;
	nnFrom = pNNFromLayer->m_neurons.begin();
	bpFrom = pBPFromLayer->m_neurons.begin();
	while(bpFrom != pBPFromLayer->m_neurons.end())
	{
		vector<GBackPropWeight>::iterator bp_w = bpFrom->m_weights.begin();
		vector<double>::iterator nn_w = nnFrom->m_weights.begin();
		bp_w->m_delta *= momentum;
		bp_w->m_delta += (learningRate * bpFrom->m_error);
		*nn_w = MAX(-1e12, MIN(1e12, *nn_w + bp_w->m_delta));
		bp_w++;
		nn_w++;
		for(vector<GNeuron>::iterator k = pNNToLayer->m_neurons.begin(); k != pNNToLayer->m_neurons.end(); k++)
		{
			bp_w->m_delta *= momentum;
			bp_w->m_delta += (learningRate * bpFrom->m_error * k->m_value);
			*nn_w = MAX(-1e12, MIN(1e12, *nn_w + bp_w->m_delta));
			bp_w++;
			nn_w++;
		}
		nnFrom++;
		bpFrom++;
	}
}

void GBackProp::adjustWeights(GNeuralNetLayer* pNNFromLayer, const double* pFeatures, GBackPropLayer* pBPFromLayer, double learningRate, double momentum)
{
	vector<GNeuron>::iterator nnFrom;
	vector<GBackPropNeuron>::iterator bpFrom;
	vector<double>::iterator nn_w;
	nnFrom = pNNFromLayer->m_neurons.begin();
	bpFrom = pBPFromLayer->m_neurons.begin();
	while(bpFrom != pBPFromLayer->m_neurons.end())
	{
		vector<GBackPropWeight>::iterator bp_w = bpFrom->m_weights.begin();
		vector<double>::iterator nn_w = nnFrom->m_weights.begin();
		bp_w->m_delta *= momentum;
		bp_w->m_delta += (learningRate * bpFrom->m_error);
		*nn_w = MAX(-1e12, MIN(1e12, *nn_w + bp_w->m_delta));
		bp_w++;
		nn_w++;
		for(const double* k = pFeatures; nn_w != nnFrom->m_weights.end(); k++)
		{
			if(*k != UNKNOWN_REAL_VALUE)
			{
				bp_w->m_delta *= momentum;
				bp_w->m_delta += (learningRate * bpFrom->m_error * (*k));
				*nn_w = MAX(-1e12, MIN(1e12, *nn_w + bp_w->m_delta));
			}
			bp_w++;
			nn_w++;
		}
		nnFrom++;
		bpFrom++;
	}
}

void GBackProp::backpropagate()
{
	size_t i = m_layers.size() - 1;
	GNeuralNetLayer* pNNPrevLayer = &m_pNN->m_layers[i];
	GBackPropLayer* pBPPrevLayer = &m_layers[i];
	for(i--; i < m_layers.size(); i--)
	{
		GNeuralNetLayer* pNNCurLayer = &m_pNN->m_layers[i];
		GBackPropLayer* pBPCurLayer = &m_layers[i];
		backPropLayer(pNNPrevLayer, pNNCurLayer, pBPPrevLayer, pBPCurLayer);
		pNNPrevLayer = pNNCurLayer;
		pBPPrevLayer = pBPCurLayer;
	}
}

void GBackProp::descendGradient(const double* pFeatures, double learningRate, double momentum)
{
	size_t i = m_layers.size() - 1;
	GNeuralNetLayer* pNNPrevLayer = &m_pNN->m_layers[i];
	GBackPropLayer* pBPPrevLayer = &m_layers[i];
	for(i--; i < m_layers.size(); i--)
	{
		GNeuralNetLayer* pNNCurLayer = &m_pNN->m_layers[i];
		GBackPropLayer* pBPCurLayer = &m_layers[i];
		adjustWeights(pNNPrevLayer, pNNCurLayer, pBPPrevLayer, learningRate, momentum);
		pNNPrevLayer = pNNCurLayer;
		pBPPrevLayer = pBPCurLayer;
	}

	// adjust the weights on the last hidden layer
	adjustWeights(pNNPrevLayer, pFeatures, pBPPrevLayer, m_pNN->learningRate(), m_pNN->momentum());
}

void GBackProp::adjustFeatures(double* pFeatures, double learningRate, size_t skip)
{
	GNeuralNetLayer& nnLayer = m_pNN->m_layers[0];
	GBackPropLayer& bpLayer = m_layers[0];
	vector<GNeuron>::iterator nn = nnLayer.m_neurons.begin();
	vector<GBackPropNeuron>::iterator bp = bpLayer.m_neurons.begin();
	while(nn != nnLayer.m_neurons.end())
	{
		double* pOut = pFeatures;
		for(vector<double>::iterator w = nn->m_weights.begin() + 1 + skip; w != nn->m_weights.end(); w++)
			*(pOut++) += learningRate * bp->m_error * (*w);
		nn++;
		bp++;
	}
}


// ----------------------------------------------------------------------

GNeuralNet::GNeuralNet(GRand* pRand)
: GIncrementalLearner(), m_pRand(pRand)
{
	m_featureDims = 0;
	m_labelDims = 0;
	m_pBackProp = NULL;
	m_learningRate = 0.1;
	m_momentum = 0.0;
	m_minImprovement = 0.002;
	m_epochsPerValidationCheck = 200;
	m_validationPortion = 0;
	m_backPropTargetFunction = squared_error;
	m_pActivationFunction = NULL;
	m_layers.resize(1);
}

GNeuralNet::GNeuralNet(GTwtNode* pNode, GRand* pRand)
 : GIncrementalLearner(pNode), m_pRand(pRand)
{
	// Create the layers
	m_pActivationFunction = NULL;
	m_featureDims = 0;
	m_labelDims = 0;
	m_layers.resize(1);
	m_pBackProp = NULL;
	m_featureDims = (int)pNode->field("featureDims")->asInt();
	GTwtNode* pLayerList = pNode->field("layers");
	size_t layerCount = pLayerList->itemCount();
	for(size_t i = 0; i < layerCount - 1; i++)
	{
		GTwtNode* pActivation = pLayerList->item(i)->fieldIfExists("af");
		if(pActivation)
			setActivationFunction(GActivationFunction::fromTwt(pActivation), true);
		else if(i == 0)
			ThrowError("The first layer is expected to specify an activation function");
		addLayer((int)pLayerList->item(i)->field("nodes")->asInt());
	}
	GTwtNode* pActivation = pLayerList->item(layerCount - 1)->fieldIfExists("af");
	if(pActivation)
		setActivationFunction(GActivationFunction::fromTwt(pActivation), true);
	else if(layerCount == 1)
		ThrowError("The first layer is expected to specify an activation function");
	m_labelDims = (int)pLayerList->item(layerCount - 1)->field("nodes")->asInt();

	// Enable training
	sp_relation pRel;
	pRel = new GUniformRelation(m_featureDims + m_labelDims, 0);
	enableIncrementalLearning(pRel, m_labelDims, NULL, NULL);

	// Set other settings
	m_learningRate = pNode->field("learningRate")->asDouble();
	m_momentum = pNode->field("momentum")->asDouble();
	m_backPropTargetFunction = (TargetFunction)pNode->field("target")->asInt();

	// Set the weights
	GTwtNode* pWeightList = pNode->field("weights");
	size_t wc = pWeightList->itemCount();
	if(wc != countWeights())
		ThrowError("Weights don't line up. (expected ", gformat(countWeights()), ", got ", gformat(wc), ".)");
	GTEMPBUF(double, pWeights, wc);
	for(size_t i = 0; i < wc; i++)
		pWeights[i] = pWeightList->item(i)->asDouble();
	setWeights(pWeights);
}

GNeuralNet::~GNeuralNet()
{
	releaseTrainingJunk();
	for(vector<GActivationFunction*>::iterator it = m_activationFunctions.begin(); it != m_activationFunctions.end(); it++)
		delete(*it);
}

// virtual
GTwtNode* GNeuralNet::toTwt(GTwtDoc* pDoc)
{
	if(m_labelDims == 0)
		ThrowError("The network has not been trained");
	GTwtNode* pNode = baseTwtNode(pDoc, "GNeuralNet");

	// Add the layer sizes
	pNode->addField(pDoc, "featureDims", pDoc->newInt(m_featureDims));
	GTwtNode* pLayerList = pNode->addField(pDoc, "layers", pDoc->newList(m_layers.size()));
	GActivationFunction* pPrevSF = NULL;
	for(size_t i = 0; i < m_layers.size(); i++)
	{
		GTwtNode* pLayerObj = pLayerList->setItem(i, pDoc->newObj());
		pLayerObj->addField(pDoc, "nodes", pDoc->newInt(m_layers[i].m_neurons.size()));
		GAssert(i != m_layers.size() - 1 || m_layers[i].m_neurons.size() == (size_t)m_labelDims);
		if(m_layers[i].m_pActivationFunction != pPrevSF)
		{
			pPrevSF = m_layers[i].m_pActivationFunction;
			pLayerObj->addField(pDoc, "af", m_layers[i].m_pActivationFunction->toTwt(pDoc));
		}
	}

	// Add other settings
	pNode->addField(pDoc, "learningRate", pDoc->newDouble(m_learningRate));
	pNode->addField(pDoc, "momentum", pDoc->newDouble(m_momentum));
	pNode->addField(pDoc, "target", pDoc->newInt(m_backPropTargetFunction));

	// Add the weights
	{
		int wc = countWeights();
		GTEMPBUF(double, pWeights, wc);
		weights(pWeights);
		GTwtNode* pWeightList = pNode->addField(pDoc, "weights", pDoc->newList(wc));
		for(int i = 0; i < wc; i++)
			pWeightList->setItem(i, pDoc->newDouble(pWeights[i]));
	}

	return pNode;
}

void GNeuralNet::setActivationFunction(GActivationFunction* pSF, bool hold)
{
	m_pActivationFunction = pSF;
	if(hold)
		m_activationFunctions.push_back(pSF);
}

// virtual
int GNeuralNet::featureDims()
{
	if(m_labelDims < 1)
		ThrowError("not yet trained");
	return m_featureDims;
}

// virtual
int GNeuralNet::labelDims()
{
	if(m_labelDims < 1)
		ThrowError("not yet trained");
	return m_labelDims;
}

void GNeuralNet::releaseTrainingJunk()
{
	delete(m_pBackProp);
	m_pBackProp = NULL;
}

void GNeuralNet::addLayer(int nodeCount)
{
	if(m_labelDims != 0)
		ThrowError("Changing the network structure after some training has begun is not yet supported.");
	if(nodeCount < 1)
		ThrowError("Cannot add a layer with fewer than 1 node");

	// Add a new layer to be the new output layer
	size_t i = m_layers.size();
	GAssert(i > 0); // There should already be an output layer
	m_layers.resize(i + 1);

	// Turn the old output layer into the new hidden layer
	GNeuralNetLayer& newLayer = m_layers[i - 1];
	newLayer.m_neurons.resize(nodeCount);
	if(!m_pActivationFunction)
		setActivationFunction(new GActivationLogistic(), true);
	newLayer.m_pActivationFunction = m_pActivationFunction;

	// Give each node in the previous layer a weight for the bias, plus a weight for each node in this layer
	if(i > 1)
	{
		size_t weightCount = 1 + m_layers[i - 2].m_neurons.size(); // bias, plus a connection to each previous node
		for(size_t j = 0; j < newLayer.m_neurons.size(); j++)
			newLayer.m_neurons[j].m_weights.resize(weightCount);
	}
}

void GNeuralNet::addNode(size_t layer)
{
	if(layer >= m_layers.size())
		ThrowError("layer index out of range");

	// Add a new neuron to this layer
	GNeuralNetLayer& l = m_layers[layer];
	size_t n = l.m_neurons.size();
	l.m_neurons.resize(n + 1);
	GNeuron& neuron = l.m_neurons[n];
	neuron.m_weights.resize(l.m_neurons[0].m_weights.size());
	neuron.resetWeights(m_pRand, l.m_pActivationFunction->center());

	// Add another weight to each node in the next layer
	if(layer < m_layers.size() - 1)
	{
		GNeuralNetLayer& layerNext = m_layers[layer + 1];
		for(vector<GNeuron>::iterator it = layerNext.m_neurons.begin(); it != layerNext.m_neurons.end(); it++)
			it->m_weights.push_back(0.05 * m_pRand->normal());
	}
}

void GNeuralNet::dropNode(size_t layer, size_t node)
{
	if(layer >= m_layers.size())
		ThrowError("layer index out of range");
	GNeuralNetLayer& l = m_layers[layer];
	if(node >= l.m_neurons.size())
		ThrowError("node index out of range");
	if(l.m_neurons.size() == 1)
		ThrowError("The layer must have at least one node in it");

	// Drop the neuron from this layer
	l.m_neurons.erase(l.m_neurons.begin() + node);

	// Remove the corresponding weight from each node in the next layer
	if(layer < m_layers.size() - 1)
	{
		GNeuralNetLayer& layerNext = m_layers[layer + 1];
		for(vector<GNeuron>::iterator it = layerNext.m_neurons.begin(); it != layerNext.m_neurons.end(); it++)
			it->m_weights.erase(it->m_weights.begin() + node + 1);
	}
}

size_t GNeuralNet::countWeights()
{
	if(m_labelDims == 0)
		ThrowError("enableIncrementalLearning must be called before this method");
	size_t wc = 0;
	for(vector<GNeuralNetLayer>::iterator layer = m_layers.begin(); layer != m_layers.end(); layer++)
		wc += layer->m_neurons.size() * layer->m_neurons.begin()->m_weights.size(); // We assume that every node in a layer has the same number of weights
	return wc;
}

void GNeuralNet::weights(double* pOutWeights)
{
	if(m_labelDims == 0)
		ThrowError("enableIncrementalLearning must be called before this method");
	for(vector<GNeuralNetLayer>::iterator layer = m_layers.begin(); layer != m_layers.end(); layer++)
	{
		for(vector<GNeuron>::iterator neuron = layer->m_neurons.begin(); neuron != layer->m_neurons.end(); neuron++)
		{
			for(vector<double>::iterator weight = neuron->m_weights.begin(); weight != neuron->m_weights.end(); weight++)
			{
				*pOutWeights = *weight;
				pOutWeights++;
			}
		}
	}
}

void GNeuralNet::setWeights(const double* pWeights)
{
	if(m_labelDims == 0)
		ThrowError("enableIncrementalLearning must be called before this method");
	for(vector<GNeuralNetLayer>::iterator layer = m_layers.begin(); layer != m_layers.end(); layer++)
	{
		for(vector<GNeuron>::iterator neuron = layer->m_neurons.begin(); neuron != layer->m_neurons.end(); neuron++)
		{
			for(vector<double>::iterator weight = neuron->m_weights.begin(); weight != neuron->m_weights.end(); weight++)
			{
				*weight = *pWeights;
				pWeights++;
			}
		}
	}
}

void GNeuralNet::copyWeights(GNeuralNet* pOther)
{
	if(m_labelDims == 0 || pOther->m_labelDims == 0)
		ThrowError("enableIncrementalLearning must be called on both networks before this method");
	GAssert(m_layers.size() == pOther->m_layers.size());
	vector<GNeuralNetLayer>::iterator layerOther = pOther->m_layers.begin();
	for(vector<GNeuralNetLayer>::iterator layer = m_layers.begin(); layer != m_layers.end(); layer++)
	{
		GAssert(layer->m_neurons.size() == layerOther->m_neurons.size());
		vector<GNeuron>::iterator neuronOther = layerOther->m_neurons.begin();
		for(vector<GNeuron>::iterator neuron = layer->m_neurons.begin(); neuron != layer->m_neurons.end(); neuron++)
		{
			GAssert(neuron->m_weights.size() == neuronOther->m_weights.size());
			vector<double>::iterator weightOther = neuronOther->m_weights.begin();
			for(vector<double>::iterator weight = neuron->m_weights.begin(); weight != neuron->m_weights.end(); weight++)
				*weight = *(weightOther++);
			neuronOther++;
		}
		layerOther++;
	}
}

void GNeuralNet::copyStructure(GNeuralNet* pOther)
{
	if(pOther->m_labelDims == 0)
		ThrowError("enableIncrementalLearning must be called before this method");
	releaseTrainingJunk();
	for(vector<GActivationFunction*>::iterator it = m_activationFunctions.begin(); it != m_activationFunctions.end(); it++)
		delete(*it);
	m_pActivationFunction = NULL;
	m_layers.resize(pOther->m_layers.size());
	for(size_t i = 0; i < m_layers.size(); i++)
	{
		if(pOther->m_layers[i].m_pActivationFunction != m_pActivationFunction)
		{
			setActivationFunction(pOther->m_layers[i].m_pActivationFunction->clone(), true);
			m_layers[i].m_pActivationFunction = m_pActivationFunction;
			setActivationFunction(pOther->m_layers[i].m_pActivationFunction, false);
		}
		else
			m_layers[i].m_pActivationFunction = m_layers[i - 1].m_pActivationFunction;
		m_layers[i].m_neurons.resize(pOther->m_layers[i].m_neurons.size());
		for(size_t j = 0; j < m_layers[i].m_neurons.size(); j++)
			m_layers[i].m_neurons[j].m_weights.resize(pOther->m_layers[i].m_neurons[j].m_weights.size());
	}
	setActivationFunction(m_layers[m_layers.size() - 1].m_pActivationFunction, false);
	m_featureDims = pOther->m_featureDims;
	m_labelDims = pOther->m_labelDims;
	m_learningRate = pOther->m_learningRate;
	m_momentum = pOther->m_momentum;
	m_validationPortion = pOther->m_validationPortion;
	m_minImprovement = pOther->m_minImprovement;
	m_epochsPerValidationCheck = pOther->m_epochsPerValidationCheck;
	m_backPropTargetFunction = pOther->m_backPropTargetFunction;
	if(pOther->m_pBackProp)
		m_pBackProp = new GBackProp(this);
}

void GNeuralNet::perturbAllWeights(double deviation)
{
	if(m_labelDims == 0)
		ThrowError("enableIncrementalLearning must be called before this method");
	for(vector<GNeuralNetLayer>::iterator layer = m_layers.begin(); layer != m_layers.end(); layer++)
	{
		for(vector<GNeuron>::iterator neuron = layer->m_neurons.begin(); neuron != layer->m_neurons.end(); neuron++)
			for(vector<double>::iterator weight = neuron->m_weights.begin(); weight != neuron->m_weights.end(); weight++)
				(*weight) += (m_pRand->normal() * deviation);
	}
}

void GNeuralNet::clipWeights(double max)
{
	if(m_labelDims == 0)
		ThrowError("enableIncrementalLearning must be called before this method");
	for(vector<GNeuralNetLayer>::iterator layer = m_layers.begin(); layer != m_layers.end(); layer++)
	{
		for(vector<GNeuron>::iterator neuron = layer->m_neurons.begin(); neuron != layer->m_neurons.end(); neuron++)
		{
			for(vector<double>::iterator weight = neuron->m_weights.begin() + 1; weight != neuron->m_weights.end(); weight++)
				(*weight) = MAX(-max, MIN(max, *weight));
		}
	}
}

void GNeuralNet::decayWeights(double lambda, double gamma)
{
	if(m_labelDims == 0)
		ThrowError("enableIncrementalLearning must be called before this method");
	for(vector<GNeuralNetLayer>::iterator layer = m_layers.begin(); layer != m_layers.end(); layer++)
	{
		double d = (1.0 - lambda * m_learningRate);
		for(vector<GNeuron>::iterator neuron = layer->m_neurons.begin(); neuron != layer->m_neurons.end(); neuron++)
		{
			for(vector<double>::iterator weight = neuron->m_weights.begin() + 1; weight != neuron->m_weights.end(); weight++)
				*weight *= d;
		}
		lambda *= gamma;
	}
}

void GNeuralNet::forwardProp(const double* pRow)
{
	// Propagate from the feature vector to the first layer
	vector<GNeuralNetLayer>::iterator pLayer = m_layers.begin();
	double net;
	for(vector<GNeuron>::iterator i = pLayer->m_neurons.begin(); i != pLayer->m_neurons.end(); i++)
	{
		vector<double>::iterator j = i->m_weights.begin();
		net = *(j++); // (the first weight is the bias)
		const double* pR = pRow;
		while(j != i->m_weights.end())
			net += *(j++) * *(pR++);
		i->m_value = pLayer->m_pActivationFunction->squash(net);
	}

	// Do the rest of the hidden layers
	vector<GNeuralNetLayer>::iterator pPrevLayer = pLayer;
	for(pLayer++; pLayer != m_layers.end(); pLayer++)
	{
		for(vector<GNeuron>::iterator i = pLayer->m_neurons.begin(); i != pLayer->m_neurons.end(); i++)
		{
			vector<double>::iterator j = i->m_weights.begin();
			net = *(j++); // (the first weight is the bias)
			vector<GNeuron>::iterator k = pPrevLayer->m_neurons.begin();
			while(k != pPrevLayer->m_neurons.end())
			{
				net += *(j++) * k->m_value;
				k++;
			}
			i->m_value = pLayer->m_pActivationFunction->squash(net);
		}
		pPrevLayer = pLayer;
	}
}

// virtual
void GNeuralNet::predictDistribution(const double* pIn, GPrediction* pOut)
{
	if(m_labelDims == 0)
		ThrowError("enableIncrementalLearning must be called before this method");

	// Do the evaluation
	forwardProp(pIn);

	// Convert outputs to external data
	GNeuralNetLayer& outputLayer = m_layers[m_layers.size() - 1];
	for(vector<GNeuron>::iterator i = outputLayer.m_neurons.begin(); i != outputLayer.m_neurons.end(); i++)
	{
		GNormalDistribution* pNorm = pOut->makeNormal();
		pNorm->setMeanAndVariance(i->m_value, 1.0);
		pOut++;
	}
}

void GNeuralNet::copyPrediction(double* pOut)
{
	GNeuralNetLayer& outputLayer = m_layers[m_layers.size() - 1];
	for(vector<GNeuron>::iterator i = outputLayer.m_neurons.begin(); i != outputLayer.m_neurons.end(); i++)
		*(pOut++) = i->m_value;
}

double GNeuralNet::sumSquaredPredictionError(const double* pTarget)
{
	GNeuralNetLayer& outputLayer = m_layers[m_layers.size() - 1];
	double sse = 0.0;
	for(vector<GNeuron>::iterator i = outputLayer.m_neurons.begin(); i != outputLayer.m_neurons.end(); i++)
	{
		double d = *(pTarget++) - i->m_value;
		sse += (d * d);
	}
	return sse;
}

// virtual
void GNeuralNet::predict(const double* pIn, double* pOut)
{
	if(m_labelDims == 0)
		ThrowError("enableIncrementalLearning must be called before this method");
	forwardProp(pIn);
	copyPrediction(pOut);
}

void GNeuralNet::train(GData* pData, int labelDims)
{
	int validationRows = (int)(m_validationPortion * pData->rows());
	if(validationRows > 0)
	{
		GData dataValidate(pData->relation(), pData->heap());
		dataValidate.reserve(validationRows);
		pData->splitBySize(&dataValidate, validationRows);
		trainWithValidation(pData, &dataValidate, labelDims);
		pData->mergeVert(&dataValidate);
	}
	else
		trainWithValidation(pData, pData, labelDims);
}

// virtual
void GNeuralNet::trainSparse(GSparseMatrix* pData, int labelDims)
{
	sp_relation pRel = new GUniformRelation(0, pData->cols());
	enableIncrementalLearning(pRel, labelDims, NULL, NULL);

	GTEMPBUF(size_t, indexes, pData->rows());
	GIndexVec::makeIndexVec(indexes, pData->rows());
	GTEMPBUF(double, pFullRow, pData->cols());
	for(int epochs = 0; epochs < 50; epochs++) // todo: need a better stopping criterion
	{
		GIndexVec::shuffle(indexes, pData->rows(), m_pRand);
		for(unsigned int i = 0; i < pData->rows(); i++)
		{
			pData->fullRow(pFullRow, indexes[i]);
			forwardProp(pFullRow);
			setErrorOnOutputLayer(pFullRow + pData->cols() - labelDims, m_backPropTargetFunction);
			m_pBackProp->backpropagate();
			m_pBackProp->descendGradient(pFullRow, m_learningRate, m_momentum);
		}
	}
}

double GNeuralNet::validationSquaredError(GData* pValidationData)
{
	int i, nIndex;
	double* pRow;
	double d;
	double dError = 0;
	GNeuralNetLayer& outputLayer = m_layers[m_layers.size() - 1];
	size_t nCount = pValidationData->rows();
	for(size_t n = 0; n < nCount; n++)
	{
		pRow = pValidationData->row(n);
		forwardProp(pRow);
		nIndex = m_featureDims;
		for(i = 0; i < m_labelDims; i++)
		{
			d = pRow[nIndex++] - outputLayer.m_neurons[i].m_value;
			d *= d;
			dError += d;
		}
	}
	return dError;
}

void GNeuralNet::trainEpoch(GData* pTrainingData)
{
	pTrainingData->shuffle(m_pRand);
	size_t nRowCount = pTrainingData->rows();
	for(size_t n = 0; n < nRowCount; n++)
	{
		const double* pRow = pTrainingData->row(n);
		forwardProp(pRow);
		setErrorOnOutputLayer(pRow + m_featureDims, m_backPropTargetFunction);
		m_pBackProp->backpropagate();
		m_pBackProp->descendGradient(pRow, m_learningRate, m_momentum);
	}
}

int GNeuralNet::trainWithValidation(GData* pTrainingData, GData* pValidationData, int labelDims)
{
	enableIncrementalLearning(pTrainingData->relation(), labelDims, NULL, NULL);

	// Do the epochs
	int nEpochs;
	double dBestError = 1e308;
	int nEpochsSinceValidationCheck = 0;
	double dSumSquaredError;
	for(nEpochs = 0; true; nEpochs++)
	{
		// Do an epoch
		trainEpoch(pTrainingData);

		// Check for termination condition
		if(nEpochsSinceValidationCheck >= m_epochsPerValidationCheck)
		{
			nEpochsSinceValidationCheck = 0;
			dSumSquaredError = validationSquaredError(pValidationData);
			if(1.0 - dSumSquaredError / dBestError < m_minImprovement)
				break;
			if(dSumSquaredError < dBestError)
				dBestError = dSumSquaredError;
		}
		else
			nEpochsSinceValidationCheck++;
	}

	releaseTrainingJunk();
	return nEpochs;
}

// virtual
void GNeuralNet::enableIncrementalLearning(sp_relation& pRelation, int labelDims, double* pMins, double* pRanges)
{
	if(labelDims < 1)
		ThrowError("labelDims must be at least 1");
	if(pRelation->size() < labelDims)
		ThrowError("fewer total attributes than label dims");
	if(!pRelation->areContinuous(0, pRelation->size()))
		ThrowError("Only continuous values are supported.");

	// Adjust the size of the output layer
	m_labelDims = labelDims;
	m_featureDims = pRelation->size() - m_labelDims;
	GNeuralNetLayer& layerOut = m_layers[m_layers.size() - 1];
	layerOut.m_neurons.resize(labelDims);
	if(!m_pActivationFunction)
		setActivationFunction(new GActivationLogistic(), true);
	layerOut.m_pActivationFunction = m_pActivationFunction;
	if(m_layers.size() > 1)
	{
		// Establish the number of weights on the output layer
		size_t weightCount = 1 + m_layers[m_layers.size() - 2].m_neurons.size();
		for(int i = 0; i < labelDims; i++)
			layerOut.m_neurons[i].m_weights.resize(weightCount);
	}

	// Establish the number of weights on the first layer
	GNeuralNetLayer& layerIn = m_layers[0];
	for(size_t i = 0; i < layerIn.m_neurons.size(); i++)
		layerIn.m_neurons[i].m_weights.resize(1 + m_featureDims);

	// Initialize the weights with small random values
	double inputCenter = 0.5; // Assume inputs have a range from 0 to 1. If this not correct, learning may be slightly slower.
	for(size_t i = 0; i < m_layers.size(); i++)
	{
		m_layers[i].resetWeights(m_pRand, inputCenter);
		inputCenter = m_layers[i].m_pActivationFunction->center();
	}

	// Make the training junk
	releaseTrainingJunk();
	m_pBackProp = new GBackProp(this);
}

// virtual
void GNeuralNet::trainIncremental(const double* pIn, const double* pOut)
{
	if(m_labelDims == 0)
		ThrowError("enableIncrementalLearning must be called before this method");
	forwardProp(pIn);
	setErrorOnOutputLayer(pOut, m_backPropTargetFunction);
	m_pBackProp->backpropagate();
	m_pBackProp->descendGradient(pIn, m_learningRate, m_momentum);
}

void GNeuralNet::setErrorOnOutputLayer(const double* pTarget, TargetFunction eTargetFunction)
{
	// Compute error on output layer
	GBackPropLayer& bpOutputLayer = m_pBackProp->layer(m_layers.size() - 1);
	GNeuralNetLayer& nnOutputLayer = m_layers[m_layers.size() - 1];
	switch(eTargetFunction)
	{
		case squared_error:
			for(size_t j = 0; j < bpOutputLayer.m_neurons.size(); j++)
			{
				if(*pTarget == UNKNOWN_REAL_VALUE)
					bpOutputLayer.m_neurons[j].m_error = 0.0;
				else
					bpOutputLayer.m_neurons[j].m_error = (*pTarget - nnOutputLayer.m_neurons[j].m_value) * nnOutputLayer.m_pActivationFunction->derivativeInverse(nnOutputLayer.m_neurons[j].m_value);
				pTarget++;
			}
			break;

		case cross_entropy:
			for(size_t j = 0; j < bpOutputLayer.m_neurons.size(); j++)
			{
				if(*pTarget == UNKNOWN_REAL_VALUE)
					bpOutputLayer.m_neurons[j].m_error = 0.0;
				else
					bpOutputLayer.m_neurons[j].m_error = *pTarget - nnOutputLayer.m_neurons[j].m_value;
				pTarget++;
			}
			break;

		case root_cubed_error:
			for(size_t j = 0; j < bpOutputLayer.m_neurons.size(); j++)
			{
				if(*pTarget == UNKNOWN_REAL_VALUE)
					bpOutputLayer.m_neurons[j].m_error = 0.0;
				else
				{
					bpOutputLayer.m_neurons[j].m_error = GMath::signedRoot(*pTarget - nnOutputLayer.m_neurons[j].m_value) * nnOutputLayer.m_pActivationFunction->derivativeInverse(nnOutputLayer.m_neurons[j].m_value);
				}
				pTarget++;
			}
			break;

		default:
			ThrowError("Unrecognized target function for back-propagation");
			break;
	}
}

#ifndef NO_TEST_CODE
const char g_neural_net_test_data[] =
	"@RELATION Y\n"
	"@ATTRIBUTE x continuous\n"
	"@ATTRIBUTE y continuous\n"
	"@ATTRIBUTE g continuous\n"
	"@DATA\n"
	"0,-0.7,1\n";

void GNeuralNet_testMath()
{
	// Parse the data
	GData* pData = GData::parseArff(g_neural_net_test_data, sizeof(g_neural_net_test_data) - 1);
	Holder<GData> hData(pData);

	// Make the Neural Network
	GRand prng(0);
	GNeuralNet nn(&prng);
	nn.setLearningRate(0.175);
	nn.setMomentum(0.9);
	nn.addLayer(3);
	nn.enableIncrementalLearning(pData->relation(), 1, NULL, NULL);
	if(nn.countWeights() != 13)
		ThrowError("Wrong number of weights");
	GNeuralNetLayer& layerOut = nn.getLayer(1);
	layerOut.m_neurons[0].m_weights[0] = 0.02; // w_0
	layerOut.m_neurons[0].m_weights[1] = -0.01; // w_1
	layerOut.m_neurons[0].m_weights[2] = 0.03; // w_2
	layerOut.m_neurons[0].m_weights[3] = 0.02; // w_3
	GNeuralNetLayer& layerHidden = nn.getLayer(0);
	layerHidden.m_neurons[0].m_weights[0] = -0.01; // w_4
	layerHidden.m_neurons[0].m_weights[1] = -0.03; // w_5
	layerHidden.m_neurons[0].m_weights[2] = 0.03; // w_6
	layerHidden.m_neurons[1].m_weights[0] = 0.01; // w_7
	layerHidden.m_neurons[1].m_weights[1] = 0.04; // w_8
	layerHidden.m_neurons[1].m_weights[2] = -0.02; // w_9
	layerHidden.m_neurons[2].m_weights[0] = -0.02; // w_10
	layerHidden.m_neurons[2].m_weights[1] = 0.03; // w_11
	layerHidden.m_neurons[2].m_weights[2] = 0.02; // w_12

	bool useCrossEntropy = false;

	// Test forward prop
	double tol = 1e-12;
	double pat[3];
	GVec::copy(pat, pData->row(0), 2);
	nn.predict(pat, pat + 2);
	// Here is the math (done by hand) for why these results are expected:
	// Row: {0, -0.7, 1}
	// o_1 = squash(w_4*1+w_5*x+w_6*y) = 1/(1+exp(-(-.01*1-.03*0+.03*(-.7)))) = 0.4922506205862
	// o_2 = squash(w_7*1+w_8*x+w_9*y) = 1/(1+exp(-(.01*1+.04*0-.02*(-.7)))) = 0.50599971201659
	// o_3 = squash(w_10*1+w_11*x+w_12*y) = 1/(1+exp(-(-.02*1+.03*0+.02*(-.7)))) = 0.49150081873869
	// o_0 = squash(w_0*1+w_1*o_1+w_2*o_2+w_3*o_3) = 1/(1+exp(-(.02*1-.01*.4922506205862+.03*.50599971201659+.02*.49150081873869))) = 0.51002053349535
	if(ABS(pat[2] - 0.51002053349535) > tol) ThrowError("forward prop problem");

	// Test that the output error is computed properly
	nn.trainIncremental(pData->row(0), pData->row(0) + 2);
	GBackProp* pBP = nn.backProp();
	// Here is the math (done by hand) for why these results are expected:
	// e_0 = output*(1-output)*(target-output) = .51002053349535*(1-.51002053349535)*(1-.51002053349535) = 0.1224456672531
	if(useCrossEntropy)
	{
		// Here is the math for why these results are expected:
		// e_0 = target-output = 1-.51002053349535 = 0.4899794665046473
		if(ABS(pBP->layer(1).m_neurons[0].m_error - 0.4899794665046473) > tol) ThrowError("problem computing output error");
	}
	else
	{
		// Here is the math for why these results are expected:
		// e_0 = output*(1-output)*(target-output) = .51002053349535*(1-.51002053349535)*(1-.51002053349535) = 0.1224456672531
		if(ABS(pBP->layer(1).m_neurons[0].m_error - 0.1224456672531) > tol) ThrowError("problem computing output error");
	}

	// Test Back Prop
	if(useCrossEntropy)
	{
		if(ABS(pBP->layer(0).m_neurons[0].m_error + 0.0012246544194742083) > tol) ThrowError("back prop problem");
		// e_2 = o_2*(1-o_2)*(w_2*e_0) = 0.00091821027577176
		if(ABS(pBP->layer(0).m_neurons[1].m_error - 0.0036743168717579557) > tol) ThrowError("back prop problem");
		// e_3 = o_3*(1-o_3)*(w_3*e_0) = 0.00061205143636003
		if(ABS(pBP->layer(0).m_neurons[2].m_error - 0.002449189448583718) > tol) ThrowError("back prop problem");
	}
	else
	{
		// e_1 = o_1*(1-o_1)*(w_1*e_0) = .4922506205862*(1-.4922506205862)*(-.01*.1224456672531) = -0.00030604063598154
		if(ABS(pBP->layer(0).m_neurons[0].m_error + 0.00030604063598154) > tol) ThrowError("back prop problem");
		// e_2 = o_2*(1-o_2)*(w_2*e_0) = 0.00091821027577176
		if(ABS(pBP->layer(0).m_neurons[1].m_error - 0.00091821027577176) > tol) ThrowError("back prop problem");
		// e_3 = o_3*(1-o_3)*(w_3*e_0) = 0.00061205143636003
		if(ABS(pBP->layer(0).m_neurons[2].m_error - 0.00061205143636003) > tol) ThrowError("back prop problem");
	}

	// Test weight update
	if(useCrossEntropy)
	{
		if(ABS(layerOut.m_neurons[0].m_weights[0] - 0.10574640663831328) > tol) ThrowError("weight update problem");
		if(ABS(layerOut.m_neurons[0].m_weights[1] - 0.032208721880745944) > tol) ThrowError("weight update problem");
	}
	else
	{
		// d_0 = (d_0*momentum)+(learning_rate*e_0*1) = 0*.9+.175*.1224456672531*1
		// w_0 = w_0 + d_0 = .02+.0214279917693 = 0.041427991769293
		if(ABS(layerOut.m_neurons[0].m_weights[0] - 0.041427991769293) > tol) ThrowError("weight update problem");
		// d_1 = (d_1*momentum)+(learning_rate*e_0*o_1) = 0*.9+.175*.1224456672531*.4922506205862
		// w_1 = w_1 + d_1 = -.01+.0105479422563 = 0.00054794224635029
		if(ABS(layerOut.m_neurons[0].m_weights[1] - 0.00054794224635029) > tol) ThrowError("weight update problem");
		if(ABS(layerOut.m_neurons[0].m_weights[2] - 0.040842557664356) > tol) ThrowError("weight update problem");
		if(ABS(layerOut.m_neurons[0].m_weights[3] - 0.030531875498533) > tol) ThrowError("weight update problem");
		if(ABS(layerHidden.m_neurons[0].m_weights[0] + 0.010053557111297) > tol) ThrowError("weight update problem");
		if(ABS(layerHidden.m_neurons[0].m_weights[1] + 0.03) > tol) ThrowError("weight update problem");
		if(ABS(layerHidden.m_neurons[0].m_weights[2] - 0.030037489977908) > tol) ThrowError("weight update problem");
		if(ABS(layerHidden.m_neurons[1].m_weights[0] - 0.01016068679826) > tol) ThrowError("weight update problem");
		if(ABS(layerHidden.m_neurons[1].m_weights[1] - 0.04) > tol) ThrowError("weight update problem");
		if(ABS(layerHidden.m_neurons[1].m_weights[2] + 0.020112480758782) > tol) ThrowError("weight update problem");
		if(ABS(layerHidden.m_neurons[2].m_weights[0] + 0.019892890998637) > tol) ThrowError("weight update problem");
		if(ABS(layerHidden.m_neurons[2].m_weights[1] - 0.03) > tol) ThrowError("weight update problem");
		if(ABS(layerHidden.m_neurons[2].m_weights[2] - 0.019925023699046) > tol) ThrowError("weight update problem");
	}
}

void GNeuralNet_testInputGradient(GRand* pRand)
{
	for(int i = 0; i < 20; i++)
	{
		// Make the neural net
		GNeuralNet nn(pRand);
//		nn.addLayer(5);
//		nn.addLayer(10);
		GUniformRelation* pRelation = new GUniformRelation(15);
		sp_relation pRel = pRelation;
		nn.enableIncrementalLearning(pRel, 10, NULL, NULL);

		// Init with random weights
		size_t weightCount = nn.countWeights();
		double* pWeights = new double[weightCount + 5 + 10 + 10 + 5 + 5];
		ArrayHolder<double> hWeights(pWeights);
		double* pFeatures = pWeights + weightCount;
		double* pTarget = pFeatures + 5;
		double* pOutput = pTarget + 10;
		double* pFeatureGradient = pOutput + 10;
		double* pEmpiricalGradient = pFeatureGradient + 5;
		for(size_t j = 0; j < weightCount; j++)
			pWeights[j] = pRand->normal() * 0.8;
		nn.setWeights(pWeights);

		// Compute target output
		GVec::setAll(pFeatures, 0.0, 5);
		nn.predict(pFeatures, pTarget);

		// Move away from the goal and compute baseline error
		for(int i = 0; i < 5; i++)
			pFeatures[i] += pRand->normal() * 0.1;
		nn.predict(pFeatures, pOutput);
		double sseBaseline = GVec::squaredDistance(pTarget, pOutput, 10);

		// Compute the feature gradient
		nn.forwardProp(pFeatures);
		nn.setErrorOnOutputLayer(pTarget);
		nn.backProp()->backpropagate();
		GVec::copy(pFeatureGradient, pFeatures, 5);
		nn.backProp()->adjustFeatures(pFeatureGradient, 1.0);
		GVec::subtract(pFeatureGradient, pFeatures, 5);
		GVec::multiply(pFeatureGradient, -2.0, 5);

		// Empirically measure gradient
		for(int i = 0; i < 5; i++)
		{
			pFeatures[i] += 0.0001;
			nn.predict(pFeatures, pOutput);
			double sse = GVec::squaredDistance(pTarget, pOutput, 10);
			pEmpiricalGradient[i] = (sse - sseBaseline) / 0.0001;
			pFeatures[i] -= 0.0001;
		}

		// Check it
		double corr = GVec::correlation(pFeatureGradient, pEmpiricalGradient, 5);
		if(corr > 1.0)
			ThrowError("pathological results");
		if(corr < 0.999)
			ThrowError("failed");
	}
}

void GNeuralNet_testBinaryClassification(GRand* pRand)
{
	GUniformRelation* pRelation = new GUniformRelation(2, 2);
	sp_relation pRel = pRelation;
	GData dataTrain(pRel);
	for(size_t i = 0; i < 100; i++)
	{
		double* pRow = dataTrain.newRow();
		pRow[0] = (double)pRand->next(2);
		pRow[1] = 1.0 - pRow[0];
	}
	GNeuralNet nn(pRand);
	GFilter tl(&nn, false);
	tl.setFeatureTransform(new GNominalToCat(), true);
	tl.setLabelTransform(new GNominalToCat(), true);
	tl.train(&dataTrain, 1);
	double r;
	tl.accuracy(&dataTrain, &r);
	if(r != 1)
		ThrowError("Failed simple sanity test");
}

// static
void GNeuralNet::test()
{
	GRand prng(0);
	GNeuralNet_testMath();
	GNeuralNet_testBinaryClassification(&prng);

	// Test perceptron
	{
		GNeuralNet nn(&prng);
		GFilter tl(&nn, false);
		tl.setLabelTransform(new GNominalToCat(), true);
		tl.basicTest(0.76, &prng);
	}

	// Test NN with one hidden layer. (It should actually do slightly worse because
	// this test is designed to be best separated by a line, so the extra flexibility is
	// an overfit liability.)
	{
		GNeuralNet nn(&prng);
		nn.addLayer(3);
		GFilter tl(&nn, false);
		tl.setLabelTransform(new GNominalToCat(), true);
		tl.basicTest(0.76, &prng);
	}

	GNeuralNet_testInputGradient(&prng);
}

#endif





GNeuralNetPseudoInverse::GNeuralNetPseudoInverse(GNeuralNet* pNN, double padding)
: m_padding(padding)
{
	size_t maxNodes = 0;
	size_t i;
	for(i = 0; i < pNN->layerCount(); i++)
	{
		GNeuralNetLayer& nnLayer = pNN->getLayer(i);
		maxNodes = MAX(maxNodes, nnLayer.m_neurons.size());
		GNeuralNetInverseLayer* pLayer = new GNeuralNetInverseLayer();
		m_layers.push_back(pLayer);
		pLayer->m_pActivationFunction = nnLayer.m_pActivationFunction;
		GData weights(nnLayer.m_neurons[0].m_weights.size() - 1);
		weights.newRows(nnLayer.m_neurons.size());
		size_t r = 0;
		for(vector<GNeuron>::iterator it = nnLayer.m_neurons.begin(); it != nnLayer.m_neurons.end(); it++)
		{
			vector<double>::iterator itW = it->m_weights.begin();
			double unbias = -*itW;
			itW++;
			double* pRow = weights.row(r);
			for( ; itW != it->m_weights.end(); itW++)
			{
				*(pRow++) = *itW;
				unbias -= nnLayer.m_pActivationFunction->center() * (*itW);
			}
			pLayer->m_unbias.push_back(unbias);
			r++;
		}
		pLayer->m_pInverseWeights = weights.pseudoInverse();
	}
	m_pBuf1 = new double[2 * maxNodes];
	m_pBuf2 = m_pBuf1 + maxNodes;
}

GNeuralNetPseudoInverse::~GNeuralNetPseudoInverse()
{
	for(vector<GNeuralNetInverseLayer*>::iterator it = m_layers.begin(); it != m_layers.end(); it++)
		delete(*it);
	delete[] MIN(m_pBuf1, m_pBuf2);
}

void GNeuralNetPseudoInverse::computeFeatures(const double* pLabels, double* pFeatures)
{
	size_t inCount = 0;
	vector<GNeuralNetInverseLayer*>::iterator it = m_layers.end() - 1;
	GVec::copy(m_pBuf2, pLabels, (*it)->m_pInverseWeights->cols());
	for(; true; it--)
	{
		GNeuralNetInverseLayer* pLayer = *it;
		inCount = pLayer->m_pInverseWeights->rows();
		std::swap(m_pBuf1, m_pBuf2);

		// Invert the layer
		double* pT = m_pBuf1;
		for(vector<double>::iterator ub = pLayer->m_unbias.begin(); ub != pLayer->m_unbias.end(); ub++)
		{
			*pT = pLayer->m_pActivationFunction->inverse(*pT) + *ub;
			pT++;
		}
		pLayer->m_pInverseWeights->multiply(m_pBuf1, m_pBuf2);

		// Clip and uncenter the value
		pLayer = *it;
		double halfRange = pLayer->m_pActivationFunction->halfRange();
		double center = pLayer->m_pActivationFunction->center();
		pT = m_pBuf2;
		for(size_t i = 0; i < inCount; i++)
		{
			*pT = MAX(m_padding - halfRange, MIN(halfRange - m_padding, *pT)) + center;
			pT++;
		}

		if(it == m_layers.begin())
			break;
	}
	GVec::copy(pFeatures, m_pBuf2, inCount);
}

#ifndef NO_TEST_CODE
// static
void GNeuralNetPseudoInverse::test()
{
	GRand prng(0);
	GNeuralNet nn(&prng);
	nn.addLayer(5);
	nn.addLayer(7);
	sp_relation pRel = new GUniformRelation(15);
	nn.enableIncrementalLearning(pRel, 12, NULL, NULL);
	nn.decayWeights(-9.0 * nn.learningRate()); // multiply all non-bias weights by 10
	GNeuralNetPseudoInverse nni(&nn, 0.001);
	double labels[12];
	double features[3];
	double features2[3];
	for(size_t i = 0; i < 20; i++)
	{
		for(size_t j = 0; j < 3; j++)
			features[j] = prng.uniform() * 0.98 + 0.01;
		nn.predict(features, labels);
		nni.computeFeatures(labels, features2);
		if(GVec::squaredDistance(features, features2, 3) > 1e-8)
			ThrowError("failed");
	}
}
#endif















/*
GModerateNet::GModerateNet(GRand* pRand)
: GNeuralNet(pRand)
{
	m_lambda = 0.01;
}

// virtual
GModerateNet::~GModerateNet()
{
}

// virtual
GTwtNode* GModerateNet::toTwt(GTwtDoc* pDoc)
{
	ThrowError("Not implemented yet");
	return NULL;
}

// virtual
void GModerateNet::train(GData* pData, int labelDims)
{
	// Split into a training and test set
	GData dataTrain(pData->relation());
	GReleaseDataHolder hDataTrain(&dataTrain);
	for(size_t i = 0; i < pData->rows(); i++)
		dataTrain.takeRow(pData->row(i));
	GData dataValidate(pData->relation());
	GReleaseDataHolder hDataValidate(&dataValidate);
//	dataTrain.splitBySize(&dataValidate, dataTrain.rows() / 2);

	// Make an internal features set
	GData* pInternalSet = dataTrain.clone();
	Holder<GData> hInternalSet(pInternalSet);

	// Do the training
	GNeuralNet::enableIncrementalLearning(pData->relation(), labelDims, NULL, NULL);
	size_t featureDims = pData->cols() - labelDims;
	double dBestError = 1e308;
	int nEpochsSinceValidationCheck = 0;
	double dSumSquaredError;
	while(true)
	{
		pInternalSet->shuffle2(m_pRand, dataTrain);
		for(size_t i = 0; i < pInternalSet->rows(); i++)
		{
			// Train the weights
			double* pRow = pInternalSet->row(i);
			forwardProp(pRow);
			setErrorOnOutputLayer(pRow + featureDims);
			m_pBackProp->backpropagate();
			m_pBackProp->descendGradient(pRow, m_learningRate, m_momentum);

			// Train the inputs
			m_pBackProp->adjustFeatures(pRow, m_learningRate);

			// Decay the inputs
			GVec::multiply(pRow, 1.0 - m_lambda, featureDims);
			GVec::addScaled(pRow, m_lambda, dataTrain.row(i), featureDims);
		}

		// Check for termination condition
		if(nEpochsSinceValidationCheck >= m_epochsPerValidationCheck)
		{
			nEpochsSinceValidationCheck = 0;
			dSumSquaredError = validationSquaredError(&dataTrain);
			if(1.0 - dSumSquaredError / dBestError < m_minImprovement)
				break;
			if(dSumSquaredError < dBestError)
				dBestError = dSumSquaredError;
		}
		else
			nEpochsSinceValidationCheck++;
	}

	releaseTrainingJunk();
}

// virtual
void GModerateNet::enableIncrementalLearning(sp_relation& pRelation, int labelDims, double* pMins, double* pRanges)
{
	ThrowError("Not implemented yet");
}

// virtual
void GModerateNet::predictDistribution(const double* pIn, GPrediction* pOut)
{
	ThrowError("Not implemented yet");
}

// virtual
void GModerateNet::clear()
{
	GNeuralNet::clear();
}

// virtual
void GModerateNet::trainIncremental(const double* pIn, const double* pOut)
{
	ThrowError("Not implemented yet");
}

// virtual
void GModerateNet::trainSparse(GSparseMatrix* pData, int labelDims)
{
	ThrowError("Not implemented yet");
}
*/






} // namespace GClasses

