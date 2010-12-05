/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#ifndef __GNEURALNET_H__
#define __GNEURALNET_H__

#include "GLearner.h"
#include <vector>

namespace GClasses {

class GNeuralNet;
class GRand;
class GBackProp;
class GImage;
class GActivationFunction;



/// Represents a single neuron in a neural network
class GNeuron
{
public:
	double m_value;
	std::vector<double> m_weights; // Weight zero is always the bias weight

	void resetWeights(GRand* pRand, double inputCenter);
};

/// Represents a layer of neurons in a neural network
class GNeuralNetLayer
{
public:
	std::vector<GNeuron> m_neurons;
	GActivationFunction* m_pActivationFunction;

	void resetWeights(GRand* pRand, double inputCenter);
};

/// An internal class used by GBackProp
class GBackPropWeight
{
public:
	double m_delta;

	GBackPropWeight()
	{
		m_delta = 0;
	}
};

/// An internal class used by GBackProp
class GBackPropNeuron
{
public:
	double m_error;
	std::vector<GBackPropWeight> m_weights;
};

/// An internal class used by GBackProp
class GBackPropLayer
{
public:
	std::vector<GBackPropNeuron> m_neurons;
};

/// This class performs backpropagation on a neural network. (It is a separate
/// class, because it is only needed while training. There is no reason to waste
/// this space after training is complete, or if you choose to use a different
/// technique to train the neural network.)
class GBackProp
{
friend class GNeuralNet;
protected:
	GNeuralNet* m_pNN;
	std::vector<GBackPropLayer> m_layers;

public:
	/// This class will adjust the weights in pNN
	GBackProp(GNeuralNet* pNN);

	~GBackProp()
	{
	}

	/// Returns a layer (not a layer of the neural network, but a corresponding layer of values used for back-prop)
	GBackPropLayer& layer(int layer)
	{
		return m_layers[layer];
	}

	/// Backpropagates the error from the "from" layer to the "to" layer. (If the "to" layer has fewer units than the "from"
	/// layer, then it will begin propagating with the (fromBegin+1)th weight and stop when the "to" layer runs out of units.
	/// It would be an error if the number of units in the "from" layer is less than the number of units in the "to" layer
	/// plus fromBegin.
	static void backPropLayer(GNeuralNetLayer* pNNFromLayer, GNeuralNetLayer* pNNToLayer, GBackPropLayer* pBPFromLayer, GBackPropLayer* pBPToLayer, size_t fromBegin = 0);

	/// This is another implementation of backPropLayer. This one is somewhat more flexible, but slightly less efficient.
	/// It supports backpropagating error from one or two layers. (pNNFromLayer2 should be NULL if you are backpropagating from just one
	/// layer.) It also supports temporal backpropagation by unfolding in time and then averaging the error across all of the unfolded
	/// instantiations. "pass" specifies how much of the error for this pass to accept. 1=all of it, 2=half of it, 3=one third, etc.
	static void backPropLayer2(GNeuralNetLayer* pNNFromLayer1, GNeuralNetLayer* pNNFromLayer2, GNeuralNetLayer* pNNToLayer, GBackPropLayer* pBPFromLayer1, GBackPropLayer* pBPFromLayer2, GBackPropLayer* pBPToLayer, int pass);

	/// Adjust weights in pNNFromLayer. (The error for pNNFromLayer layer must have already been computed.) (If you are
	/// backpropagating error from two layers, you can just call this method twice, once for each previous layer.)
	static void adjustWeights(GNeuralNetLayer* pNNFromLayer, GNeuralNetLayer* pNNToLayer, GBackPropLayer* pBPFromLayer, double learningRate, double momentum);

	/// Adjust weights in pNNFromLayer. (The error for pNNFromLayer layer must have already been computed.) (If you are
	/// backpropagating error from two layers, you can just call this method twice, once for each previous layer.)
	static void adjustWeights(GNeuralNetLayer* pNNFromLayer, const double* pFeatures, GBackPropLayer* pBPFromLayer, double learningRate, double momentum);

	/// This method assumes that the error term is already set at every unit in the output layer. It uses back-propagation
	/// to compute the error term at every hidden unit. (It does not update any weights.)
	void backpropagate();

	/// This method assumes that the error term is alrady set for every network unit. It adjusts weights to descend the
	/// gradient of the error surface with respect to the weights.
	void descendGradient(const double* pFeatures, double learningRate, double momentum);

	/// This method assumes that the error term is alrady set for every network unit. It descends the gradient
	/// by adjusting the features (not the weights).
	void adjustFeatures(double* pFeatures, double learningRate, size_t skip = 0);
};



/// An artificial neural network
class GNeuralNet : public GIncrementalLearner
{
friend class GBackProp;
public:
	enum TargetFunction
	{
		squared_error, /// (default) best for regression
		cross_entropy, /// best for classification
		root_cubed_error, /// slow, but extra stable
	};

protected:
	std::vector<GNeuralNetLayer> m_layers;
	GRand* m_pRand;
	GBackProp* m_pBackProp;
	int m_featureDims, m_labelDims;
	std::vector<GActivationFunction*> m_activationFunctions;
	GActivationFunction* m_pActivationFunction;
	double m_learningRate;
	double m_momentum;
	double m_validationPortion;
	double m_minImprovement;
	int m_epochsPerValidationCheck;
	TargetFunction m_backPropTargetFunction;

public:
	GNeuralNet(GRand* pRand);

	/// Load from a text-format
	GNeuralNet(GTwtNode* pNode, GRand* pRand);

	virtual ~GNeuralNet();

#ifndef NO_TEST_CODE
	/// Performs unit tests for this class. Throws an exception if there is a failure.
	static void test();
#endif

	/// Saves the model to a text file. (This doesn't save the short-term
	/// memory used for incremental learning, so if you're doing "incremental"
	/// learning, it will wake up with amnesia when you load it again.)
	virtual GTwtNode* toTwt(GTwtDoc* pDoc);

	/// See the comment for GSupervisedLearner::featureDims
	virtual int featureDims();

	/// See the comment for GSupervisedLearner::labelDims
	virtual int labelDims();

	/// Sets the activation function to use with all subsequently added
	/// layers. (Note that the activation function for the output layer is
	/// set when train or enableIncrementalLearning is called, so if you
	/// only wish to set the squshing function for the output layer, call
	/// this method after all hidden layers have been added, but before you call train.)
	/// If hold is true, then the neural network will hold on to this instance
	/// of the activation function and delete it when the neural network is deleted.
	void setActivationFunction(GActivationFunction* pSF, bool hold);

	/// Adds a hidden layer to the network. (The first hidden layer
	/// that you add will be adjacent to the input features. The last
	/// hidden layer that you add will be adjacent to the output
	/// layer.)
	void addLayer(int nNodes);

	/// Returns the number of layers in this neural network. (Every network has
	/// at least one output layer, plus all of the hidden layers that you add by calling
	/// addLayer.)
	size_t layerCount() { return m_layers.size(); }

	/// Returns a reference to the specified layer.
	GNeuralNetLayer& getLayer(int layer) { return m_layers[layer]; }

	/// Adds a new node at the end of the specified layer. (The new node is initialized
	/// with small weights, so this operation should initially have little impact on
	/// predictions.)
	void addNode(size_t layer);

	/// Removes the specified node from the specified layer. (An exception will be thrown
	/// the layer only has one node.)
	void dropNode(size_t layer, size_t node);

	/// Returns the backprop object associated with this neural net (if there is one)
	GBackProp* backProp() { return m_pBackProp; }

	/// Set the portion of the data that will be used for validation. If the
	/// value is 0, then all of the data is used for both training and validation.
	void setValidationPortion(double d) { m_validationPortion = d; }

	/// Counts the number of weights in the network. (This value is not cached, so
	/// you should cache it rather than frequently call this method.)
	size_t countWeights();

	/// Perturbs all weights in the network by a random normal offset with the
	/// specified deviation.
	void perturbAllWeights(double deviation);

	/// Clips all non-bias weights to fall within the range [-max, max].
	void clipWeights(double max);

	/// Multiplies all non-bias weights by (1.0 - (learning_rate * lambda)),
	/// starting with the output layer, and ending with the first hidden layer.
	/// Typical values for lambda are small (like 0.001.)
	/// After each layer, the value of lambda is multiplied by gamma.
	/// (If gamma is greater than 1.0, then weights in hidden layers will decay
	/// faster, and if gamma is less than 1.0, then weights in hidden layers will
	/// decay slower.) It may be significant to note that if a regularizing
	/// penalty is added to the error of lambda times the sum-squared values of
	/// non-bias weights, then on-line weight updating works out to the same as
	/// decaying the weights after each application of back-prop.
	void decayWeights(double lambda, double gamma = 1.0);

	/// Returns the current learning rate
	double learningRate() { return m_learningRate; }

	/// Set the learning rate
	void setLearningRate(double d) { m_learningRate = d; }

	/// Performs a single epoch of training
	void trainEpoch(GData* pTrainingData);

	/// Returns the current momentum value
	double momentum() { return m_momentum; }

	/// Momentum has the effect of speeding convergence and helping
	/// the gradient descent algorithm move past some local minimums
	void setMomentum(double d) { m_momentum = d; }

	/// Specifies the minimum improvement (as a ratio) that must be
	/// made since the last validation check for trainingn to continue.
	/// (For example, if the mean squared error at the previous validation check
	/// was 50, and the mean squared error at the current validation check
	/// is 49, then training will stop if d is > 0.02.)
	void setMinImprovement(double d) { m_minImprovement = d; }

	/// Sets the number of iterations that will be performed before
	/// each time the network is tested again with the validation set
	/// to determine if we have a better best-set of weights, and
	/// whether or not it's achieved the termination condition yet.
	/// (An iteration is defined as a single pass through all rows in
	/// the training set.)
	void setIterationsPerValidationCheck(int n) { m_epochsPerValidationCheck = n; }

	/// Specify the target function to use for back-propagation. The default is squared_error.
	/// cross_entropy tends to be faster, and is well-suited for classification tasks.
	void setBackPropTargetFunction(TargetFunction eTF) { m_backPropTargetFunction = eTF; }

	/// Splits the provided data into a training and validation set and trains
	/// the network. To set the ratio, use SetTrainingPortion.
	virtual void train(GData* pData, int labelDims);

	/// See the comment for GIncrementalLearner::trainSparse
	/// Assumes all attributes are continuous.
	virtual void trainSparse(GSparseMatrix* pData, int labelDims);

	/// See the comment for GIncrementalLearner::enableIncrementalLearning
	virtual void enableIncrementalLearning(sp_relation& pRelation, int labelDims, double* pMins, double* pRanges);

	/// See the comment for GIncrementalLearner::trainIncremental
	virtual void trainIncremental(const double* pIn, const double* pOut);

	/// See the comment for GSupervisedLearner::predictDistribution
	virtual void predictDistribution(const double* pIn, GPrediction* pOut);

	/// See the comment for GSupervisedLearner::predict
	virtual void predict(const double* pIn, double* pOut);

	/// See the comment for GSupervisedLearner::clear
	virtual void clear() {}

	/// Train the network until the termination condition is met.
	/// Returns the number of epochs required to train it.
	int trainWithValidation(GData* pTrainingData, GData* pValidationData, int labelDims);

	/// Some extra junk is allocated when training to make it efficient.
	/// This method is called when training is done to get rid of that
	/// extra junk.
	void releaseTrainingJunk();

	/// Gets the internal training data set
	GData* internalTrainingData();

	/// Gets the internal validation data set
	GData* internalValidationData();

	/// Sets all the weights from an array of doubles. The number of
	/// doubles in the array can be determined by calling countWeights().
	void setWeights(const double* pWeights);

	/// Copy the weights from pOther. It is assumed (but not checked) that
	/// pOther has the same network structure as this neural network.
	void copyWeights(GNeuralNet* pOther);

	/// Copies the layers, nodes, and settings from pOther (but not the
	/// weights). enableIncrementalLearning must have been called on pOther
	/// so that it has a complete structure.
	void copyStructure(GNeuralNet* pOther);

	/// Serializes the network weights into an array of doubles. The
	/// number of doubles in the array can be determined by calling
	/// countWeights().
	void weights(double* pOutWeights);

	/// Convert all the input values to the internal representation
	void inputsToInternal(const double* pExternal, double* pInternal);

	/// Convert all the output values to the internal representation
	/// pExternal points to the first external output value (so just the
	/// output part of the row).
	/// pInternal points to the first internal input value (so it's the
	/// whole internal row).
	void outputsToInternal(const double* pExternal, double* pInternal);

	/// Convert the internal output values to the external representation
	void outputsToExternal(GPrediction* pExternal);

	/// Returns the pseudo-random number generator associated with this object
	GRand* getRand() { return m_pRand; }

	/// Evaluates a feature vector. (The results will be in the nodes of the output layer.)
	void forwardProp(const double* pRow);

	/// This method assumes forwardProp has been called. It copies the predicted vector into pOut.
	void copyPrediction(double* pOut);

	/// This method assumes forwardProp has been called. It computes the sum squared prediction error
	/// with the specified target vector.
	double sumSquaredPredictionError(const double* pTarget);

	/// This method assumes that forwardProp has already been called. (Note that
	/// the predict method calls forwardProp). It computes the error
	/// values at each node in the output layer. After calling this method,
	/// it is typical to call backProp()->backpropagate(), to compute the error on
	/// the hidden nodes, and then to call backProp()->descendGradient to update
	/// the weights. pTarget contains the target values for the ouptut nodes.
	void setErrorOnOutputLayer(const double* pTarget, TargetFunction eTargetFunction = squared_error);

protected:
	/// Measures the sum squared error against the specified dataset
	double validationSquaredError(GData* pValidationData);
};


/// A helper class used by GNeuralNetPseudoInverse
class GNeuralNetInverseLayer
{
public:
	GActivationFunction* m_pActivationFunction;
	std::vector<double> m_unbias;
	GData* m_pInverseWeights;

	~GNeuralNetInverseLayer()
	{
		delete(m_pInverseWeights);
	}
};

/// Computes the pseudo-inverse of a neural network.
class GNeuralNetPseudoInverse
{
protected:
	double m_padding;
	std::vector<GNeuralNetInverseLayer*> m_layers;
	double* m_pBuf1;
	double* m_pBuf2;

public:
	/// padding specifies a margin in which label values will be clipped inside
	/// the activation function output range to avoid extreme feature values (-inf, inf, etc.).
	GNeuralNetPseudoInverse(GNeuralNet* pNN, double padding = 0.01);
	~GNeuralNetPseudoInverse();

	/// Computes the input features from the output labels. In cases of
	/// under-constraint, the feature vector with the minimum magnitude is chosen.
	/// In cases of over-constraint, the feature vector is chosen with a corresponding
	/// label vector that minimizes sum-squared error with the specified label
	/// vector.
	void computeFeatures(const double* pLabels, double* pFeatures);

#ifndef NO_TEST_CODE
	static void test();
#endif
};

/*
/// This is an experimental neural network that has the ability to adjust features (inputs) as well as weights
/// in order to make good predictions. The idea is that this should prevent outliers from having too much
/// influence on the model. The value of this idea, however, is not yet well-established.
class GModerateNet : public GNeuralNet
{
protected:
	double m_lambda;

public:
	GModerateNet(GRand* pRand);
	virtual ~GModerateNet();
	double lambda() { return m_lambda; }
	void setLambda(double d) { m_lambda = d; }
	virtual GTwtNode* toTwt(GTwtDoc* pDoc);
	virtual void train(GData* pData, int labelDims);
	virtual void predictDistribution(const double* pIn, GPrediction* pOut);
	virtual void clear();
	virtual void enableIncrementalLearning(sp_relation& pRelation, int labelDims, double* pMins, double* pRanges);
	virtual void trainIncremental(const double* pIn, const double* pOut);
	virtual void trainSparse(GSparseMatrix* pData, int labelDims);
};
*/
} // namespace GClasses

#endif // __GNEURALNET_H__

