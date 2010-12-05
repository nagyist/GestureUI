/*
	Copyright (C) 2010, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#ifndef __GACTIVATION_H__
#define __GACTIVATION_H__

#include <math.h>
#include "GMacros.h"

namespace GClasses {

class GTwtNode;
class GTwtDoc;


/// The base class for activation functions. Typically, this are
/// sigmoid-shaped functions used to "squash" the output of a network
/// node. These are typically used in conjunction with the GNeuralNet class.
class GActivationFunction
{
public:
	/// Returns the name of this activation function
	virtual const char* name() = 0;

	/// The activation function
	virtual double squash(double x) = 0;

	/// The derivative of the activation function
	virtual double derivative(double x) = 0;

	/// The inverse of the activation function
	virtual double inverse(double y) = 0;

	/// The center output value. This should return the value of squash(0.0).
	virtual double center() = 0;

	/// The absolute difference between the max (or min) output value and the center
	virtual double halfRange() = 0;

	/// Returns a clone of this object
	virtual GActivationFunction* clone() = 0;

	/// The derivative function applied to the result of the inverse function.
	/// This method exists for efficiency reasons, because in some cases, such
	/// as GActivationLogistic, the combination of these two methods can be simplified
	/// algebraically, thus reducing computational overhead.
	virtual double derivativeInverse(double y) { return derivative(inverse(y)); }

	/// Serialize this object
	GTwtNode* toTwt(GTwtDoc* pDoc);

	/// Deserialize this object
	static GActivationFunction* fromTwt(GTwtNode* pNode);
};

/// The logistic activation function
class GActivationLogistic : public GActivationFunction
{
public:
	/// Returns the name of this activation function
	virtual const char* name() { return "logistic"; }

	/// The logistic function. Returns 1.0/(e^(-x)+1.0)
	virtual double squash(double x)
	{
		if(x >= 700.0) // Don't trigger a floating point exception
			return 1.0;
		if(x < -700.0) // Don't trigger a floating point exception
			return 0.0;
		return 1.0 / (exp(-x) + 1.0);
	}

	/// Returns d*(1.0-d), where d=squash(x)
	virtual double derivative(double x) { double d = squash(x); return d * (1.0 - d); }

	/// The logit function. Returns log(y)-log(1.0-y)
	virtual double inverse(double y)
	{
		if(y >= 1.0)
			return 700.0;
		if(y <= 0.0)
			return -700.0;
		return (log(y) - log(1.0 - y));
	}

	/// Returns y*(1.0-y)
	virtual double derivativeInverse(double y) { return y * (1.0 - y); }

	/// Returns 0.5
	virtual double center() { return 0.5; }

	/// Returns 0.5
	virtual double halfRange() { return 0.5; }

	/// See the comment for GActivationFunction::clone
	virtual GActivationFunction* clone() { return new GActivationLogistic(); }
};

/// The arctan activation function
class GActivationArcTan : public GActivationFunction
{
public:
	/// Returns the name of this activation function
	virtual const char* name() { return "arctan"; }

	/// Returns atan(x). The result will be in the range -PI/2 <= y <= PI/2
	virtual double squash(double x) { return atan(x); }

	/// Returns 1/(x*x+1.0)
	virtual double derivative(double x) { return 1.0 / (x * x + 1.0); }

	/// Returns tan(y), where -PI/2 <= y <= PI/2
	virtual double inverse(double y) { return tan(y); }

	/// Returns 0
	virtual double center() { return 0.0; }

	/// Returns PI / 2
	virtual double halfRange();

	/// See the comment for GActivationFunction::clone
	virtual GActivationFunction* clone() { return new GActivationArcTan(); }
};

/// The hyperbolic tangent activation function
class GActivationTanH : public GActivationFunction
{
public:
	/// Returns the name of this activation function
	virtual const char* name() { return "tanh"; }

	/// Returns tanh(x). The result is in the range -1 <= y <= 1
	virtual double squash(double x) { return tanh(x); }

	/// Returns sech(x)*sech(x)
	virtual double derivative(double x)
	{
		double d = 2.0 / (exp(x) + exp(-x)); // sech(x)
		return d * d;
	}

	/// Returns atanh(y), where -1 <= y <= 1
	virtual double inverse(double y)
	{
#ifdef WIN32
		return 0.5 * log((1.0 + y) / (1.0 - y));
#else
		return atanh(y);
#endif
	}

	/// Returns 0.0
	virtual double center() { return 0.0; }

	/// Returns 1.0
	virtual double halfRange() { return 1.0; }

	/// See the comment for GActivationFunction::clone
	virtual GActivationFunction* clone() { return new GActivationTanH(); }
};

/// The hyperbolic tangent activation function
class GActivationAlgebraic : public GActivationFunction
{
public:
	/// Returns the name of this activation function
	virtual const char* name() { return "algebraic"; }

	/// Returns x/(sqrt(x*x+1.0). The result is in the range -1 <= y <= 1
	virtual double squash(double x) { return x / (sqrt(x * x + 1.0)); }

	/// Returns 1.0/(sqrt(x*x+1))-(x*x)/pow(x*x+1,1.5)
	virtual double derivative(double x) { return 1.0 / (sqrt(x * x + 1)) - (x * x) / pow(x * x + 1, 1.5); }

	/// Returns y / (sqrt(1.0 - (y * y)))
	virtual double inverse(double y) { return y / (sqrt(1.0 - (y * y))); }

	/// Returns 0.0
	virtual double center() { return 0.0; }

	/// Returns 1.0
	virtual double halfRange() { return 1.0; }

	/// See the comment for GActivationFunction::clone
	virtual GActivationFunction* clone() { return new GActivationAlgebraic(); }
};

/// Use this function when you do not want to squash the net. For example,
/// using this activation function with a network that has no hidden layers
/// makes a perceptron model. Also, it is common to use this activation
/// function on the output layer for regression problems.
class GActivationIdentity : public GActivationFunction
{
public:
	/// Returns the name of this activation function
	virtual const char* name() { return "identity"; }

	/// Returns x
	virtual double squash(double x) { return x; }

	/// Returns 1.0
	virtual double derivative(double x) { return 1.0; }

	/// Returns y
	virtual double inverse(double y) { return y; }

	/// Returns 1.0
	virtual double derivativeInverse(double y) { return 1.0; }

	/// Returns 0.0
	virtual double center() { return 0.0; }

	/// Returns 1e308
	virtual double halfRange() { return 1e308; }

	/// See the comment for GActivationFunction::clone
	virtual GActivationFunction* clone() { return new GActivationIdentity(); }
};

/// This provides an alternative to using GActivationIdentity on the output layer
/// for regression problems. It may add more power because it is non-linear, but
/// like the identity function, its co-domain is the same as its domain.
class GActivationBend : public GActivationFunction
{
public:
	/// Returns the name of this activation function
	virtual const char* name() { return "bend"; }

	/// Returns x+log(exp(x)+1.0)
	virtual double squash(double x)
	{
		if(x >= 700.0) // Don't trigger a floating point exception
			return x + x;
		if(x < -700.0) // Don't trigger a floating point exception
			return x;
		return x + log(exp(x) + 1.0);
	}

	/// Returns the logistic function of x, plus 1
	virtual double derivative(double x)
	{
		if(x >= 700.0) // Don't trigger a floating point exception
			return 2.0;
		if(x < -700.0) // Don't trigger a floating point exception
			return 1.0;
		return 1.0 / (exp(-x) + 1.0) + 1.0;
	}

	/// Returns log(0.5*(sqrt(4.0*exp(y)+1.0)-1.0))
	virtual double inverse(double y)
	{
		if(y >= 1000.0)
			return 0.5 * y;
		if(y < -500.0)
			return y;
		return log(0.5 * (sqrt(4.0 * exp(y) + 1.0) - 1.0));
	}

	/// Returns 0.0
	virtual double center() { return 0.0; }

	/// Returns 1e308
	virtual double halfRange() { return 1e308; }

	/// See the comment for GActivationFunction::clone
	virtual GActivationFunction* clone() { return new GActivationBend(); }
};

/// This is an output-layer activation function shaped
/// like a sigmoid, but with both a co-domain and domain
/// that spans the continuous values.
class GActivationBiDir : public GActivationFunction
{
public:
	/// Returns the name of this activation function
	virtual const char* name() { return "bidir"; }

	virtual double squash(double x)
	{
		double d = sqrt(x * x + 1.0);
		return sqrt(d + x) - sqrt(d - x);
	}

	virtual double derivative(double x)
	{
		if(ABS(x) > 1e7)
			return 0.0;
		double d = sqrt(x * x + 1.0);
		return (x / d + 1.0) / (2.0 * sqrt(d + x)) - (x / d - 1.0) / (2.0 * sqrt(d - x));
	}

	virtual double inverse(double y)
	{
		double d = y * y;
		if(y >= 0.0)
			return 0.5 * sqrt(d * d + 4.0 * d);
		else
			return -0.5 * sqrt(d * d + 4.0 * d);
	}

	/// Returns 0.0
	virtual double center() { return 0.0; }

	/// Returns 1e308
	virtual double halfRange() { return 1e308; }

	/// See the comment for GActivationFunction::clone
	virtual GActivationFunction* clone() { return new GActivationBiDir(); }
};

/// This is an experimental activation function intended to
/// reduce the required computation involved in inverting neural networks.
class GActivationPiecewise : public GActivationFunction
{
public:
	/// Returns the name of this activation function
	virtual const char* name() { return "piecewise"; }

	virtual double squash(double x);

	virtual double derivative(double x)
	{
		ThrowError("Not implemented yet");
		return 0;
	}

	virtual double inverse(double y)
	{
		ThrowError("Not implemented yet");
		return 0;
	}

	/// Returns y*(1.0-y)
	virtual double derivativeInverse(double y) { return y * (1.0 - y); }

	/// Returns 0.5
	virtual double center() { return 0.5; }

	/// Returns 0.5
	virtual double halfRange() { return 0.5; }

	/// See the comment for GActivationFunction::clone
	virtual GActivationFunction* clone() { return new GActivationPiecewise(); }
};


} // namespace GClasses

#endif // __GACTIVATION_H__

