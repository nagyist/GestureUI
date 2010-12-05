/*
	Copyright (C) 2006, Mike Gashler

	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.

	see http://www.gnu.org/copyleft/lesser.html
*/

#ifndef __GPOLYNOMIAL_H__
#define __GPOLYNOMIAL_H__

#include "GLearner.h"

namespace GClasses {

/// This regresses a multi-dimensional polynomial to fit the data
class GPolynomial : public GSupervisedLearner
{
protected:
	int m_featureDims;
	int m_nControlPoints;
	int m_nCoefficients;
	double* m_pCoefficients;

public:
	/// It will have the same number of control points in every dimension
	/// This class assumes there is only one output value.
	GPolynomial(int nControlPoints);

	/// Load from a text-based format
	GPolynomial(GTwtNode* pNode);

	virtual ~GPolynomial();

#ifndef NO_TEST_CODE
	/// Performs unit tests for this class. Throws an exception if there is a failure.
	static void test();
#endif // NO_TEST_CODE

	/// Save to a text-based format
	virtual GTwtNode* toTwt(GTwtDoc* pDoc);

	/// See the comment for GSupervisedLearner::featureDims
	virtual int featureDims();

	/// See the comment for GSupervisedLearner::labelDims
	virtual int labelDims() { return 1; }

	/// Specify the number of input and output features
	void init(GRelation* pRelation);

	/// Returns the total number of coefficients in this polynomial
	int coefficientCount() { return m_nCoefficients; }

	/// Returns the number of control points (per dimension)
	int controlPointCount() { return m_nControlPoints; }

	/// Returns the coefficient at the specified coordinates. pCoords should
	/// be an array of size m_nDimensions, and each value should be from 0 to m_nControlPoints - 1
	double coefficient(int* pCoords);

	/// Returns the full array of coefficients
	double* coefficientArray() { return m_pCoefficients; }

	/// Sets the coefficient at the specified coordinates. pCoords should
	/// be an array of size m_nDimensions, and each value should be from 0 to m_nControlPoints - 1
	void setCoefficient(int* pCoords, double dVal);

	/// Copies pOther into this polynomial. Both polynomials must have the
	/// same dimensionality, and this polynomial must have >= 
	void copy(GPolynomial* pOther);

	/// Regress the polynomial using all the samples in a collection
	virtual void train(GData* pData, int labelDims);

	/// See the comment for GSupervisedLearner::predictDistribution
	virtual void predictDistribution(const double* pIn, GPrediction* pOut);

	/// See the comment for GSupervisedLearner::predict
	virtual void predict(const double* pIn, double* pOut);

	/// See the comment for GSupervisedLearner::clear
	virtual void clear();

	/// Sets all the coefficients. pVector must be of size GetCoefficientCount()
	void setCoefficients(const double* pVector);

	/// Converts to a multi-dimensional Bezier curve
	void toBezierCoefficients();

	/// Converts from a multi-dimensional Bezier curve
	void fromBezierCoefficients();

	/// Differentiates the polynomial with respect to every dimension
	void differentiate();

	/// Integrates the polynomial in every dimension. This assumes the
	/// constant of integration is always zero. It also assumes that all
	/// of the highest-order coefficients are zero. If that isn't true,
	/// this polynomial won't be big enough to hold the answer, and the
	/// highest-order coefficients will be dropped. The best way to ensure
	/// that doesn't happen is to copy into a bigger (one more control point)
	/// polynomial before integrating.
	void integrate();

protected:
	/// This converts from control-point-lattice coordinates to an array index.
	/// The array is stored with lattice position (0, 0, 0, ...) (the constant coefficient)
	/// in array position 0, and lattice position (1, 0, 0, ...) in position 1, etc.
	int calcIndex(int* pDegrees);
};


} // namespace GClasses

#endif // __GPOLYNOMIAL_H__
